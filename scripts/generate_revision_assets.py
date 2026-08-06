#!/usr/bin/env python3
"""Generate revision statistics and figures from the deposited LAMMPS data.

The script is intentionally self-contained so that the manuscript figures and the
reported numerical values can be regenerated directly from the repository.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"
ANALYSIS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

DT_PS = 0.0005  # 0.5 fs in ps


def parse_lammps_thermo(path: Path) -> dict[str, np.ndarray]:
    """Parse all custom thermo blocks echoed in a LAMMPS log."""
    records: list[dict[str, float]] = []
    columns: list[str] | None = None
    aliases = {
        "Step": "step",
        "Temp": "temp",
        "Press": "press",
        "Density": "density",
        "Volume": "vol",
        "Vol": "vol",
        "PotEng": "pe",
        "KinEng": "ke",
        "TotEng": "etotal",
        "Lx": "lx",
        "Ly": "ly",
        "Lz": "lz",
    }

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped:
                continue
            tokens = stripped.split()
            if tokens and tokens[0] == "Step" and "Temp" in tokens:
                columns = [aliases.get(tok, tok.lower()) for tok in tokens]
                continue
            if columns is None:
                continue
            if stripped.startswith("Loop time") or stripped.startswith("ERROR"):
                columns = None
                continue
            try:
                values = [float(tok) for tok in tokens]
            except ValueError:
                continue
            if len(values) != len(columns):
                continue
            records.append(dict(zip(columns, values)))

    if not records:
        raise RuntimeError(f"No thermo records found in {path}")

    # Keep the last occurrence of duplicated step values at run boundaries.
    by_step: dict[int, dict[str, float]] = {}
    for record in records:
        by_step[int(record["step"])] = record
    ordered = [by_step[key] for key in sorted(by_step)]
    keys = sorted({key for record in ordered for key in record})
    return {key: np.asarray([record.get(key, np.nan) for record in ordered]) for key in keys}


def load_two_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1]


def downsample_indices(length: int, maximum: int = 450) -> np.ndarray:
    if length <= maximum:
        return np.arange(length)
    return np.unique(np.linspace(0, length - 1, maximum, dtype=int))


def sample_mean_sd(values: np.ndarray) -> tuple[float, float]:
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    return float(np.mean(clean)), float(np.std(clean, ddof=1))


def linear_fit(time_ps: np.ndarray, msd: np.ndarray, start: float, end: float) -> dict[str, float]:
    mask = (time_ps >= start) & (time_ps <= end)
    x = time_ps[mask]
    y = msd[mask]
    slope, intercept = np.polyfit(x, y, 1)
    prediction = slope * x + intercept
    ss_res = np.sum((y - prediction) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    diffusion = slope / 6.0 * 1.0e-4
    return {
        "start_ps": start,
        "end_ps": end,
        "slope_A2_per_ps": float(slope),
        "intercept_A2": float(intercept),
        "D_cm2_per_s": float(diffusion),
        "R2": float(r2),
    }


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    short = parse_lammps_thermo(ROOT / "log.nvt_short.lammps")
    npt = parse_lammps_thermo(ROOT / "log.npt_eq.lammps")
    direct = parse_lammps_thermo(ROOT / "log.nvt_direct_controlled.lammps")
    prod = parse_lammps_thermo(ROOT / "log.nvt_prod_long.lammps")

    msd_step, msd_water = load_two_columns(ROOT / "msd_waterO.dat")
    msd_caf_step, msd_caf = load_two_columns(ROOT / "msd_caffeine.dat")
    msd_time = msd_step * DT_PS
    msd_caf_time = msd_caf_step * DT_PS

    fits = [
        linear_fit(msd_time, msd_water, 100.0, 500.0),
        linear_fit(msd_time, msd_water, 150.0, 500.0),
        linear_fit(msd_time, msd_water, 200.0, 500.0),
    ]
    primary_fit = fits[0]

    # Comparable 100 ps windows.
    direct_100 = direct["step"] <= 200000
    prod_100 = prod["step"] <= 200000
    npt_second_half = npt["step"] >= 50000

    summary = {
        "time_step_fs": 0.5,
        "short_nvt_ps": 25.0,
        "npt_ps": 50.0,
        "production_ps": 500.0,
        "full_protocol_ps": 575.0,
        "direct_100ps": {
            "temperature_K": sample_mean_sd(direct["temp"][direct_100]),
            "pressure_atm": sample_mean_sd(direct["press"][direct_100]),
            "density_g_cm3": sample_mean_sd(direct["density"][direct_100]),
            "volume_A3": sample_mean_sd(direct["vol"][direct_100]),
        },
        "post_npt_first_100ps": {
            "temperature_K": sample_mean_sd(prod["temp"][prod_100]),
            "pressure_atm": sample_mean_sd(prod["press"][prod_100]),
            "density_g_cm3": sample_mean_sd(prod["density"][prod_100]),
            "volume_A3": sample_mean_sd(prod["vol"][prod_100]),
        },
        "npt_second_half": {
            "temperature_K": sample_mean_sd(npt["temp"][npt_second_half]),
            "pressure_atm": sample_mean_sd(npt["press"][npt_second_half]),
            "density_g_cm3": sample_mean_sd(npt["density"][npt_second_half]),
            "volume_A3": sample_mean_sd(npt["vol"][npt_second_half]),
        },
        "production_500ps": {
            "temperature_K": sample_mean_sd(prod["temp"]),
            "pressure_atm": sample_mean_sd(prod["press"]),
            "density_g_cm3": sample_mean_sd(prod["density"]),
            "volume_A3": sample_mean_sd(prod["vol"]),
        },
        "diffusion_fits": fits,
        "reported_diffusion": {
            "window_ps": [100.0, 500.0],
            "D_cm2_per_s": primary_fit["D_cm2_per_s"],
            "note": "Primary 100-500 ps fit. Alternative windows are sensitivity checks and are not combined into a formal uncertainty.",
        },
        "settings": {
            "pair_modify": "mix geometric tail no",
            "minimize": "1.0e-6 1.0e-8 5000 50000",
            "shake": "1.0e-4 200 0 b 26 a 44",
            "thermostat_damping_fs": 200.0,
            "barostat_damping_fs": 1000.0,
            "barostat_coupling": "isotropic",
            "finite_size_diffusion_correction": "not applied",
        },
    }
    (ANALYSIS / "revision_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Compact data for independent regeneration outside GitHub Actions.
    compact = {"stages": []}
    for name, data, offset in [
        ("NVT curto", short, 0.0),
        ("NPT", npt, 25.0),
        ("NVT produção", prod, 75.0),
    ]:
        idx = downsample_indices(len(data["step"]))
        compact["stages"].append(
            {
                "name": name,
                "time_ps": (offset + data["step"][idx] * DT_PS).round(6).tolist(),
                "temperature_K": data["temp"][idx].round(6).tolist(),
                "pressure_atm": data["press"][idx].round(6).tolist(),
                "density_g_cm3": data["density"][idx].round(8).tolist(),
            }
        )
    (ANALYSIS / "figure6_compact.json").write_text(
        json.dumps(compact, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
    )

    # Figure 5: controlled comparison.
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.3))
    labels = ["NVT direto", "Protocolo corrigido"]
    x = np.arange(2)
    axes[0].bar(x, [0.8983, 0.9881])
    axes[0].axhline(0.997, linestyle="--", linewidth=1.2)
    axes[0].set_ylim(0.85, 1.02)
    axes[0].set_ylabel("Densidade (g cm⁻³)")
    axes[0].set_title("(a) Densidade")
    for i, value in enumerate([0.8983, 0.9881]):
        axes[0].text(i, value + 0.004, f"{value:.4f}", ha="center")

    axes[1].bar(x, [-1555.0, -77.3], yerr=[388.0, 250.8], capsize=5)
    axes[1].axhline(0.0, linewidth=0.8)
    axes[1].set_ylim(-2100, 500)
    axes[1].set_ylabel("Pressão média (atm)")
    axes[1].set_title("(b) Pressão média")
    for i, value in enumerate([-1555.0, -77.3]):
        axes[1].text(i, value - 120 if i == 0 else value + 90, f"{value:g}", ha="center")

    axes[2].bar(x, [27000.0, 24546.1])
    axes[2].set_ylim(23000, 28000)
    axes[2].set_ylabel("Volume (Å³)")
    axes[2].set_title("(c) Volume")
    for i, value in enumerate([27000.0, 24546.1]):
        axes[2].text(i, value + 140, f"{value:.0f}", ha="center")

    for axis in axes:
        axis.set_xticks(x, labels, rotation=15)
    fig.suptitle("Comparação entre os protocolos NVT direto e corrigido")
    fig.tight_layout()
    save_figure(fig, "Fig5_comparacao_protocolos_600dpi")

    # Figure 6: full 575 ps protocol.
    fig, axes = plt.subplots(3, 1, figsize=(11.2, 8.6), sharex=True)
    stage_specs = [(short, 0.0), (npt, 25.0), (prod, 75.0)]
    for data, offset in stage_specs:
        time = offset + data["step"] * DT_PS
        idx = downsample_indices(len(time), 1800)
        axes[0].plot(time[idx], data["temp"][idx], linewidth=0.6)
        axes[1].plot(time[idx], data["press"][idx], linewidth=0.6)
        axes[2].plot(time[idx], data["density"][idx], linewidth=0.8)

    for axis in axes:
        axis.axvline(25.0, linestyle="--", linewidth=1.0)
        axis.axvline(75.0, linestyle="--", linewidth=1.0)
        axis.axvspan(75.0, 175.0, alpha=0.08)
        axis.set_xlim(0.0, 575.0)
    axes[0].axhline(summary["production_500ps"]["temperature_K"][0], linestyle="--", linewidth=1.0)
    axes[1].axhline(summary["production_500ps"]["pressure_atm"][0], linestyle="--", linewidth=1.0)
    axes[2].axhline(0.997, linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("Temperatura (K)")
    axes[1].set_ylabel("Pressão (atm)")
    axes[2].set_ylabel("Densidade (g cm⁻³)")
    axes[2].set_xlabel("Tempo do protocolo (ps)")
    axes[0].set_title("(a) Temperatura")
    axes[1].set_title("(b) Pressão")
    axes[2].set_title("(c) Densidade")
    for axis in axes:
        ymax = axis.get_ylim()[1]
        axis.text(12.5, ymax, "NVT curto", ha="center", va="bottom")
        axis.text(50.0, ymax, "NPT", ha="center", va="bottom")
        axis.text(325.0, ymax, "NVT de produção", ha="center", va="bottom")
    fig.suptitle("Evolução temporal ao longo do protocolo de simulação corrigido")
    fig.tight_layout()
    save_figure(fig, "Fig6_evolucao_protocolo_600dpi")

    # Figure 7: production temperature and MSD.
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 6.6), sharex=True)
    prod_time = prod["step"] * DT_PS
    idx_prod = downsample_indices(len(prod_time), 1600)
    axes[0].plot(prod_time[idx_prod], prod["temp"][idx_prod], linewidth=0.65, label="temperatura")
    axes[0].axhline(298.0, linestyle="--", linewidth=1.0, label="alvo: 298 K")
    axes[0].axhline(summary["production_500ps"]["temperature_K"][0], linestyle=":", linewidth=1.0, label="média")
    axes[0].set_ylabel("Temperatura (K)")
    axes[0].set_title("(a) Estabilidade térmica na produção")
    axes[0].legend(ncol=3, fontsize=8)

    axes[1].plot(msd_time, msd_water, linewidth=1.1, label="água (O)")
    axes[1].plot(msd_caf_time, msd_caf, linewidth=0.9, label="cafeína")
    fit_mask = (msd_time >= 100.0) & (msd_time <= 500.0)
    fit_y = primary_fit["slope_A2_per_ps"] * msd_time[fit_mask] + primary_fit["intercept_A2"]
    axes[1].plot(msd_time[fit_mask], fit_y, linestyle="--", linewidth=1.1, label="ajuste linear da água")
    axes[1].axvline(100.0, linestyle=":", linewidth=1.0)
    axes[1].set_xlim(0.0, 500.0)
    axes[1].set_xlabel("Tempo da produção (ps)")
    axes[1].set_ylabel("MSD (Å²)")
    axes[1].set_title("(b) Mobilidade molecular e regime difusivo")
    axes[1].legend(fontsize=8)
    axes[1].text(
        0.98,
        0.08,
        f"ajuste 100–500 ps: a = {primary_fit['slope_A2_per_ps']:.3f} Å² ps⁻¹\n"
        f"D = {primary_fit['D_cm2_per_s']:.2e} cm² s⁻¹",
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )
    fig.tight_layout()
    save_figure(fig, "Fig7_msd_difusao_600dpi")


if __name__ == "__main__":
    main()
