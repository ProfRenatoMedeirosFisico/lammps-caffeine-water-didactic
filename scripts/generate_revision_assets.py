#!/usr/bin/env python3
"""Regenerate the statistics and Figures 5-7 used in the revised manuscript.

The script uses the deposited thermodynamic table and MSD files. It therefore does
not depend on ``log.nvt_prod_long.lammps`` and remains reproducible even when only
the processed production table and exact trajectory are available.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
ANALYSIS = ROOT / "analysis"
FIGURES.mkdir(exist_ok=True)
ANALYSIS.mkdir(exist_ok=True)
DT_PS = 0.0005  # 0.5 fs


def load_thermo(path: Path) -> dict[str, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    return {
        "step": data[:, 0],
        "temp": data[:, 1],
        "press": data[:, 2],
        "density": data[:, 3],
        "pe": data[:, 4],
        "ke": data[:, 5],
        "etotal": data[:, 6],
        "vol": data[:, 7],
    }


def load_msd(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    return data[:, 0] * DT_PS, data[:, 1]


def mean_sd(values: np.ndarray) -> tuple[float, float]:
    return float(np.mean(values)), float(np.std(values, ddof=1))


def fit_diffusion(time_ps: np.ndarray, msd: np.ndarray, start: float, end: float) -> dict[str, float]:
    mask = (time_ps >= start) & (time_ps <= end)
    slope, intercept = np.polyfit(time_ps[mask], msd[mask], 1)
    fitted = slope * time_ps[mask] + intercept
    ss_res = float(np.sum((msd[mask] - fitted) ** 2))
    ss_tot = float(np.sum((msd[mask] - np.mean(msd[mask])) ** 2))
    return {
        "window_ps": [start, end],
        "slope_A2_ps": float(slope),
        "D_cm2_s": float(slope / 6.0 * 1.0e-4),
        "R2": 1.0 - ss_res / ss_tot,
    }


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def figure5() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.35))
    x = np.arange(2)
    labels = ["NVT direto", "Protocolo corrigido"]

    axes[0].bar(x, [0.8983, 0.9881])
    axes[0].axhline(0.997, linestyle="--", linewidth=1.2, label="SPC/E: 0,997")
    axes[0].set_ylim(0.85, 1.02)
    axes[0].set_ylabel("Densidade (g cm⁻³)")
    axes[0].set_title("(a) Densidade")
    axes[0].legend(fontsize=8)
    for index, value in enumerate([0.8983, 0.9881]):
        axes[0].text(index, value + 0.004, f"{value:.4f}".replace(".", ","), ha="center")

    axes[1].bar(x, [-1555.0, -77.3], yerr=[388.0, 250.8], capsize=5)
    axes[1].axhline(0.0, linewidth=0.8)
    axes[1].set_ylim(-2100, 500)
    axes[1].set_ylabel("Pressão média (atm)")
    axes[1].set_title("(b) Pressão média")
    axes[1].text(0, -1730, "−1555", ha="center")
    axes[1].text(1, 50, "−77,3", ha="center")

    axes[2].bar(x, [27000.0, 24546.1])
    axes[2].set_ylim(23000, 28000)
    axes[2].set_ylabel("Volume (Å³)")
    axes[2].set_title("(c) Volume")
    for index, value in enumerate([27000.0, 24546.1]):
        axes[2].text(index, value + 150, f"{value:.0f}", ha="center")

    for axis in axes:
        axis.set_xticks(x, labels, rotation=15)
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("Comparação entre os protocolos NVT direto e corrigido")
    fig.tight_layout()
    save(fig, "Fig5_comparacao_protocolos_600dpi")


def figure6(prod_stats: dict[str, tuple[float, float]]) -> None:
    """Synthesis of the 575 ps corrected protocol without reusing direct-NVT statistics.

    The previous schematic incorrectly used the controlled direct-NVT pressure
    (-1555 ± 388 atm) as if it characterized the short-NVT stage. The revised
    figure reports pressure only for phases for which the manuscript gives a
    consistent stationary summary: the second half of NPT and the 500 ps
    production. The NPT density change is shown as a schematic transition.
    """
    fig, axes = plt.subplots(3, 1, figsize=(11.3, 8.8), sharex=True)

    # (a) Protocol temperature: 5 ps ramp, 20 ps hold, then reported NPT and production summaries.
    axes[0].plot([0, 5], [150, 298], linewidth=1.5, label="rampa 150→298 K")
    axes[0].plot([5, 25], [298, 298], linewidth=1.5, label="manutenção NVT: 298 K")
    axes[0].plot([25, 75], [298.5, 298.5], linewidth=1.5)
    axes[0].fill_between([50, 75], 298.5 - 6.3, 298.5 + 6.3, alpha=0.16,
                         label="NPT, 2ª metade: 298,5 ± 6,3 K")
    axes[0].plot([75, 575], [prod_stats["temperature_K"][0]] * 2, linewidth=1.5)
    axes[0].fill_between(
        [75, 575],
        prod_stats["temperature_K"][0] - prod_stats["temperature_K"][1],
        prod_stats["temperature_K"][0] + prod_stats["temperature_K"][1],
        alpha=0.16,
        label=f"produção: {prod_stats['temperature_K'][0]:.2f} ± {prod_stats['temperature_K'][1]:.2f} K".replace(".", ","),
    )
    axes[0].axhline(298.0, linestyle="--", linewidth=1.0)
    axes[0].set_ylim(145, 312)
    axes[0].set_ylabel("Temperatura (K)")
    axes[0].set_title("(a) Protocolo térmico e faixas de flutuação")
    axes[0].legend(fontsize=7, ncol=2, loc="lower right")

    # (b) Pressure: only stationary summaries actually reported in the manuscript.
    axes[1].fill_between([50, 75], 2.2 - 492.8, 2.2 + 492.8, alpha=0.18)
    axes[1].plot([50, 75], [2.2, 2.2], linewidth=1.5, label="NPT, 2ª metade: 2,2 ± 492,8 atm")
    axes[1].fill_between(
        [75, 575],
        prod_stats["pressure_atm"][0] - prod_stats["pressure_atm"][1],
        prod_stats["pressure_atm"][0] + prod_stats["pressure_atm"][1],
        alpha=0.18,
    )
    axes[1].plot(
        [75, 575], [prod_stats["pressure_atm"][0]] * 2, linewidth=1.5,
        label=f"produção: {prod_stats['pressure_atm'][0]:.2f} ± {prod_stats['pressure_atm'][1]:.2f} atm".replace(".", ","),
    )
    axes[1].axhline(0.0, linestyle="--", linewidth=1.0)
    axes[1].text(12.5, -560, "NVT curto:\ncaixa subdensa\nρ = 0,8983 g cm⁻³", ha="center", va="center", fontsize=7)
    axes[1].set_ylim(-700, 650)
    axes[1].set_ylabel("Pressão (atm)")
    axes[1].set_title("(b) Pressão nas fases estacionárias")
    axes[1].legend(fontsize=8, loc="lower right")

    # (c) Density: initial fixed box, schematic NPT relaxation, final fixed-volume production.
    axes[2].plot([0, 25], [0.8983, 0.8983], linewidth=1.5, label="NVT curto: 0,8983 g cm⁻³")
    axes[2].plot([25, 75], [0.8983, 0.9881], linewidth=1.5, label="NPT: relaxação do volume")
    axes[2].plot([75, 575], [prod_stats["density_g_cm3"][0]] * 2, linewidth=1.5,
                 label="produção: 0,9881 g cm⁻³")
    axes[2].axhline(0.997, linestyle="--", linewidth=1.0, label="referência SPC/E ≈ 0,997")
    axes[2].set_ylim(0.88, 1.015)
    axes[2].set_ylabel("Densidade (g cm⁻³)")
    axes[2].set_xlabel("Tempo acumulado do protocolo (ps)")
    axes[2].set_title("(c) Ajuste de densidade e volume fixo na produção")
    axes[2].legend(fontsize=8, ncol=2, loc="lower right")

    for axis in axes:
        axis.axvline(25, linestyle=":", linewidth=1.0)
        axis.axvline(75, linestyle=":", linewidth=1.0)
        axis.axvspan(75, 175, alpha=0.05)
        axis.set_xlim(0, 575)
        axis.grid(alpha=0.18)
        ytop = axis.get_ylim()[1]
        axis.text(12.5, ytop, "NVT curto\n25 ps", ha="center", va="top", fontsize=7)
        axis.text(50, ytop, "NPT\n50 ps", ha="center", va="top", fontsize=7)
        axis.text(325, ytop, "NVT de produção\n500 ps", ha="center", va="top", fontsize=7)

    fig.suptitle("Síntese temporal do protocolo corrigido (duração total: 575 ps)")
    fig.tight_layout()
    save(fig, "Fig6_protocolo_575ps_600dpi")


def figure7(thermo: dict[str, np.ndarray], water_time: np.ndarray, water_msd: np.ndarray, caf_time: np.ndarray, caf_msd: np.ndarray, primary: dict[str, float]) -> None:
    prod_time = thermo["step"] * DT_PS
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 6.6), sharex=True)
    axes[0].plot(prod_time, thermo["temp"], linewidth=0.65, label="temperatura")
    axes[0].axhline(298.0, linestyle="--", linewidth=1.0, label="alvo: 298 K")
    axes[0].axhline(np.mean(thermo["temp"]), linestyle=":", linewidth=1.0, label="média")
    axes[0].set_ylabel("Temperatura (K)")
    axes[0].set_title("(a) Estabilidade térmica na produção")
    axes[0].legend(ncol=3, fontsize=8)

    axes[1].plot(water_time, water_msd, linewidth=1.1, label="água (O)")
    axes[1].plot(caf_time, caf_msd, linewidth=0.9, label="cafeína")
    mask = (water_time >= 100.0) & (water_time <= 500.0)
    intercept = np.mean(water_msd[mask] - primary["slope_A2_ps"] * water_time[mask])
    axes[1].plot(water_time[mask], primary["slope_A2_ps"] * water_time[mask] + intercept, linestyle="--", linewidth=1.1, label="ajuste linear da água")
    axes[1].axvline(100.0, linestyle=":", linewidth=1.0)
    axes[1].set_xlim(0.0, 500.0)
    axes[1].set_xlabel("Tempo da produção (ps)")
    axes[1].set_ylabel("MSD (Å²)")
    axes[1].set_title("(b) Mobilidade molecular e regime difusivo")
    axes[1].legend(fontsize=8)
    axes[1].text(0.98, 0.08, f"ajuste 100–500 ps: a = {primary['slope_A2_ps']:.3f} Å² ps⁻¹\nD = {primary['D_cm2_s']:.2e} cm² s⁻¹", transform=axes[1].transAxes, ha="right", va="bottom")
    fig.tight_layout()
    save(fig, "Fig7_msd_difusao_600dpi")


def main() -> None:
    thermo = load_thermo(ROOT / "thermo_avg_prod_long.dat")
    water_time, water_msd = load_msd(ROOT / "msd_waterO.dat")
    caf_time, caf_msd = load_msd(ROOT / "msd_caffeine.dat")

    fits = [fit_diffusion(water_time, water_msd, start, 500.0) for start in (100.0, 150.0, 200.0)]
    production = {
        "temperature_K": mean_sd(thermo["temp"]),
        "pressure_atm": mean_sd(thermo["press"]),
        "density_g_cm3": mean_sd(thermo["density"]),
        "volume_A3": mean_sd(thermo["vol"]),
    }
    summary = {
        "protocol_duration_ps": {"short_nvt": 25.0, "npt": 50.0, "production": 500.0, "total": 575.0},
        "controlled_100ps": {"direct_pressure_atm": [-1555.0, 388.0], "post_npt_pressure_atm": [-77.3, 250.8]},
        "npt_second_half": {"temperature_K": [298.5, 6.3], "pressure_atm": [2.2, 492.8], "density_g_cm3": [0.9932, 0.0089], "volume_A3": [24423.0, 218.0]},
        "production_500ps": production,
        "diffusion_primary": fits[0],
        "diffusion_sensitivity": fits[1:],
        "method_settings": {
            "pair_modify": "mix geometric tail no",
            "minimize": "1.0e-6 1.0e-8 5000 50000",
            "shake": "1.0e-4 200 0 b 26 a 44",
            "thermostat_damping_fs": 200.0,
            "barostat": "iso 1.0 1.0 1000.0",
            "finite_size_correction": "not applied",
        },
    }
    (ANALYSIS / "revision_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    figure5()
    figure6(production)
    figure7(thermo, water_time, water_msd, caf_time, caf_msd, fits[0])


if __name__ == "__main__":
    main()
