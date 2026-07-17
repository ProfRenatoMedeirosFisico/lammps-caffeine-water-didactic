# Protocolo didático de dinâmica molecular: sistema cafeína–água com LAMMPS

Arquivos de entrada, dados e logs associados ao artigo:

> **Uso do LAMMPS no ensino de dinâmica molecular: protocolo didático para o sistema cafeína–água**
> Submetido à *Revista Brasileira de Ensino de Física* (RBEF), 2026.

## Descrição

Repositório com o material completo para reproduzir o protocolo de simulação descrito no artigo. O sistema consiste em uma molécula de cafeína solvatada por 800 moléculas de água em uma caixa cúbica de ~30 × 30 × 30 Å³, campo de força OPLS-AA (cafeína) + SPC/E rígido via SHAKE (água).

O protocolo segue cinco etapas com funções físicas distintas:

| Etapa | Ensemble | Duração | Finalidade |
|-------|----------|---------|-----------|
| 1. Minimização | — | convergência | Remover contatos desfavoráveis |
| 2. NVT curto | NVT 150→298 K | 50 ps | Acomodação térmica com rampa |
| 3. NPT | NPT 298 K, 1 atm | 500 ps | Relaxação volumétrica e ajuste de densidade |
| 4. NVT produção | NVT 298 K | 500 ps | Trajetória principal; cálculo de MSD e D |
| 5. NVT direto (controlado) | NVT 298 K | 100 ps | Comparação isolada sem NPT (mesma semente, mesmo Tdamp) |

## Requisitos

- **LAMMPS** (testado com `20230802.3`): https://www.lammps.org
- **Packmol**: https://m3g.ims.unicamp.br/packmol
- **Python 3** com `numpy` e `matplotlib` (análise pós-simulação)

## Campo de força

| Componente | Modelo | Fonte |
|---|---|---|
| Cafeína | OPLS-AA (24 tipos atômicos) | LigParGen — Dodda et al., *JPCB* **121**, 3864 (2017) |
| Água | SPC/E rígida (SHAKE) | Berendsen et al., *J. Phys. Chem.* **91**, 6269 (1987) |

Parâmetros globais (em `in.common`):

- Combinação: geométrica (`pair_modify mix geometric`)
- Escala 1-4: `special_bonds lj/coul 0.0 0.0 0.5`
- Cutoff: 10 Å; eletrostática de longo alcance: PPPM (10⁻⁴)
- Passo de integração: 0,5 fs; SHAKE: tol=10⁻⁴, 200 iter, ligação O-H e ângulo H-O-H

## Arquivos do repositório

### Entradas e topologias

| Arquivo | Descrição |
|---------|-----------|
| `in.common` | Parâmetros globais incluídos por todos os inputs |
| `packmol_caffeine_water.inp` | Input do Packmol para montagem da caixa |
| `caffeine_clean.pdb` | Estrutura PDB da cafeína |
| `water.pdb` | Estrutura PDB da água SPC/E |
| `system.data` | Topologia inicial gerada pelo Packmol + fftool/LigParGen |
| `system_final_opls_spce.data` | Topologia corrigida (tipos atômicos OPLS-AA finais) |
| `in.nvt_prod_long_v2` | Input da produção NVT longa (500 ps + MSD) — **novo** |
| `in.nvt_direct_controlled` | Input do NVT direto controlado (comparação sem NPT) — **novo** |

### Dados de reinício (checkpoints)

| Arquivo | Etapa de origem |
|---------|----------------|
| `minimized_opls_spce.data` | Pós-minimização |
| `nvt_short_final.data` | Pós-NVT curto |
| `npt_eq_final.data` | Pós-NPT |
| `nvt_prod_final.data` | Pós-produção original (100 ps) |
| `nvt_prod_long_final.data` | Pós-produção longa (500 ps) — **novo** |
| `system_initial.pdb` | Configuração inicial (pós-Packmol) |

### Logs termodinâmicos

| Arquivo | Etapa |
|---------|-------|
| `log.min_serial.lammps` | Minimização |
| `log.minimization.lammps` | Minimização (versão alternativa) |
| `log.nvt_short.lammps` | NVT curto |
| `log.nvt_serial.lammps` | NVT curto (execução serial de referência) |
| `log.npt_eq.lammps` | NPT equilíbrio |
| `log.nvt_prod.lammps` | NVT produção original (100 ps) |
| `log.nvt_prod_long.lammps` | NVT produção longa (500 ps) — **novo** |
| `log.nvt_direct_controlled.lammps` | NVT direto controlado (100 ps) — **novo** |

## Como reproduzir

### 1. Montagem da caixa (opcional — `.data` já incluído)

```bash
packmol < packmol_caffeine_water.inp
```

### 2. Minimização

```bash
lmp -in in.minimization -log log.minimization.lammps
```

### 3. NVT curto (aquecimento 150 → 298 K)

```bash
lmp -in in.nvt_short -log log.nvt_short.lammps
```

### 4. Equilíbrio NPT (298 K, 1 atm)

```bash
lmp -in in.npt_eq -log log.npt_eq.lammps
```

### 5. Produção NVT longa (500 ps + MSD)

```bash
lmp -in in.nvt_prod_long_v2 -log log.nvt_prod_long.lammps
```

Gera: `msd_waterO.dat`, `msd_caffeine.dat`, `nvt_prod_long_final.data`

### 6. NVT direto controlado (comparação sem NPT)

Parte de `nvt_short_final.data`, mesmo Tdamp=200 fs, mesma temperatura alvo.

```bash
lmp -in in.nvt_direct_controlled -log log.nvt_direct_controlled.lammps
```

### 7. Réplica independente

Repetir as etapas 3–5 alterando a semente em `in.nvt_short` (linha `velocity all create 150.0 SEMENTE`). A réplica do artigo usou semente `8452197`.

## Resultados obtidos

- **Coeficiente de difusão da água (SPC/E):** D = 2,65 × 10⁻⁵ cm² s⁻¹ (janela 100–500 ps; referência SPC/E: ~2,5 × 10⁻⁵ cm² s⁻¹)
- **Densidade pós-NPT:** ρ = 0,997 ± 0,003 g cm⁻³ (referência SPC/E: 0,997 g cm⁻³)
- **Temperatura de produção:** T = 298,2 ± 0,5 K
- **Pressão NVT direto (sem NPT):** P = −1555 ± 388 atm (evidência da necessidade do NPT)

## Referência

Se utilizar este material, por favor cite:

> R. A. Medeiros et al., Uso do LAMMPS no ensino de dinâmica molecular: protocolo didático para o sistema cafeína–água, *Revista Brasileira de Ensino de Física* (2026). [DOI a ser atribuído após publicação]

## Licença

[MIT](LICENSE)
