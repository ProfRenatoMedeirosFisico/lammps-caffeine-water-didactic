# Protocolo didÃ¡tico de dinÃ¢mica molecular: sistema cafeÃ­naâ€“Ã¡gua com LAMMPS

Arquivos de entrada, dados e logs associados ao artigo:

> **Uso do LAMMPS no ensino de dinÃ¢mica molecular: protocolo didÃ¡tico para o sistema cafeÃ­naâ€“Ã¡gua**
> Submetido Ã  *Revista Brasileira de Ensino de FÃ­sica* (RBEF), 2026.

## DescriÃ§Ã£o

RepositÃ³rio com o material completo para reproduzir o protocolo de simulaÃ§Ã£o descrito no artigo. O sistema consiste em uma molÃ©cula de cafeÃ­na solvatada por 800 molÃ©culas de Ã¡gua em uma caixa cÃºbica de ~30 Ã— 30 Ã— 30 Ã…Â³, campo de forÃ§a OPLS-AA (cafeÃ­na) + SPC/E rÃ­gido via SHAKE (Ã¡gua).

O protocolo segue cinco etapas com funÃ§Ãµes fÃ­sicas distintas:

| Etapa | Ensemble | DuraÃ§Ã£o | Finalidade |
|-------|----------|---------|-----------|
| 1. MinimizaÃ§Ã£o | â€” | convergÃªncia | Remover contatos desfavorÃ¡veis |
| 2. NVT curto | NVT 150â†’298 K | 50 ps | AcomodaÃ§Ã£o tÃ©rmica com rampa |
| 3. NPT | NPT 298 K, 1 atm | 500 ps | RelaxaÃ§Ã£o volumÃ©trica e ajuste de densidade |
| 4. NVT produÃ§Ã£o | NVT 298 K | 500 ps | TrajetÃ³ria principal; cÃ¡lculo de MSD e D |
| 5. NVT direto (controlado) | NVT 298 K | 100 ps | ComparaÃ§Ã£o isolada sem NPT (mesma semente, mesmo Tdamp) |

## Requisitos

- **LAMMPS** (testado com `20230802.3`): https://www.lammps.org
- **Packmol**: https://m3g.ims.unicamp.br/packmol
- **Python 3** com `numpy` e `matplotlib` (anÃ¡lise pÃ³s-simulaÃ§Ã£o)

## Campo de forÃ§a

| Componente | Modelo | Fonte |
|---|---|---|
| CafeÃ­na | OPLS-AA (24 tipos atÃ´micos) | LigParGen â€” Dodda et al., *JPCB* **121**, 3864 (2017) |
| Ãgua | SPC/E rÃ­gida (SHAKE) | Berendsen et al., *J. Phys. Chem.* **91**, 6269 (1987) |

ParÃ¢metros globais (em `in.common`):

- CombinaÃ§Ã£o: geomÃ©trica (`pair_modify mix geometric`)
- Escala 1-4: `special_bonds lj/coul 0.0 0.0 0.5`
- Cutoff: 10 Ã…; eletrostÃ¡tica de longo alcance: PPPM (10â»â´)
- Passo de integraÃ§Ã£o: 0,5 fs; SHAKE: tol=10â»â´, 200 iter, ligaÃ§Ã£o O-H e Ã¢ngulo H-O-H

## Arquivos do repositÃ³rio

### Entradas e topologias

| Arquivo | DescriÃ§Ã£o |
|---------|-----------|
| `in.common` | ParÃ¢metros globais incluÃ­dos por todos os inputs |
| `packmol_caffeine_water.inp` | Input do Packmol para montagem da caixa |
| `caffeine_clean.pdb` | Estrutura PDB da cafeÃ­na |
| `water.pdb` | Estrutura PDB da Ã¡gua SPC/E |
| `system.data` | Topologia inicial gerada pelo Packmol + fftool/LigParGen |
| `system_final_opls_spce.data` | Topologia corrigida (tipos atÃ´micos OPLS-AA finais) |
| `in.nvt_prod_long_v2` | Input da produÃ§Ã£o NVT longa (500 ps + MSD) â€” **novo** |
| `in.nvt_direct_controlled` | Input do NVT direto controlado (comparaÃ§Ã£o sem NPT) â€” **novo** |

### Dados de reinÃ­cio (checkpoints)

| Arquivo | Etapa de origem |
|---------|----------------|
| `minimized_opls_spce.data` | PÃ³s-minimizaÃ§Ã£o |
| `nvt_short_final.data` | PÃ³s-NVT curto |
| `npt_eq_final.data` | PÃ³s-NPT |
| `nvt_prod_final.data` | PÃ³s-produÃ§Ã£o original (100 ps) |
| `nvt_prod_long_final.data` | PÃ³s-produÃ§Ã£o longa (500 ps) â€” **novo** |
| `system_initial.pdb` | ConfiguraÃ§Ã£o inicial (pÃ³s-Packmol) |

### Logs termodinÃ¢micos

| Arquivo | Etapa |
|---------|-------|
| `log.min_serial.lammps` | MinimizaÃ§Ã£o |
| `log.minimization.lammps` | MinimizaÃ§Ã£o (versÃ£o alternativa) |
| `log.nvt_short.lammps` | NVT curto |
| `log.nvt_serial.lammps` | NVT curto (execuÃ§Ã£o serial de referÃªncia) |
| `log.npt_eq.lammps` | NPT equilÃ­brio |
| `log.nvt_prod.lammps` | NVT produÃ§Ã£o original (100 ps) |
| `log.nvt_prod_long.lammps` | NVT produÃ§Ã£o longa (500 ps) â€” **novo** |
| `log.nvt_direct_controlled.lammps` | NVT direto controlado (100 ps) â€” **novo** |

## Como reproduzir

### 1. Montagem da caixa (opcional â€” `.data` jÃ¡ incluÃ­do)

```bash
packmol < packmol_caffeine_water.inp
```

### 2. MinimizaÃ§Ã£o

```bash
lmp -in in.minimization -log log.minimization.lammps
```

### 3. NVT curto (aquecimento 150 â†’ 298 K)

```bash
lmp -in in.nvt_short -log log.nvt_short.lammps
```

### 4. EquilÃ­brio NPT (298 K, 1 atm)

```bash
lmp -in in.npt_eq -log log.npt_eq.lammps
```

### 5. ProduÃ§Ã£o NVT longa (500 ps + MSD)

```bash
lmp -in in.nvt_prod_long_v2 -log log.nvt_prod_long.lammps
```

Gera: `msd_waterO.dat`, `msd_caffeine.dat`, `nvt_prod_long_final.data`

### 6. NVT direto controlado (comparaÃ§Ã£o sem NPT)

Parte de `nvt_short_final.data`, mesmo Tdamp=200 fs, mesma temperatura alvo.

```bash
lmp -in in.nvt_direct_controlled -log log.nvt_direct_controlled.lammps
```

### 7. RÃ©plica independente

Repetir as etapas 3â€“5 alterando a semente em `in.nvt_short` (linha `velocity all create 150.0 SEMENTE`). A rÃ©plica do artigo usou semente `8452197`.

## Resultados obtidos

- **Coeficiente de difusÃ£o da Ã¡gua (SPC/E):** D = 2,65 Ã— 10â»âµ cmÂ² sâ»Â¹ (janela 100â€“500 ps; referÃªncia SPC/E: ~2,5 Ã— 10â»âµ cmÂ² sâ»Â¹)
- **Densidade pÃ³s-NPT:** Ï = 0,997 Â± 0,003 g cmâ»Â³ (referÃªncia SPC/E: 0,997 g cmâ»Â³)
- **Temperatura de produÃ§Ã£o:** T = 298,2 Â± 0,5 K
- **PressÃ£o NVT direto (sem NPT):** P = âˆ’1555 Â± 388 atm (evidÃªncia da necessidade do NPT)

## ReferÃªncia

Se utilizar este material, por favor cite:

> R. A. Medeiros et al., Uso do LAMMPS no ensino de dinÃ¢mica molecular: protocolo didÃ¡tico para o sistema cafeÃ­naâ€“Ã¡gua, *Revista Brasileira de Ensino de FÃ­sica* (2026). [DOI a ser atribuÃ­do apÃ³s publicaÃ§Ã£o]

## LicenÃ§a

[MIT](LICENSE)

<!-- FULL_500PS_DATA_START -->
### Arquivos completos da produção NVT de 500 ps

A trajetória e os dados abaixo correspondem à execução usada nas análises de MSD, difusão e nas Figuras 5 e 6. A produção empregou `1 000 000` de passos com `timestep = 0,5 fs`, totalizando `500 ps`.

| Arquivo | Conteúdo |
|---|---|
| `log.nvt_prod_long.lammps` | Log completo da produção, encerrado normalmente no passo 1 000 000 |
| `thermo_avg_prod_long.dat` | Médias termodinâmicas a cada 1 000 passos |
| `msd_waterO.dat` | MSD dos oxigênios da água |
| `msd_caffeine.dat` | MSD da cafeína |
| `nvt_prod_long.zip` | Trajetória exata `nvt_prod_long.lammpstrj`, compactada |

SHA-256:

```text
log.nvt_prod_long.lammps
9bce39671764b8e9650e8215cc67c1505ec41bbac0b6b604fa2c6524ed61b035

thermo_avg_prod_long.dat
7fea3bbe178a6fe3f771cd4bc9f99799afb8516c127d35d957861b1c75d90c42

nvt_prod_long.zip
932066e7963a3bea0729bc2a106f68911cb3c1cd2f8e40e2a3ab74e72b5fa53a
```
<!-- FULL_500PS_DATA_END -->

