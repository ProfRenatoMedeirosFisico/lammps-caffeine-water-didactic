# Protocolo didático de dinâmica molecular: sistema cafeína–água com LAMMPS

Arquivos de entrada, estruturas, dados processados, trajetória e scripts associados ao manuscrito RBEF-2026-0227.

## Sistema e objetivo

O sistema contém uma molécula de cafeína e 800 moléculas de água em uma caixa cúbica inicial de aproximadamente 30 × 30 × 30 Å³. A cafeína usa parâmetros OPLS-AA obtidos com LigParGen e a água usa o modelo SPC/E rígido por SHAKE. O material foi organizado para ensinar a diferença entre estabilidade numérica e consistência física, com ênfase no papel dos ensembles NVT e NPT.

## Protocolo efetivamente utilizado

| Etapa | Ensemble | Passos | Duração | Finalidade |
|---|---|---:|---:|---|
| Minimização | — | até 5000 iterações | — | Remover contatos desfavoráveis |
| NVT curto | NVT, 150→298 K e manutenção | 50000 | 25 ps | Acomodação térmica |
| NPT | NPT isotrópico, 298 K e 1 atm | 100000 | 50 ps | Relaxação do volume e da densidade |
| NVT produção | NVT, 298 K | 1000000 | 500 ps | Produção, MSD e estimativa de D |
| NVT direto controlado | NVT, 298 K | 200000 | 100 ps | Controle sem NPT, iniciado na mesma condição pós-NVT curto |

O protocolo corrigido possui duração acumulada de **575 ps**: 25 ps de NVT curto, 50 ps de NPT e 500 ps de produção.

## Configurações relevantes

- LAMMPS: `units real`, `atom_style full`, condições periódicas.
- Interações não ligadas: `lj/cut/coul/long 10.0 10.0` e PPPM com precisão `1.0e-4`.
- Mistura: `pair_modify mix geometric tail no`.
- Não foi aplicada correção analítica de cauda de Lennard-Jones.
- Escala 1–4: `special_bonds lj/coul 0.0 0.0 0.5`.
- Minimização: `minimize 1.0e-6 1.0e-8 5000 50000` com gradiente conjugado.
- SHAKE da água: `fix ... shake 1.0e-4 200 0 b 26 a 44`.
- Passo de integração: 0,5 fs. Somente a água é rígida; as ligações C–H da cafeína permanecem flexíveis.
- Termostato de Nosé–Hoover: `Tdamp = 200 fs`.
- Barostato isotrópico: `iso 1.0 1.0 1000.0`, com `Pdamp = 1000 fs`.
- Inicialização do protocolo principal: 150 K, distribuição gaussiana e `mom yes` para remover o momento linear líquido.

A combinação OPLS-AA/SPC/E por regra geométrica de mistura é uma aproximação prática; os parâmetros cruzados cafeína–água não foram ajustados conjuntamente.

## Arquivos de entrada

| Arquivo | Função |
|---|---|
| `in.common` | Parâmetros globais |
| `in.minimization` | Minimização por gradiente conjugado |
| `in.nvt_short` | Rampa 150→298 K e manutenção em NVT |
| `in.npt_eq` | Equilibração NPT isotrópica de 50 ps |
| `in.nvt_prod_long_v2` | Produção NVT de 500 ps, MSD e dados termodinâmicos |
| `in.nvt_direct_controlled` | Comparação controlada de 100 ps sem NPT |

## Dados e estruturas

- `system_final_opls_spce.data`: topologia inicial parametrizada.
- `minimized_opls_spce.data`: estrutura pós-minimização.
- `nvt_short_final.data`: estrutura e velocidades após NVT curto.
- `npt_eq_final.data`: último frame da etapa NPT usado na produção.
- `nvt_prod_long_final.data`: estrutura após a produção de 500 ps.
- `nvt_prod_long.zip`: trajetória exata da produção longa, armazenada por Git LFS.
- `thermo_avg_prod_long.dat`: série termodinâmica processada da produção.
- `msd_waterO.dat` e `msd_caffeine.dat`: deslocamentos quadráticos médios.
- `log.minimization.lammps`: log completo da minimização usada para preparar `minimized_opls_spce.data`.

## Resultados usados no manuscrito revisado

### Comparação em janelas equivalentes de 100 ps

- NVT direto: ρ = 0,8983 g cm⁻³; V = 27000 Å³; P = −1555 ± 388 atm.
- NVT após NPT: ρ = 0,9881 g cm⁻³; V = 24546,1 Å³; P = −77,3 ± 250,8 atm.

### Segunda metade do NPT

- T = 298,5 ± 6,3 K.
- P = 2,2 ± 492,8 atm.
- ρ = 0,9932 ± 0,0089 g cm⁻³.
- V = 24423 ± 218 Å³.

### Produção de 500 ps

O manuscrito usa um único conjunto consistente de estatísticas para a produção completa:

- T = 297,97 ± 5,91 K.
- P = −90,19 ± 439,15 atm.
- ρ = 0,9881 g cm⁻³.
- V = 24546,1 Å³.

### Difusão

O valor final é obtido pela janela principal de 100–500 ps:

- inclinação: a = 1,590 Å² ps⁻¹;
- D = 2,65 × 10⁻⁵ cm² s⁻¹.

As janelas 150–500 e 200–500 ps são usadas apenas como testes de sensibilidade e não são combinadas como uma incerteza estatística, pois compartilham a maior parte dos mesmos dados. O valor de D é **não corrigido para tamanho finito**; não foi aplicada a correção de Yeh–Hummer.

## Figuras e análise reprodutível

O script `scripts/generate_revision_assets.py` recalcula as estatísticas da produção, ajusta o MSD e gera as Figuras 5–7 em três formatos:

- PNG a 600 dpi;
- PDF;
- SVG.

A Figura 6 foi corrigida para **não reutilizar** a pressão do ramo NVT direto (−1555 ± 388 atm) como se fosse a média do NVT curto. O painel de pressão mostra somente os resumos estacionários reportados para a segunda metade do NPT e para a produção de 500 ps.

Execução:

```bash
python -m pip install numpy matplotlib
python scripts/generate_revision_assets.py
```

Os resultados numéricos são gravados em `analysis/revision_summary.json`.

## Reprodução das etapas

```bash
lmp -in in.minimization -log log.minimization.lammps
lmp -in in.nvt_short -log log.nvt_short.lammps
lmp -in in.npt_eq -log log.npt_eq.lammps
lmp -in in.nvt_prod_long_v2 -log log.nvt_prod_long.lammps
lmp -in in.nvt_direct_controlled -log log.nvt_direct_controlled.lammps
```

## Arquivos históricos que não definem a comparação revisada

Alguns arquivos antigos permanecem no histórico do repositório porque foram gerados durante a construção do protocolo. Eles **não devem ser usados para reconstruir a comparação final do manuscrito**:

- `log.nvt_serial.lammps`: execução diagnóstica anterior, com inicialização direta a 298 K e `Tdamp = 100 fs`;
- `log.nvt_prod.lammps`: produção anterior de 100 ps;
- `log.min_serial.lammps`: cópia histórica da minimização; o log canônico correspondente é agora `log.minimization.lammps`.

A comparação final deve usar `in.nvt_direct_controlled`/`log.nvt_direct_controlled.lammps` e os primeiros 100 ps da produção de 500 ps.

## Réplica independente

O manuscrito também discute uma segunda réplica estatística iniciada com semente distinta. Os valores reportados são mantidos como verificação complementar. Os arquivos brutos originais dessa réplica não foram recuperados neste snapshot do repositório; por isso, os arquivos históricos da raiz **não devem ser apresentados como substitutos da Réplica 2**. Quando os arquivos originais forem recuperados, devem ser depositados em uma pasta própria (`replica_2/`) para preservar a proveniência dos números reportados.

## Licença

MIT — consulte `LICENSE`.
