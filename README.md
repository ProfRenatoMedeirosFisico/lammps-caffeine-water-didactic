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
- Inicialização: 150 K, distribuição gaussiana e `mom yes` para remover o momento linear líquido.

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

Cada réplica foi iniciada a partir do último frame de sua própria etapa NPT. Como o volume instantâneo flutua, os volumes finais das réplicas diferem ligeiramente. O uso do volume médio equilibrado e posterior reescalonamento da célula é uma alternativa mais controlada para futuras aplicações.

## Resultados usados no manuscrito revisado

### Comparação em janelas equivalentes de 100 ps

- NVT direto: ρ = 0,8983 g cm⁻³; V = 27000 Å³; P = −1555 ± 388 atm.
- NVT após NPT: ρ = 0,9881 g cm⁻³; V = 24546,1 Å³; P = −77,3 ± 250,8 atm.

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

- `figures/Fig5_comparacao_protocolos.svg`: comparação corrigida com −1555 e −77,3 atm.
- `figures/Fig6_protocolo_575ps.svg`: eixo acumulado correto até 575 ps.
- `scripts/generate_revision_assets.py`: recalcula estatísticas, ajusta o MSD e exporta as Figuras 5–7 em PNG de 600 dpi e PDF vetorial.

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

## Licença

MIT — consulte `LICENSE`.
