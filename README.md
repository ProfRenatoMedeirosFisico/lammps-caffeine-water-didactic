# Protocolo didático de dinâmica molecular: sistema cafeína–água com LAMMPS

Arquivos de entrada, dados e scripts de análise associados ao artigo:

> **Uso do LAMMPS no ensino de dinâmica molecular: protocolo didático para o sistema cafeína–água**
>
> Submetido à *Revista Brasileira de Ensino de Física* (RBEF).

## Descrição

Este repositório contém o material necessário para reproduzir integralmente o protocolo de simulação de dinâmica molecular descrito no artigo. O sistema consiste em uma molécula de cafeína solvatada por 800 moléculas de água em uma caixa cúbica de aproximadamente 30 × 30 × 30 Å³.

O protocolo segue cinco etapas com funções físicas distintas:

1. **Minimização** — remoção de contatos desfavoráveis na configuração inicial.
2. **NVT curto** — acomodação térmica com aquecimento gradual de 150 a 298 K.
3. **NPT de equilíbrio** — relaxação volumétrica e ajuste de densidade a 298 K e 1 atm.
4. **NVT de produção** — trajetória principal para análise de observáveis.
5. **Réplica independente** — verificação de reprodutibilidade com semente diferente.


## Requisitos

- **LAMMPS** (versão 22 Jul 2025 ou posterior): https://www.lammps.org
- **Packmol**: https://m3g.ims.unicamp.br/packmol
- **Python 3** com `numpy` e `matplotlib` (para análise)

## Campo de força

| Componente | Modelo | Fonte |
|---|---|---|
| Cafeína | OPLS-AA (24 tipos atômicos) | LigParGen (Dodda et al., *Nucleic Acids Res.* **45**, W331, 2017) |
| Água | SPC/E (rígida, SHAKE) | Berendsen et al., *J. Phys. Chem.* **91**, 6269, 1987 |

Parâmetros adicionais:

- Regras de combinação: geométrica (`pair_modify mix geometric`)
- Escala 1-4: 0,5 para LJ e Coulomb (`special_bonds lj/coul 0.0 0.0 0.5`)
- Cutoff: 10 Å (LJ e real-space Coulomb)
- Eletrostática de longo alcance: PPPM (precisão 10⁻⁴)
- Correção de cauda LJ: desativada (`tail no`)
- Passo de integração: 0,5 fs

## Como reproduzir

### 1. Montagem da caixa (opcional — o arquivo `.data` já está incluído)

```bash
cd packmol/
packmol < pack_system.inp
```

### 2. Minimização

```bash
cd lammps/
lammps -in in.minimization
```

Gera: `minimized_opls_spce.data`

### 3. NVT curto

```bash
lammps -in in.nvt_short
```

Gera: `nvt_short_final.data`

### 4. Equilíbrio NPT

```bash
lammps -in in.npt_eq
```

Gera: `npt_eq_final.data`

### 5. Produção NVT

```bash
lammps -in in.nvt_prod
```

Gera: `nvt_prod_final.data`, trajetória e logs termodinâmicos.

### 6. Réplica independente

Repetir as etapas 3–5 alterando a semente de velocidades no arquivo `in.nvt_short` (linha `velocity all create 150.0 SEMENTE`). No artigo, a réplica usou semente `8452197`.

### 7. Análise

```bash
cd analysis/
python plot_observables.py
```

## Referência

Se utilizar este material, por favor cite:

> [Autor(es)], Uso do LAMMPS no ensino de dinâmica molecular: protocolo didático para o sistema cafeína–água, *Revista Brasileira de Ensino de Física* (2026). [DOI a ser atribuído]

## Licença

Este material é distribuído sob a licença [MIT](LICENSE).
