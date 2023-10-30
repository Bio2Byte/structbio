# StructBio - Structural Biology course files

<img src="https://www.vub.be/themes/custom/ocelot_vub/assets/vectors/vub-logo-icon.svg" width=200 />

Document written by [Adrián Diaz](mailto:adrian.diaz@vub.be) & [David Bickel](mailto:david.bickel@vub.be).

## Course organization
This practical course consists in 3 sessions of 2 hours each based in Jupyter Notebook files. 

- Day 1: AlphaFold predictions
  - [structbio2022_day01_part01.ipynb](https://github.com/Bio2Byte/structbio/blob/main/1_alphafold/structbio2022_day01_part01.ipynb)
  - [structbio2022_day01_part02.ipynb](https://github.com/Bio2Byte/structbio/blob/main/1_alphafold/structbio2022_day01_part02.ipynb)
- Day 2 & 3: MD Simulations
  - 1_prepare_input
    - [preparation.ipynb](https://github.com/Bio2Byte/structbio/blob/main/2_moldynamics/1_prepare_input/preparation.ipynb)
  - 2_equilibration
    - [equilibration.ipynb](https://github.com/Bio2Byte/structbio/blob/main/2_moldynamics/2_equilibration/equilibration.ipynb)
  - 3_production
  - 4_analysis
    - [analysis1.ipynb](https://github.com/Bio2Byte/structbio/blob/main/2_moldynamics/4_analysis/analysis1.ipynb)
    - [analysis2.ipynb](https://github.com/Bio2Byte/structbio/blob/main/2_moldynamics/4_analysis/analysis2.ipynb)

### Getting started

To get a copy of this repository, please follow these steps inside a terminal session:

```sh
# Change the current directory to your VSC SCRATCH directory
cd $VSC_SCRATCH

# Download course files from GitHub
git clone https://github.com/Bio2Byte/structbio.git

# Create a softlink of the directory inside your VSC HOME directory
sh structbio/hydra_setup.sh
```

### Contents

```
.
├── 1_alphafold
│   ├── examples
│   │   ├── ccdb_ccdb
│   │   │   ├── 1VUB_1_CCDB_CCDB_PAE.png
│   │   │   ├── 1VUB_1_CCDB_CCDB_plddt.png
│   │   │   └── dimer_ccdB_ccbB_relaxed_rank_1_model_5.pdb
│   │   ├── ccdb_ccdb_ccda
│   │   │   ├── ccdb_ccdb_ccda_41a4b_relaxed_rank_1_model_4.pdb
│   │   │   ├── multimer_ccdB_ccdB_ccdA_PAE.png
│   │   │   └── multimer_ccdB_ccdB_ccdA_pLDDT.png
│   │   └── ccdb_ccdb_gyr_gyr
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA.a3m
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_coverage.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA.done.txt
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_PAE.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_plddt.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_predicted_aligned_error_v1.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_1_model_3.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_2_model_1.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_3_model_2.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_4_model_5.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_5_model_4.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_template_domain_names.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_1_model_3.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_1_model_3_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_2_model_1.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_2_model_1_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_3_model_2.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_3_model_2_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_4_model_5.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_4_model_5_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_5_model_4.pdb
│   │       └── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_5_model_4_scores.json
│   ├── structbio2022_day01_part01.ipynb
│   └── structbio2022_day01_part02.ipynb
├── 2_moldynamics
│   ├── 1_prepare_input
│   │   ├── add_ions.mdp
│   │   ├── charmm36-feb2021.ff
│   │   │   ├── atomtypes.atp
│   │   │   ├── cmap.itp
│   │   │   ├── ffbonded.itp
│   │   │   ├── ffnonbonded.itp
│   │   │   ├── forcefield.doc
│   │   │   ├── forcefield.itp
│   │   │   ├── gb.itp
│   │   │   ├── ions.itp
│   │   │   ├── merged.arn
│   │   │   ├── merged.c.tdb
│   │   │   ├── merged.hdb
│   │   │   ├── merged.n.tdb
│   │   │   ├── merged.r2b
│   │   │   ├── merged.rtp
│   │   │   ├── merged.vsd
│   │   │   ├── nbfix.itp
│   │   │   ├── old_c36_cmap.itp
│   │   │   ├── spce.itp
│   │   │   ├── spc.itp
│   │   │   ├── tip3p.itp
│   │   │   ├── tip4p.itp
│   │   │   └── watermodels.dat
│   │   └── preparation.ipynb
│   ├── 2_equilibration
│   │   ├── 01_min_all.mdp
│   │   ├── 02_nvt_heat.mdp
│   │   ├── 03_npt_equi.mdp
│   │   ├── equilibration.ipynb
│   │   ├── mdrun_02_nvt_heat.slrm
│   │   └── mdrun_03_npt_equi.slrm
│   ├── 3_production
│   │   ├── mdrun.slrm
│   │   └── npt_pme.mdp
│   └── 4_analysis
│       ├── analysis1.ipynb
│       ├── analysis2.ipynb
│       └── dimer_analysis
│           ├── 1x75_superimposed.pdb
│           ├── 3hpw_superimposed.pdb
│           ├── blosum62.py
│           ├── dimer_analysis.py
│           ├── __init__.py
│           └── ndx_parser.py
├── hydra_setup.sh
└── README.md
```
