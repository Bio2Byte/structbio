# StructBio
<img src="https://www.vub.be/themes/custom/ocelot_vub/assets/vectors/vub-logo-icon.svg" width=200 />
Structural Biology course files

Document written by Adrián Diaz ([@agdiaz](https://github.com/agdiaz/)) and David Bickel ([@bickeld](https://github.com/bickeld/))

## Course organization
This practical course consists in 3 sessions of 2 hours each based in Jupyter Notebook files. 

- Day 1: AlphaFold predictions
  - [structbio2022_day01_part01.ipynb](https://github.com/Bio2Byte/structbio/blob/main/day_01/structbio2022_day01_part01.ipynb)
  - [structbio2022_day01_part02.ipynb](https://github.com/Bio2Byte/structbio/blob/main/day_01/structbio2022_day01_part02.ipynb)
- Day 2: MD Simulations
  - 1_prepare_input
    - [preparation.ipynb](https://github.com/Bio2Byte/structbio/blob/main/day_02/1_prepare_input/preparation.ipynb)
  - 2_equilibration
    - [equilibration.ipynb](https://github.com/Bio2Byte/structbio/blob/main/day_02/2_equilibration/equilibration.ipynb)
- Day 3: MD Simulations
  - TBD

### Quick start
To get a copy of this repository, please follow these steps inside a terminal session:

```console
# Change the current directory to your VSC user's SCRATCH directory
cd $VSC_SCRATCH

git clone https://github.com/Bio2Byte/structbio.git

# Create a shortcut of that copy inside your VSC user's HOME directory
ln -s $VSC_SCRATCH/structbio $VSC_HOME/structbio
```

### Structure

```bash
.
├── README.md
├── day_01
│   ├── examples
│   │   ├── ccdb_ccdb
│   │   │   ├── 1VUB_1_CCDB_CCDB_PAE.png
│   │   │   ├── 1VUB_1_CCDB_CCDB_plddt.png
│   │   │   └── dimer_ccdB_ccbB_relaxed_rank_1_model_5.pdb
│   │   ├── ccdb_ccdb_ccda
│   │   │   ├── ccdb_ccdb_ccda_41a4b_relaxed_rank_1_model_4.pdb
│   │   │   ├── multimer_ccdB_ccdB_ccdA_PAE.png
│   │   │   └── multimer_ccdB_ccdB_ccdA_pLDDT.png
│   │   └── ccdb_ccdb_gyr_gyr
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA.a3m
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA.done.txt
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_PAE.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_coverage.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_plddt.png
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_predicted_aligned_error_v1.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_1_model_3.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_2_model_1.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_3_model_2.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_4_model_5.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_relaxed_rank_5_model_4.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_template_domain_names.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_1_model_3.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_1_model_3_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_2_model_1.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_2_model_1_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_3_model_2.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_3_model_2_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_4_model_5.pdb
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_4_model_5_scores.json
│   │       ├── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_5_model_4.pdb
│   │       └── multimer_ccdB_ccdB_gyrA_gyrA_unrelaxed_rank_5_model_4_scores.json
│   ├── structbio2022_day01_part01.ipynb
│   └── structbio2022_day01_part02.ipynb
└── day_02
    ├── 1_prepare_input
    │   ├── add_ions.mdp
    │   ├── charmm36-feb2021.ff
    │   │   ├── atomtypes.atp
    │   │   ├── cmap.itp
    │   │   ├── ffbonded.itp
    │   │   ├── ffnonbonded.itp
    │   │   ├── forcefield.doc
    │   │   ├── forcefield.itp
    │   │   ├── gb.itp
    │   │   ├── ions.itp
    │   │   ├── merged.arn
    │   │   ├── merged.c.tdb
    │   │   ├── merged.hdb
    │   │   ├── merged.n.tdb
    │   │   ├── merged.r2b
    │   │   ├── merged.rtp
    │   │   ├── merged.vsd
    │   │   ├── nbfix.itp
    │   │   ├── old_c36_cmap.itp
    │   │   ├── spc.itp
    │   │   ├── spce.itp
    │   │   ├── tip3p.itp
    │   │   ├── tip4p.itp
    │   │   └── watermodels.dat
    │   └── preparation.ipynb
    ├── 2_equilibration
    │   ├── 01_min_all.mdp
    │   ├── 02_nvt_heat.mdp
    │   ├── 03_npt_equi.mdp
    │   ├── equilibration.ipynb
    │   ├── mdrun_02_nvt_heat.slrm
    │   └── mdrun_03_npt_equi.slrm
    └── 3_production
        ├── mdrun.slrm
        └── npt_pme.mdp
```
