# Data

This directory contains the reaction systems used to run tau-leaping experiments, including 
those used to produce the results in the paper.


## File format

Each reaction system with $N$ species and $M$ reactions is stored as a headerless CSV file with 
$M + 1$ rows and $N + 1$ columns:

- **First row:** the first entry is unused (set to 0), followed by the initial counts of the $N$ 
  species.
- **Remaining $M$ rows:** one per reaction. The first entry is the reaction rate constant, 
  followed by the net stoichiometry vector, i.e. the change in each species count when the 
  reaction fires.

For example, the row `0.0001,-1,-1,1,0` describes a reaction with rate constant $10^{-4}$ that 
consumes one molecule each of species 1 and 2 and produces one molecule of species 3. Files can 
be read and written using `read_csv` and `write_to_csv` in `utils.py`.


## Sources

Some of the CSVs are obtained from previously published research, with the original sources 
listed below.

| File | Source |
| --- | --- |
| `monosaccharide_casein_system.csv` | Brands and van Boekel (2002) |
| `oscillatory_mapk_cascades.csv` | Kholodenko (2000) |
| `mapk_scaffold.csv` | Levchenko et al. (2000) |
| `caspase_activation.csv` | Legewie et al. (2006) |
| `arachidonic_acid.csv` | Yang et al. (2007) |
| `michaelis_menten.csv` | Standard Michaelis-Menten enzyme kinetics |

The remaining datasets were generated synthetically using the functions in 
`synthetic_experiments.py` with the following arguments:

| File | Call |
| --- | --- |
| `protein_aggregation.csv` | `protein_aggregation(X_0=5000, max_length=150, params=np.array([1e-5, 0.05, 0.01]), ...)` |
| `morphogen_gradient.csv` | `morphogen_gradient(X_0=0, num_bins=200, params=np.array([200, 5, 0.01]), ...)` |

where `params` gives the nucleation, elongation and depolymerisation rates, and the production, 
diffusion and degradation rates respectively.


```bibtex
@article{brands2002,
    title={Kinetic modelling of reactions in heated monosaccharide-casein systems},
    author={Brands, Carline M. and van Boekel, Martinus A.},
    journal={J. Agric. Food Chem.},
    volume={50},
    number={23},
    pages={6725--6739},
    year={2002}
}

@article{kholodenko2000,
    title={Negative feedback and ultrasensitivity can bring about oscillations in the mitogen-activated protein kinase cascades},
    author={Kholodenko, Boris N.},
    journal={Euro. J. Biochem.},
    volume={267},
    number={6},
    pages={1583--1588},
    year={2000}
}

@article{levchenko2000,
    title={Scaffold proteins may biphasically affect the levels of mitogen-activated protein kinase signaling and reduce its threshold properties},
    author={Levchenko, Andre and Bruck, Jehoshua and Sternberg, Paul W.},
    journal={Proc. Natl. Acad. Sci.},
    volume={97},
    number={11},
    pages={5818--5823},
    year={2000}
}

@article{legewie2006,
    title={Mathematical modeling identifies inhibitors of apoptosis as mediators of positive feedback and bistability},
    author={Legewie, Stefan and Bl{\"u}thgen, Nils and Herzel, Hanspeter},
    journal={PLoS Comp. Bio.},
    volume={2},
    number={9},
    year={2006}
}

@article{yang2007,
    title={Dynamic simulations on the arachidonic acid metabolic network},
    author={Yang, Kun and Ma, Wenzhe and Liang, Huanhuan and Ouyang, Qi and Tang, Chao and Lai, Luhua},
    journal={PLoS Comp. Bio.},
    volume={3},
    number={3},
    year={2007}
}
```
