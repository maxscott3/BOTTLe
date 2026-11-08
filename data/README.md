# Data Sources

This directory contains datasets that can be used to run tau-leaping experiments, in particular 
one can recreate the results from the paper. 

Some of the CSVs are obtained from previously published research. 
The original sources are listed below.

| File | Source |
| --- | --- |
| `monosaccharide_casein_system.csv` | Brands and Boekel (2002) |
| `oscillatory_mapk_cascades.csv` | Kholodenko (2000) |
| `mapk_scaffold.csv` | Levchenko et al. (2000) |
| `caspase_activation.csv` | Legewie et al. (2006) |
| `arachidonic_acid.csv` | Yang et al. (2007) |

Additional sources can also be generated synthetically, see `synthetic_experiments.py`, or uploaded manually.


```bibtext

@article{brands2002,
	title={Kinetic modelling of reactions in heated monosaccharide−casein systems},
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
    author={Legewie, Stefan and Blüthgen, Nils and Herzel, Hanspeter},
    journal={PLoS Comp. Bio.},
    volume={2},
    number={9}
    year={2006}
}

@article{yang2007,
    title={Dynamic simulations on the arachidonic acid metabolic network},
    author={Yang, Kun and Ma, Wenzhe and Liang, Huanhuan and Ouyang, Qi and Tang, Chao and Lai, Luhua},
    journal={PLos Comp. Bio.},
    volume={3},
    number={3},
    year={2007}
}

