[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

⚠️ Note ⚠️ that this is a fork of the original repo. For the original repo please visit https://github.com/zlisysu/Weighted_cc

Please ensure that your input file contains at least three columns: The first two columns are the names of the mutation ligands, and the third colunm are their pair-wise energy. 

# Installation:
```bash
python -m pip install git+https://github.com/David-Araripe/Weighted_cc.git
```

Example:


1. Example with only three data columns: ligand1, ligand2 and pair-wise energy(Input file: Example/example_without_w).

Usage: `wccc -f Example/example_without_w -r 3A -e -8.83`

output:
```text
Node    dG_cc       path_dependent_error     path_independent_error
 3K    -9.6490             0.6800                    0.3917
 4M    -9.2534             0.6800                    0.3917
 3N    -7.2975             0.4600                    0.2542
 3A    -8.8300             0.0000                    0.2126
 3M    -9.3731             0.4600                    0.2542
 4P    -9.5708             0.9300                    0.3917
 4L    -8.2269             0.9300                    0.3917
 4N    -8.4780             0.8500                    0.2580
 4I    -7.7850             0.6800                    0.2542
 ```


2. Example with five data columns: ligand1, ligand2, pair-wise energy, bar_std, slide_std(Input file: Example/example_with_w)

Usage:
`wccc -f Example/example_with_w -r 3A -e -8.83`

output:
```text
Node    dG_cc       dG_wcc1      dG_wcc2      path_dependent_error     path_independent_error
 3K    -9.6490      -9.8057      -9.9030             0.6800                    0.3917
 4M    -9.2534      -9.1101      -9.1422             0.6800                    0.3917
 3N    -7.2975      -7.3467      -7.4097             0.4600                    0.2542
 3A    -8.8300      -8.8300      -8.8300             0.0000                    0.2126
 3M    -9.3731      -9.3372      -9.3444             0.4600                    0.2542
 4P    -9.5708      -9.6788      -9.7445             0.9300                    0.3917
 4L    -8.2269      -8.2082      -8.3992             0.9300                    0.3917
 4N    -8.4780      -8.5057      -8.5536             0.8500                    0.2580
 4I    -7.7850      -7.8019      -7.9035             0.6800                    0.2542
```

3. If you want obtain pair-wise results after calculation, add "-p yes" option

Usage: 
`wccc -f Example/example_with_w -r 3A -e -8.83 -p yes`

output:
Printing Pairwise Energies:
```text
Pair     ddG_cc      ddG_wcc1   ddG_wcc2  pair_error
3K-4M    0.3956       0.6956     0.7607    0.5600
3N-3A   -1.5325      -1.4833    -1.4203    0.4600
3M-3A    0.5431       0.5072     0.5144    0.4600
3K-3N    2.3515       2.4590     2.4932    0.5000
4P-3K   -0.0782      -0.1269    -0.1585    0.6300
4M-4P   -0.3174      -0.5688    -0.6022    0.6300
3K-4L    1.4221       1.5975     1.5038    0.6300
4M-3M   -0.1197      -0.2271    -0.2022    0.5000
4N-3K   -1.1710      -1.3001    -1.3494    0.5100
4I-3N    0.4875       0.4552     0.4938    0.5000
4I-3M   -1.5881      -1.5354    -1.4409    0.5000
4M-4L    1.0264       0.9018     0.7432    0.6300
4M-4N    0.7753       0.6043     0.5888    0.5100
****************************************************************************************************
Node    dG_cc       dG_wcc1      dG_wcc2      path_dependent_error     path_independent_error
 3K    -9.6490      -9.8057      -9.9030             0.6800                    0.3917
 4M    -9.2534      -9.1101      -9.1422             0.6800                    0.3917
 3N    -7.2975      -7.3467      -7.4097             0.4600                    0.2542
 3A    -8.8300      -8.8300      -8.8300             0.0000                    0.2126
 3M    -9.3731      -9.3372      -9.3444             0.4600                    0.2542
 4P    -9.5708      -9.6788      -9.7445             0.9300                    0.3917
 4L    -8.2269      -8.2082      -8.3992             0.9300                    0.3917
 4N    -8.4780      -8.5057      -8.5536             0.8500                    0.2580
 4I    -7.7850      -7.8019      -7.9035             0.6800                    0.2542
```
References
----------
Li, Yishui & Liu, Runduo & Liu, Jie & Luo, Haibin & Wu, Chengkun & Li, Zhe. (2022). An Open Source Graph-Based Weighted Cycle Closure Method for Relative Binding Free Energy Calculations. Journal of Chemical Information and Modeling. 
Refer to the publication above for a detailed description of the wcc method and the parameters  and please cite it to support our work if you use this software in your research.
