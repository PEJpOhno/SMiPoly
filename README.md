# SMiPoly

![license](https://anaconda.org/conda-forge/smipoly/badges/license.svg)  
|PePy.tech (pip)|conda-forge:|  
|:---:|:---:|  
|[![Downloads](https://static.pepy.tech/badge/smipoly)](https://pepy.tech/project/smipoly)|![conda](https://anaconda.org/conda-forge/smipoly/badges/downloads.svg)|  

## 1. What is SMiPoly?  
"SMiPoly (**S**mall **M**olecules **i**nto **Poly**mers)" is rule-based virtual library generator for discovery of functional polymers. It is consist of two submodules, "monc.py" and "polg.py".  
"monc.py" is a monomer classifier from a list of small molecules, and "polg.py" is a polymer repeating unit generator from the classified monomer list.  

## 2. Current version and requirements
current version = 0.2.1  
requirements
  - pyhon 3.7, 3.8, 3.9, 3.10, 3.11, 3.12  
  - rdkit >= 2020.09.1.0 #(2019.09.3 is unavailable)  
  - numpy >= 1.20.2  
  - pandas >= 1.2.4  

## 3. Installation and usage
### 3-1. Installatin  
SMiPoly can be installed with pip or conda. 
### 3-1-1. Install with pip  
Create new virtual environment and activate it.
To install this package, run as follows.

```sh
$pip install smipoly
```
### 3-1-2. Install with conda  

Add the channel "conda-forge" if it have not been enable.  

```sh
$conda config --add channels conda-forge
```

Create a new environment. 
```sh
$conda create -n "YOUR_NEW_ENVIRONMNT_NAME" python  
or 
$conda create -n "YOUR_NEW_ENVIRONMNT_NAME" python="required version (ex. 3.10)"
```
Then activate it. 
```sh
$conda activate "YOUR_NEW_ENVIRONMNT_NAME"
```
And install SMiPoly. 
```sh
$conda install smipoly
```

Or after create and activate a new environment, 
```sh
$conda install -c conda-forge smipoly
```

### 3-2. Quick start
Download 'sample_data/202207_smip_monset.csv' and 'sample_script/sample_smip_demo2.ipynb' from [SMiPoly repository](https://github.com/PEJpOhno/SMiPoly) to the same directry on your computer.
Then run sample_smip_demo.ipynb. To run this demo script, Jupyter Notebook is required.

## 4. Module contents  
### 4-1. monc.py  
The functions of monc.py is as follows.  
  - extract monomers from a list of small molecules.
  - classify extracted monomers into each monomer class.

The chemical structure of the small molecule compounds should be expressed in simplified molecular input line entry system (SMILES) and given as pandas DataFrame.  

**Functions**  
smip.**monc.moncls**(*df, smiColn, minFG = 2, maxFG = 4, dsp_rsl=False*)  
smip.**monc.olecls**(*df, smiColn, minFG = 1, maxFG = 4, dsp_rsl=False*) 

ARGUMENTS:  

  - df: name of the object DataFrame  
  - smicoln: The column label of the SMILES column, given as a *str*.  
  - minFG: minimum number of the polymerizable functional groups in the monomer for successive polymerization (default for moncls, 2: 2 or more; for olecls, 1: 1 or more)  
  - maxFG: maxmum nimber of the polymerizable functional groups in the monomer for successive polymerization (default 4: 4 or less)  
  - dsp_rsl: display classified result (default False)  

**Defined monomer class**  
By the function "moncls"  
  - vinylidene  
  - cyclic olefin  
  - epoxide and diepoxide  
  - lactone  
  - lactam  
  - hydroxy carboxylic acid  
  - amino acid  
  - cyclic carboxylic acid anhydride and bis(cyclic carboxylic acid anhydride)  
  - hindered phenol  
  - dicarboxylic acid and acid halide  
  - diol  
  - diamine and primary diamine  
  - diisocyanate  
  - bis(halo aryl)sulfone  
  - bis(fluoro aryl)ketone  

By the function "olecls"  
(The following class of compounds are also belong to the class "vinylidene" and / or "cyclic olefin".)  
  - acryl  
  - styryl  
  - allyl  
  - conjugated dienes  
  - vinyl ether  
  - vinyl ester  
  - maleic imide derivatives  

### 4-2. polg.py  
The library "polg.py" has two functions, "bipolym" and "ole_copolym".  

#### 4-2-1. Function "bipolym"    
The function "bipolym" gives all synthesizable polymer repeating units starting from the classified monomer list generated by "monc.moncls".  
For chain polymerization (polyolefins and some polyether), it gives homo and binary-copolymers. For successive (or step) polymerization,  it gives homopolymer only.

smip.**polg.biplym**(*df, targ = \['all'\], Pmode = 'a', dsp_rsl=False*)  

ARGUMENTS:    
  - df: name of the DataFrame of classified monomers generated by *monc.moncls*.  
  - targ: targetted polymer class. When present, it can be a list of *str*. The selectable elements are 'polyolefin', 'polyester', 'polyether', 'polyamide', 'polyimide', 'polyurethane', 'polyoxazolidone' and 'all' (default = ['all'])  
  - Pmod: generate all isomers of the polymer repeating unit ('a') or the polymer repeating unit of its representation ('r'). (default = 'a')   
  **Future warning:**  Mode 'r' will be deprecated and merged into mode 'a'. Use mode 'a' instead.  
  - dsp_rsl: display the DataFrame of the generated polymers. (default False)  

**Defined polymer class**  
  - polyolefin, polycyclic olefin and their binary copolymers  
  - polyester (from lactone, hydroxy carboxylic acid, dicarboxylic acid + diol, diol + CO and cyclic carboxylic acid anhydride + epoxide)  
  - polyether (from epoxide, hindered phenol, bis(halo aryl)sulfone + diol and bis(fluoro aryl)ketone + diol)  
  - polyamide (from lactam, amino acid and dicarboxylic acid + diamine)  
  - polyimide (bis(cyclic carboxylic acid anhydride + primary diamine)  
  - polyurethane (diisocyanate + diol)  
  - polyoxazolidone (diepoxide + diisocyanate)  

#### 4-2-2. Function "ole_copolym"  
The function "ole_copolym" gives olefinic (co)polymer repeating units starting from the classified monomer list generated by "monc.olecls". The combination of the olefinic monomer class and the number of component are given as the arguments.   

smip.**polg.ole_copolym**(df, targ = *\[ \]*, ncomp = *1*, dsp_rsl = None, drop_dupl = None)  

ARGUMENTS:  
  - df: name of the DataFrame of classified monomers generated by *monc.olecls*.  
  - targ: list of monomer class (es) to copolymerize.    
  - ncomp: Number of components given as Int. (default 1) 
  - dsp_rsl: display the DataFrame of the generated polymers. (default False)  
  - drop_dupl: drop duplicated copolymer (default False)  


### 4-3 Sample data
The sample dataset './sample_data/202207_smip_monset.csv' includes common 1,083 monomers collected from published documents such as scientific articles, catalogues and so on.

### 4-4. Utilities  
By using the files in './utilities' directory, one can modify or add the definition of monomers, the rules of polymerization reactions and polymer classes.  
To apply the new rule(s), replace the old './smipoly/rules' directory by the new one. The files must be run according to the number assigned the head of the each filename.  

  - 1_MonomerDefiner.ipynb: definitions of monomers  
  - 2_Ps_rxnL.ipynb: rules of polymerization reactions    
  - 3_Ps_GenL.ipynb: definitions of polymer classes with combinations of starting monomer(s) and polymerization reaction  

## 5. Copyright and license  
Copyright (c) 2022 Mitsuru Ohno  
Released under the BSD-3 license, license that can be found in the LICENSE file.  


## 6. Publications  
SMiPoly: Generation of a Synthesizable Polymer Virtual Library Using Rule-Based Polymerization Reactions  
Mitsuru Ohno, Yoshihiro Hayashi, Qi Zhang, Yu Kaneko, and Ryo Yoshida  
*Journal of Chemical Information and Modeling* **2023** *63* (17), 5539-5548  
DOI: 10.1021/acs.jcim.3c00329  
https://doi.org/10.1021/acs.jcim.3c00329  
(version 0.0.1 was used)    

## 7. Related projects  
RadonPy (Fully automated calculation for a comprehensive set of polymer properties)  
https://github.com/RadonPy/RadonPy  

## 8. Directry configuration  

```sh
SMiPoly
├── src
│   └── smipoly
│       ├── __init__.py
│       ├── _version.py
│       ├── smip
│       │   ├── __init__.py
│       │   ├── funclib.py
│       │   ├── monc.py
│       │   └── polg.py
│       └── rules
│           ├── excl_lst.json
│           ├── mon_dic_inv.json
│           ├── mon_dic.json
│           ├── mon_lst.json
│           ├── mon_vals.json
│           ├── ps_class.json
│           ├── ps_gen.pkl
│           └── ps.rxn.pkl
├── LICENSE
├── pyproject.toml
├── setup.py
├── setup.cfg
├── README.md
├── sample_data
│   └── 202207_smip_monset.csv
├── sample_script
│   └── sample_smip_demo2.ipynb
└── utilities
    ├── 1_MonomerDefiner.ipynb
    ├── 2_Ps_rxnL.ipynb
    ├── 3_Ps_GenL.ipynb
    └── rules/
```

## Reference  
https://future-chem.com/rdkit-chemical-rxn/  
https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html  
https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html  
