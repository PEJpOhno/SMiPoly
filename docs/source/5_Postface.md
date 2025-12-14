# 5. Postface  
## References  
SMARTS:  
[https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)  

SMARTS Tutorial:  
[https://www.daylight.com/dayhtml_tutorials/languages/smarts/](https://www.daylight.com/dayhtml_tutorials/languages/smarts/)  

SMARTS Examples:  
[https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html)  

SMILES:  
[https://doi.org/10.1021/ci00057a005](https://doi.org/10.1021/ci00057a005)  

At the beginning of the project, I frequently referred to the URL below.  
[https://future-chem.com/rdkit-chemical-rxn/](https://future-chem.com/rdkit-chemical-rxn/)  

## For future works  
An example of polymerization via non-classical carobocation.  
[https://doi.org/10.1295/koron.2015-0011](https://doi.org/10.1295/koron.2015-0011)  

## Directry configuration  

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
├── README.md
├── sample_data
│   |── 202207_smip_monset.csv
│   ├── SMiPoly_001_generatedPolymer.csv.zip
│   └── tutorial_aziridine.csv
├── sample_script
│   └── sample_smip_demo4.ipynb
├── utilities  
│   ├── 1_MonomerDefiner.ipynb  
│   ├── 2_Ps_rxnL.ipynb  
│   ├── 3_Ps_GenL.ipynb  
│   └── rules/  
└── docs     # source files for documentation 

```
