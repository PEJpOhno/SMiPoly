**WARNING:** Please note that this document is located in the develop branch and is still in progress 

# 1. About SMiPoly  
"SMiPoly (**S**mall **M**olecules **i**nto **Poly**mers)" is rule-based virtual library generator for discovery of functional polymers. It is consist of two submodules, "monc.py" and "polg.py".  
"monc.py" is a monomer classifier from a list of small molecules, and "polg.py" is a polymer repeating unit generator from the classified monomer list.  

**How To Cite (publications)**    
SMiPoly: Generation of a Synthesizable Polymer Virtual Library Using Rule-Based Polymerization Reactions  
Mitsuru Ohno, Yoshihiro Hayashi, Qi Zhang, Yu Kaneko, and Ryo Yoshida  
*Journal of Chemical Information and Modeling* **2023** *63* (17), 5539-5548  
DOI: 10.1021/acs.jcim.3c00329  
<a href="https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00329">https://doi.org/10.1021/acs.jcim.3c00329</a>  
(version 0.0.1 was used)    

**Project Link**  
url = <a href="https://github.com/PEJpOhno/SMiPoly">https://github.com/PEJpOhno/SMiPoly</a>    

## Directry configuration of SMiPoly  

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
│   └── sample_smip_demo3.ipynb
└── utilities  
    ├── 1_MonomerDefiner.ipynb  
    ├── 2_Ps_rxnL.ipynb  
    ├── 3_Ps_GenL.ipynb  
    └── rules/  

```
