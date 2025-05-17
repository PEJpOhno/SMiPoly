# SMiPoly

![license](https://anaconda.org/conda-forge/smipoly/badges/license.svg)  
|PePy.tech (pip)|conda-forge:|  
|:---:|:---:|  
|[![Downloads](https://static.pepy.tech/badge/smipoly)](https://pepy.tech/project/smipoly)|![conda](https://anaconda.org/conda-forge/smipoly/badges/downloads.svg)|  

## 1. About SMiPoly  
"SMiPoly (**S**mall **M**olecules **i**nto **Poly**mers)" is rule-based virtual library generator for discovery of functional polymers. It is consist of two submodules, "monc.py" and "polg.py".  
"monc.py" is a monomer classifier from a list of small molecules, and "polg.py" is a polymer repeating unit generator from the classified monomer list.  

### How To Cite (publications)   
SMiPoly: Generation of a Synthesizable Polymer Virtual Library Using Rule-Based Polymerization Reactions  
Mitsuru Ohno, Yoshihiro Hayashi, Qi Zhang, Yu Kaneko, and Ryo Yoshida  
*Journal of Chemical Information and Modeling* **2023** *63* (17), 5539-5548  
DOI: 10.1021/acs.jcim.3c00329  
https://doi.org/10.1021/acs.jcim.3c00329  
(version 0.0.1 was used)    

## 2. Current version and requirements
current version = 1.1.0  
requirements
  - pyhon 3.9, 3.10, 3.11, 3.12  
  - rdkit >= 2023.9.1  
  - numpy >= 1.26.0  
  - pandas >= 2.1.0  

## 3. Installation  
SMiPoly can be installed with pip or conda. 

## 4. Getting start  
## SMiPoly documentation  
https://pejpohno.github.io/SMiPoly/  

## Sample script
Download 'sample_data/202207_smip_monset.csv' and 'sample_script/sample_smip_demo2.ipynb' from [SMiPoly repository](https://github.com/PEJpOhno/SMiPoly) to the same directry on your computer.
Then run sample_smip_demo.ipynb. To run this demo script, Jupyter Notebook is required.

## Sample data
The sample dataset './sample_data/202207_smip_monset.csv' includes common 1,083 monomers collected from published documents such as scientific articles, catalogues and so on.

## Utilities  
By using the files in './utilities' directory, one can modify or add the definition of monomers, the rules of polymerization reactions and polymer classes.  
To apply the new rule(s), replace the old './smipoly/rules' directory by the new one. The files must be run according to the number assigned the head of the each filename.  

  - 1_MonomerDefiner.ipynb: definitions of monomers  
  - 2_Ps_rxnL.ipynb: rules of polymerization reactions    
  - 3_Ps_GenL.ipynb: definitions of polymer classes with combinations of starting monomer(s) and polymerization reaction  

## 5. Copyright and license  
Copyright (c) 2022 Mitsuru Ohno  
Released under the BSD-3 license, license that can be found in the LICENSE file.  

## 6. Related projects  
RadonPy (Fully automated calculation for a comprehensive set of polymer properties)  
https://github.com/RadonPy/RadonPy  
