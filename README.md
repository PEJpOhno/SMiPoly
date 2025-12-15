# SMiPoly

![license](https://anaconda.org/conda-forge/smipoly/badges/license.svg)  
|PePy.tech (pip)|conda-forge:|  
|:---:|:---:|  
|[![Downloads](https://static.pepy.tech/badge/smipoly)](https://pepy.tech/project/smipoly)|![conda](https://anaconda.org/conda-forge/smipoly/badges/downloads.svg)|  


**Documentation :** https://pejpohno.github.io/SMiPoly/  

___Demo video___  
[![Demo video](https://img.youtube.com/vi/ilzYwNWvTeQ/sddefault.jpg)](https://youtu.be/ilzYwNWvTeQ)  

## 1. About SMiPoly  
"SMiPoly (**S**mall **M**olecules **i**nto **Poly**mers)" is rule-based virtual library generator for discovery of functional polymers. It is consist of two submodules, "monc.py" and "polg.py".  
"monc.py" is a monomer classifier from a list of small molecules, and "polg.py" is a polymer repeating unit generator from the classified monomer list.  

### How To Cite (publications)   
SMiPoly: Generation of a Synthesizable Polymer Virtual Library Using Rule-Based Polymerization Reactions  
Mitsuru Ohno, Yoshihiro Hayashi, Qi Zhang, Yu Kaneko, and Ryo Yoshida  
*Journal of Chemical Information and Modeling* **2023** *63* (17), 5539&#8203;-5548  
DOI: 10.1021/acs.jcim.3c00329  
https://doi.org/10.1021/acs.jcim.3c00329  
(version 0.0.1 was used)    


## 2. Current version and requirements
current version = 1.1.2  
requirements
  - pyhon 3.10, 3.11, 3.12, 3.13, 3.14  
  - rdkit >= 2023.9.1  
  - numpy >= 1.26.0  
  - pandas >= 2.1.0  

## 3. Getting start  

### Installation  
SMiPoly can be installed with pip or conda. 
```sh
$pip install smipoly
```  
or
```sh
$conda install conda-forge::smipoly
```  

### Sample script
Download 'sample_script/sample_smip_demo4.ipynb' from [SMiPoly repository](https://github.com/PEJpOhno/SMiPoly).  
To run this demo script, Jupyter Notebook is required.

### Sample data
The sample dataset './sample_data/202207_smip_monset.csv' includes common 1,083 monomers collected from published documents such as scientific articles, catalogues and so on.

## 4. Copyright and license  
Copyright (c) 2022 Mitsuru Ohno  
Released under the BSD-3 license, license that can be found in the LICENSE file.  

## 5. Related projects  
RadonPy (Fully automated calculation for a comprehensive set of polymer properties)  
https://github.com/RadonPy/RadonPy  

## 6. THIRD-PARTY LICENSES

This software includes the following third-party libraries:

---

Library: Python
Website: https://www.python.org/
License: Python Software Foundation License Version 2
Copyright © 2001–2025 Python Software Foundation

---

Library: NumPy
Website: https://numpy.org/
License: BSD 3-Clause License
Copyright (c) 2005–2025, NumPy Developers

---

Library: pandas
Website: https://pandas.pydata.org/
License: BSD 3-Clause License
Copyright (c) 2008–2025, PyData Development Team

---

Library: RDKit
Website: https://www.rdkit.org/
License: BSD 3-Clause License
Copyright (c) 2006–2025, Rational Discovery LLC, Greg Landrum, Julie Penzotti, and contributors

