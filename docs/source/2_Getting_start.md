# 2. Getting Start

This guide will walk you through the steps to prepare input files, classify monomers, and generate polymers using the SMiPoly framework.  
When using Jupyter Notebook, also refer to "sample_script/sample_smip_demo4.ipynb".  

---

## 1. Quick start  
For quick start, the following sampl script and dataset are available.  

**Sample script**  
Download ['./sample_script/sample_smip_demo4.ipynb'](https://github.com/PEJpOhno/SMiPoly/blob/main/sample_script/sample_smip_demo4.ipynb) from [SMiPoly repository](https://github.com/PEJpOhno/SMiPoly.git).  
To run this demo script, Jupyter Notebook is required.

**Sample data**  
The sample dataset './sample_data/202207_smip_monset.csv' includes common 1,083 monomers collected from published documents such as scientific articles, catalogues and so on.


## 2. Prepare and load Input File

To begin, you need to prepare an input file containing the chemical data for your compounds. The input file should be in a tabular format (e.g., CSV) with at least one column containing SMILES strings. 

### Steps:
1. Ensure the file is saved in a format readable by pandas (e.g., `.csv` or `.xlsx`).
2. The column containing SMILES strings should be clearly labeled (e.g., `SMILES`).  

### Example Input File
| Compound_ID | SMILES                    |
|-------------|---------------------------|
| CID174      | OCCO                      |
| CID7489     | C1=CC(=CC=C1C(=O)O)C(=O)O |
| CID6658     | O=C(OC)C(C)=C             |
| CID7501     | C=CC1=CC=CC=C1            |
| CID1140     | CC1=CC=CC=C1              |
| CID702      | CCO                       |
| CID7896     | CC(CCO)O                  |
| CID7837     | CC(=C)C(=O)OCC1CO1        |
| CID10352    |  C1CC2CC1C=C2             |


```python
import pandas as pd

df = pd.read_csv("input_file.csv")
```


### Using Sample dataset  
Also  'sample_data/202207_smip_monset.csv' is avalable under the Internet-connected environment.  
```python
import pandas as pd

df = pd.read_csv("https://raw.githubusercontent.com/PEJpOhno/SMiPoly/main/sample_data/202207_smip_monset.csv")  
```

---

## 3. Monomer Classification

The monomer classification step identifies and categorizes monomers based on their functional groups or olefinic properties.

### Monomer Extraction and Classification
Use the `moncls` function to extract and classify monomers from small molecule compounds.

#### Example:
```python
from smipoly.smip import monc

# Classify monomers
classified_monomer_df = monc.moncls(df, smiColn="SMILES", minFG=2, maxFG=4, dsp_rsl=True)

# Save results
classified_monomer_df.to_csv("classified_monomers.csv", index=False)
```

### Olefinic Monomer Classification
To classify olefinic monomers in more detail, use the `olecls` function.

#### Example:
```python
from smipoly.smip import monc

# Classify olefinic monomers
classified_olefins_df = monc.olecls(df, smiColn="SMILES", minFG=1, maxFG=4, dsp_rsl=True)

# Save results
classified_olefins_df.to_csv("classified_olefins.csv", index=False)
```

---

## 4. Generate Polymers

The polymer generation step involves creating polymers from the classified monomers.

### Polymer Generation
Use the `biplym` function to generate polymers.

#### Example:
```python
from smipoly.smip import polg

# Generate polymer CRU
generated_polymers_df = polg.biplym(classified_monomer_df, targ=["polyester"], dsp_rsl=True)

# Save results
generated_polymers_df.to_csv("generated_polymers_df.csv", index=False)
```

### Olefinic (Co)polymer Generation
Once you have applied the `olecls` function to classify olefin monomers, use the `ole_copolym` function to generate olefin (co)polymers. 

#### Example:
```python
from smipoly.smip import polg

# Generate olefinic (co)polymers
olefinic_polymers_df = polg.ole_copolym(classified_olefins_df, targ=["ROMP"], ncomp=1, dsp_rsl=True)

# Save results
olefinic_polymers_df.to_csv("olefinic_polymers_df.csv", index=False)
```

**Notes**  
When installed with pip or conda, the following will be executed automatically.  
- Ensure all required dependencies (e.g., RDKit, pandas) are installed in your environment.  
- The rules directory must contain the necessary configuration files (e.g., mon_vals.json, ps_rxn.pkl).