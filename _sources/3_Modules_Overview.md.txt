# 3. Modules Overview

## 1. monc.py
The monc.py module is part of the SMiPoly project and is responsible for categorizing monomers based on their functional groups (FGs) and other chemical properties. It processes a dataset of chemical compounds (in SMILES format) and classifies them into predefined monomer classes. The module also supports olefinic monomer classification and provides utilities for handling functional group patterns.  

### Loading Rules and Data
The module loads several pre-defined rules and data files from the rules directory:

- mon_vals.json:  
Categorizes the values of "mon_dic.json" into self-polymerizable monomer systems, binary-monomer systems, post-polymerization reactions, and olefinic monomer systems.
- mon_dic.json: Maps monomer class to numerical identifiers.
- mon_dic_inv.json: Inverse mapping of mon_dic.
- mon_lst.json: Contains SMARTS patterns for objective functional groups.
- excl_lst.json: Contains SMARTS patterns for incompatible functional groups.
- ps_rxn.pkl:Contains SMARTS patterns for incompatible functional groups. ContainsPython dictionary with mon_dic's key and value swapped polymerization reaction definitions.  

These files are loaded into dictionaries (mon_vals, mon_dic, mon_dic_inv, monL, exclL, Ps_rxnL) and converted to integer-keyed dictionaries (monLg, exclLg).

### Functions
- **moncls**  
Purpose: Classifies monomers into predefined classes based on functional groups.

- **olecls**  
Purpose: Classifies olefinic monomers into predefined classes.

### Key Features
- Monomer Classification:
    - Uses SMARTS patterns to identify functional groups in molecules.
    - Supports both mono-functionalized and poly-functionalized monomers.

- Olefin Classification:
    - Handles olefinic monomers, including conjugated dienes.
    - Refines classification results using polymerization reactions.

- Integration with Rules:
    - Relies on external JSON and pickle files for functional group definitions and reaction rules.
    - Modular design allows easy updates to rules and patterns.
- DataFrame Processing:
    - Processes chemical data stored in pandas DataFrames.
    - Appends classification results as new columns.  


## 2. polg.py
The module polg.py is a Python script designed for generating polymers from classified monomers. It is part of a larger project and relies on external data files and helper functions from other modules, such as funclib.py.  
The primary purpose of polg.py is to generate polymers based on input monomer data and predefined polymerization rules.  Below is a deContains SMARTS patterns for incompatible functional groups.tailed explanation of its components and functionality:

### Key Components
#### Data Loading
The module loads several JSON and pickle files from a rules directory. These files contain essential data for polymer generation:

- mon_vals.json: Contains categorized monomer values.
- mon_dic.json and mon_dic_inv.json: Dictionaries mapping monomer classes and their inverse mappings.
- mon_lst.json and excl_lst.json: Lists of monomers and exclusion rules.
- ps_rxn.pkl: Contains reaction rules for polymerization.
- ps_class.json: Defines polymer classes.
- ps_gen.pkl: Contains polymer generation rules.

The data is loaded into variables like mon_vals, mon_dic, monL, exclL, Ps_rxnL, Ps_classL, and Ps_GenL.

#### Helper Variables
The module processes the loaded data into dictionaries like monLg and exclLg, which map integer keys to monomer and exclusion lists.

### Functions
**biplym(df, targ=None, dsp_rsl=None)**  
This function generates polymer CRU based on an input DataFrame and specified target polymer classes.

- Key Steps:
    1. Filters and processes the input DataFrame to identify valid monomers.
    2. Generates polymers based on predefined polymerization rules and target classes.
    3. Removes duplicate polymerization reactions and adjusts the resulting DataFrame.
    4. Optionally displays the number of polymerization reactions and generated polymers.

- Notes:
    - Uses helper functions like bipolymA and homopolymA from funclib.py for polymer generation.
    - Handles both binary polymerization and homopolymerization.

**ole_copolym(df, targ=None, ncomp=None, dsp_rsl=None, drop_dupl=None)**  
This function generates olefinic (co)polymers based on olefin classes and parameters.

- Key Steps:
    1. Validates the input olefin classes and parameters.
    2. Reconstructs the list of classified olefin monomers for polymerization.
    3. Generates (co)polymers by combining monomers based on the specified number of components (ncomp).
    4. Optionally drops duplicate copolymers and displays the results summary.

- Special Handling:  
    - Handles specific olefin classes like ROMP, ROMPH, and COC with unique logic.
    - Uses helper functions like coord_polym from funclib.py.

### Key Features
1. Rule-based Polymer Generation:
    - The module relies heavily on external data files (JSON, pickle) to define monomer structure and polymerization rules.

2. Integration with RDKit:
    - RDKit is used for chemical structure handling, including molecule generation, substructure matching, and reaction processing.

3. Customizability:
    - Users can specify target polymer classes, olefin classes, and other parameters to customize the polymer generation process.  
    - User can custermize the rules in advanced usage. 
