# 4-1. 1_MonomerDefiner.ipynb
The notebook 1_MonomerDefiner.ipynb is a utility script for defining and categorizing monomers based on their functional groups (FGs). It creates dictionaries and lists that define monomer types, their objective functional groups, and incompatible functional groups. These definitions are exported as JSON files for use in other parts of the project.

## Key Components
**1. Metadata**  
The notebook includes metadata about its purpose and references:

- Purpose: Define monomers and their functional groups.
- References: Links to SMARTS theory and examples for understanding chemical patterns.  

**2. Functional Group (FG) Definitions**  
The notebook defines objective and incompatible functional groups for each monomer type. These definitions are stored in two dictionaries:

- monL: Contains SMARTS patterns for objective functional groups.
- exclL: Contains SMARTS patterns for incompatible functional groups.

**3. Monomer Dictionaries**  
The notebook defines dictionaries to map monomer names to numerical identifiers:

- mon_dic: 
Python dictionary of the monomer classes. Categorizes the values of "mon_dic.json" into self-polymerizable monomer systems, binary-monomer systems, post-polymerization reactions, and olefinic monomer systems.
    - key : monomer class  
    - value : Allocated int.  
        - 1 - 50 : for self-polymerizable system; *eg.* addition, ring-opening polymerization and self-condensation.
        - 51- 100 : for the polymerization system that require bi-monomer classes; *eg.* poly condensation, polyadditin and addition-condensation.  
        - 1000 - : detailed classification for olefinic monomers. 
- mon_dic_inv: Inverse mapping of mon_dic (from integers to names).
- mon_vals: A tuple categorized mon_dic values according to the corresponding polymerization reactions (see Table 5-1, 2, 3). 

**4. Monomer Definitions**  
Each monomer type is defined with:  
- mon_lst.json (monL)  
List of the definition of the polymerization site for each monomer class. Contains SMARTS patterns for objective functional groups.  
    - 1 - 50 : for self-polymerizable system; *eg.* addition, ring-opening polymerization and self-condensation.
    - 51- 100 : for the polymerization system that require bi-monomer classes; *eg.* poly condensation, polyadditin and addition-condensation.  
    - 200 -  post-polymerization reactions of residual polymerizable functional grous(s). 
    - 1000 - for olefinic monomer systems. 

- excl_lst.json (exclL)   
List of the definition of functional groups that should not coexist in a monomer molecule for each monomer class. Contains SMARTS patterns for incompatible functional groups.  

**5. Exporting Definitions**  
The monomer definitions are exported as JSON files for use in other scripts. 

## Correspondence Table of Numbers and Content
In following tables, 
- Allocated numbers are corresponding to 
    - values of mon_dic except 200-999
    - key of mon_dic_inv except 200-999
    - key of mon_lst and excl_lst
 
- Monomer classes are corresponding to corresponding to key of mon_dic and values of mon_dic_inv. 

**Table 4-1. Allocated numbers and defined monomer classes for self-polymerizable system** 
|No.    |monomer class| compounds                                       |
|-------|-------------|----------------------------------------------------|
| 1     | vinyl       | vinylidene                                         |
| 2     | epo         | epoxide                                            |
| 3     | cOle        | cyclic olefin                                      |
| 4     | lactone     | lactone                                            |
| 5     | lactam      | lactam                                             |
| 6     | hydCOOH     | hydroxy carboxylic acid                            |
| 7     | aminCOOH    | amino acid                                         |
| 8     | hindPhenol  | hindered phenol                                    |
| 9     | cAnhyd      | cyclic carboxylic acid anhydride                   |
| 10    | CO          | carbon monoxide. Forced addition to synthesize carbonates|
| 11    | HCHO        | form aldehyde. Forced addition for addition-condensation (future works) |
| 12    | sfonediX    | bis(p-halo aryl)sulfone                              |
| 13    | BzodiF      | bis(p-fluoro aryl)ketone                             | 


**Table 4-2. Allocated numbers and defined monomer classes for the polymerization system that require bi-monomer classes**  
|No.    |monomer class| compounds                                                    |
|-------|-------------|--------------------------------------------------------------|
| 51    | diepo       | di/polyepoxide                                               |
| 52    | diCOOH      | di/polycarboxylic acid and acid halide                       |
| 53    | diol        | primary and secodary di/polyol, di/polythiol and /polyphenol |
| 54    | diamin      | di/polyamine                                                 |
| 55    | diNCO       | di/polyisocyanate                                            |
| 56    | dicAnhyd    | bis/poly(cyclic carboxylic acid anhydride)                   |
| 57    | pridiamin   | di/polyprimary diamine                                       |
| 58    | diol_b      | primary and secodary di/polyol, di/polyphenol                |


**Table 4-3. Allocated numbers and defined monomer classes for olefinic monomers** 
|No.    |monomer class| compounds                                          |
|-------|-------------|----------------------------------------------------|
| 1001  | acryl       | acrylate                                           |
| 1002  | bEWole      | beta-ectron withdrawing group substituted olefin   |
| 1003  | styryl      | styryl                                             |
| 1004  | allyl       | allyl                                              |
| 1005  | haloCH      | halogenated olefin                                 |
| 1006  | vinylester  | vinyl ester                                        |
| 1007  | malei       | maleic imide derivatives                           |
| 1020  | conjdiene   | conjugated dienes                                  |
| 1030  | vinylether  | vinyl ether                                        |
| 1031  | tertcatCH   | beta-disubstituded aliphatic olefin                |
| 1050  | cycCH       | alycyclic olefin                                   |
| 1052  | aliphCH     | aliphatic olefin                                   |


**Table 4-4. post-polymerization reactions of residual polymerizable functional grous(s)**  
|No.  |SMARTS       | Residual functional group                          |
|-----|-------------|----------------------------------------------------|
| 200 | \[CX3\]=\[CX3\]                | olefin                          |
| 201 | \[CX4;R\]1\[OX2;R\]\[CX4;R\]1  | epoxide                         | 
| 202 | \[CX3\](=\[O\])\[OX2H1,F,Cl,Br,I\] | carboxylic acid and acyl halide |
| 203 | \[C,c\]\[OX2,SX2;H1;!$([O,S]C=*)\] | hydroxyl in alcohol and phenol  |
| 204 | \[C,c\]\[NX3;H2;!$(N\[C,S\]=*)\] | amine                         | 
| 205 | \[NX2\]=\[CX2\]=\[OX1,SX1\]    | isocyanate                      |
| 206 | \[C,c\]\[CX3,c;R\](=\[OX1\])\[OX2,o;R\]\[CX3,c;R\](=\[OX1\])\[C,c\]| carboxylic acid anhydride |


