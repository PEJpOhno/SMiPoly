# 4-2. 2_Ps_rxnL.ipynb  

The notebook 2_Ps_rxnL.ipynb is a utility script for defining a list of polymerization reactions using RDKit's SMARTS-based reaction definitions. It creates a dictionary Ps_rxnL that maps monomer types or reaction types to their corresponding chemical reactions. These reactions are later serialized into a pickle file (ps_rxn.pkl) for use in other parts of the project.

## Key Components
**Metadata**  
The notebook includes metadata about its purpose and references:

- Purpose: Define a list of polymerization reactions.
- References: Links to RDKit documentation and SMARTS theory for understanding reaction definitions.

**Imports**  
The following libraries are imported:

- json: For loading monomer dictionary (mon_dic.json).
- pickle: For saving the reaction dictionary (Ps_rxnL) to a file.
- rdkit.Chem.AllChem: For creating chemical reactions from SMARTS strings.

**Loading Monomer Dictionary**  
The monomer dictionary (mon_dic.json) is loaded to map monomer types (e.g., vinyl, epo) to numerical identifiers:

**Reaction Definitions**  
The notebook defines a series of polymerization reactions using RDKit's AllChem.ReactionFromSmarts function. These reactions are stored in the Ps_rxnL dictionary, where the keys are either numerical identifiers or descriptive names, and the values are RDKit reaction objects.　　
Note that polymerization reactions SMARS using detailed olefinic monomer classification except ROMP(H) and COC were automatically generated from corrsponding monomer by the function "smipoly.smip.funclib.ole_rxnsmarts_gen".  

**Reaction Categories**  
1. Self-polymerizable Monomer Systems:
- Self-polymerization reactions to form homopolymers for monomers like vinyl, epoxide, cyclic olefins, lactones, lactams, etc.

2. Binary Monomer Systems:

- Olefin copolymerization reactions to form altanating copolymers by combinations of monomers (e.g., vinyl-vinyl, vinyl-cyclic olefin).  
- Polymerization reaction which require essentialy bi-monomer classes (e.g. dicarboxylic acid + diol)

3. Sequential Reactions:

- Post-polymerization reactions for residual specific functional groups (e.g., epoxide, hydroxyl, carboxyl, amine).

4. Olefinic Monomers:

- Polymerization reactions using detailed olefinic monomer classification  
- Reactions for cyclic olefins (e.g., ROMP, COC).

**Saving the Reaction Dictionary**  
The Ps_rxnL dictionary is serialized and saved to a pickle file (ps_rxn.pkl) for later use. 

## Correspondence Table of Numbers and Content  

The alocated number of polymerization reaction are key of Ps_rxnL.  
The numbers within the square brackets in the monomer class correspond to the values of mon_dic and key of mon_dic_inv.  

**Table 4-5. Assigned polymerization reactions and corresponding monomer classes.**   

| No. of Polymerization Reaction | Monomer Class 1        | Monomer Class 2         | Reaction Type                  | Product |
|--------------------------------|------------------------|------------------------|---------------------------------|---------|
| 1                              | Vinyl \[1\]            |                         | Addition Chain Polymerization  | Homopolymer          |
| 3                              | Cyclic Olefin \[3\]    |                         | Addition Chain Polymerization  | Homopolymer           |
| 101                            | Vinyl \[1\]            | Vinyl \[1\]             | Addition Chain Polymerization  | Alternating copolymer |
| 102                            | Vinyl \[1\]            | Cyclic Olefin \[3\]     | Addition Chain Polymerization  | Alternating copolymer |
| 103                            | Cyclic Olefin \[3\]    | Cyclic Olefin \[3\]     | Addition Chain Polymerization  | Alternating copolymer |
| 4                              | Lactone \[4\]          |                         | Ring-Opening Chain Polymerization | Homopolymer        |
| 6                              | Hydroxy Carboxylic Acid \[6\] |                  | Polycondensation               | Homopolymer           |
| 105                            | Hydroxy Carboxylic Acid \[6\] | Hydroxy Carboxylic Acid \[6\] | Polycondensation  | Alternating copolymer |
| 104                            | Di/Polycarboxylic Acid \[52\] | Di/Polyol \[53\] | Polycondensation               | Homopolymer           |
| 106                            | Di/Polyol \[53\]       | Carbon Monoxide \[10\] | Polycondensation                | Homopolymer           |
| 112                            | Cyclic Anhydride \[9\]   | Epoxide \[2\]        | Ring-Opening Chain Polymerization | Homopolymer         |
| 2                              | Epoxide \[2\]          |                      | Ring-Opening Chain Polymerization | Homopolymer         |
| 8                              | Hindered Phenol \[8\]  |                      | Polycondensation                  | Homopolymer         |
| 114                            | Bis(p-Halogenated Aryl)Sulfone \[12\] | Di/Polyol (Without Thiol) \[58\] | Polycondensation | Homopolymer            |
| 115                            | Bis(p-Fluoroaryl)Ketone \[13\] | Di/Polyol (Without Thiol) \[58\] | Polycondensation  | Homopolymer     |
| 5                              | Lactam \[5\]           |                      | Ring-Opening Chain Polymerization | Homopolymer         |
| 7                              | Amino Acid \[7\]       |                      | Polycondensation               | Homopolymer            |
| 109                            | Amino Acid \[7\]       | Amino Acid \[7\]     | Polycondensation               | Alternating copolymer  |
| 108                            | Di/Polycarboxylic Acid \[52\] | Di/Polyamine \[54\] | Polycondensation         | Homopolymer            |
| 110                            | Di/Polycyclic Anhydride \[56\] | Primary Di/Polyamine \[57\] | Polycondensation | Homopolymer           |
| 111                            | Di/Polyisocyanate \[55\] | Di/Polyol \[53\] | Polyaddition                     | Homopolymer            |
| 113                            | Di/Polyepoxide \[51\]  | Di/Polyisocyanate \[55\] | Polyaddition               | Homopolymer            |


**Table 4-6. Post-polymerization reactions of residual polymerizable functional grous(s) for the Function "bipolym".**  

| No. of Applied Polymerization Reaction | Reaction Class | Residual Functional Group           |
|----------------------------------------|----------------|-------------------------------------|
| 1, 3, 101, 102, 103                    | 200            | olefin                              |
| 112                                    | 201            | epoxide                             |
| 4, 5, 6, 7, 104, 105, 108, 109         | 202            | carboxylic acid and acyl halide     |
| 4, 5, 6, 104, 105, 111, 114, 115       | 203            | hydroxyl in alcohol and phenol      |
| 4, 5, 7                                | 204            | amine                               |
| 111                                    | 205            | isocyanate                          |
| 112                                    | 206            | carboxylic acid anhydride           |
| 110                                    | 207            | carboxylic acid anhydride           |
| 113                                    | 208            | isocyanate                          |


**Table 4-7. Assigned polymerization reactions for detail classified olefinic monomer.**

| Applied Polymerization Reaction of Olefinic Monomer | Reaction Class | Monomer Class / Residual Functional Group  |
|-----------------------------------------------------|---------------|---------------------------------------------|
| ROMP                                               | 1050           | cycCH \[1050\]                              |
| COC                                                | 1051           | cycCH \[1050\] in COC                       |
| COC                                                | 1051           | aliphCH \[1052\] in COC                     |
| Detail Classified Polyolefin Containing Diene Monomer | 209         | isomerize 1,2-added diene to 1,4-addition   |
| ROMPH                                              | 210            | olefin hydrogenation on ROMPH               |
