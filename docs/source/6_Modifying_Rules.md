# 6. Modifying rules  

## Common Operations
### Pre-Rule Change Operations (Preparation)  
1. Clone the [SMiPoly GitHub repository](https://github.com/PEJpOhno/SMiPoly.git) and note the path on your local machine (e.g., [path_to_cloned_repo]).
2. Delete the 8 files within the ```utilities/rules``` directory of the cloned repository, leaving an empty directory.

Modify the rules using the notebooks in the ```utility``` directory: ```1_MonomerDefiner.ipynb```, ```2_Ps_rxnL.ipynb```, and ```3_Ps_GenL.ipynb```.  

### Rule Change Operations (Applying the Changed Rules)  
1. Run the following files in the utilities directory in this specific order: ```1_MonomerDefiner.ipynb```, ```2_Ps_rxnL.ipynb```, and ```3_Ps_GenL.ipynb```. Make sure to run all of them, regardless of whether any changes have been made to the notebooks.  
2. Copy the 8 files from ```utilities/rules```.  
3. Delete the 8 files from ```src/smipoly/rules```.  
4. Paste the 8 files you copied in step 2 into the now-empty ```src/smipoly/rules directory```.  
5. Create and activate a virtual environment to run the modified SMiPoly.  
6. Install with the following command:  

```
$ pip install [path_to_cloned_repo]
```
7. Verify the operation. If there are no issues, you can delete the cloned repository.  


## Modifying Rules  
Monomer and CRU substructures were represented using [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html), while the polymerization reactions were described using [reaction SMARTS](https://daylight.com/meetings/summerschool00/course/basics/smirks.html).  

### Modifying an Existing Monomer Class  
Open ```1_MonomerDefiner.ipynb``` and add or remove the corresponding monomer definitions in ```monL``` and/or ```exclL```.  
Since certain polymers are associated with multiple polymerization reactions, care should be taken regarding their impact scope.  

**Example 1.**
- Exclude vinyl monomers without F at the terminal methylidene of the polymerization site.  　　
```
# origine
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', '[CX3](F)(F)=[CX3]', '[CX3;H1](F)=[CX3]',
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', '[CX3](Cl)(F)=[CX3]')

# update
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', 
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', )
```
**Example 2.**  
- Add terminal dibromoolefins.  
```
# origine
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', '[CX3](F)(F)=[CX3]', '[CX3;H1](F)=[CX3]',
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', '[CX3](Cl)(F)=[CX3]')

# update
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', '[CX3](F)(F)=[CX3]', '[CX3;H1](F)=[CX3]',
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', '[CX3](Cl)(F)=[CX3]',
         '[CX3]([Br])([Br])=[CX3]',)
```
**Example 3.**  
- Exclude compounds containing fluorine (F) in the molecule using ```exclL[n]```.　　
```
# origine
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', '[CX3](F)(F)=[CX3]', '[CX3;H1](F)=[CX3]',
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', '[CX3](Cl)(F)=[CX3]')
exclL[n]=('[CX3H1]=[O]', '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]',
         '[CX3;!R](=[OX1])[OX2;!R][CX3;!R](=[OX1])', '[CX3;!R](=[OX1])[NX3;!R][CX3;!R](=[OX1])')

# update
n=mon_dic['vinyl']
monL[n]=('[CX3H2]=[CX3]', 
         '[CX3](Cl)(Cl)=[CX3]', '[CX3;H1](Cl)=[CX3]', )
exclL[n]=('[F]',
         '[CX3H1]=[O]', '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]',
         '[CX3;!R](=[OX1])[OX2;!R][CX3;!R](=[OX1])', '[CX3;!R](=[OX1])[NX3;!R][CX3;!R](=[OX1])')
```

### Defining a New Monomer Class  
Open ```1_MonomerDefiner.ipynb```.  
Define aziridine as a new monomer class. For simplicity, the definition of aziridine is as follows:  

- A three-membered ring containing one N and two C atoms.
- In aziridine, the ring carbons are must be CH₂ group.    
- The substituent on the nitrogen is a hydrogen atom, alkylcarbonyl or arylsulfonyl.
- The molecule does not contain any acyclic primary or secondary amines.
- Aziridines undergoe ring-opening polymerization to give polyamines.

1. Add its name as akey to ```mon_dic```. For this case, assign an unused integer below 50 as the value. The example uses 14 as the value.  
```
# origine
mon_dic = {"vinyl":1, "epo":2, "diepo":51, "cOle":3, "lactone":4, "lactam":5, "hydCOOH":6, "aminCOOH":7,
           "hindPhenol":8, "cAnhyd":9,  "CO":10, "HCHO":11, "sfonediX":12, "BzodiF":13,
           "diCOOH":52, "diol":53, "diamin":54, "diNCO":55, "dicAnhyd":56, "pridiamin":57, "diol_b":58,
           "acryl":1001, "bEWole":1002, "styryl":1003, "allyl":1004, "haloCH":1005, "vinylester":1006,
           "malei":1007, "conjdiene":1020, "vinylether":1030, "tertcatCH":1031, "cycCH":1050, "aliphCH":1052, }

# update
mon_dic = {"vinyl":1, "epo":2, "diepo":51, "cOle":3, "lactone":4, "lactam":5, "hydCOOH":6, "aminCOOH":7,
           "hindPhenol":8, "cAnhyd":9,  "CO":10, "HCHO":11, "sfonediX":12, "BzodiF":13,
           "aziridine":14,
           "diCOOH":52, "diol":53, "diamin":54, "diNCO":55, "dicAnhyd":56, "pridiamin":57, "diol_b":58,
           "acryl":1001, "bEWole":1002, "styryl":1003, "allyl":1004, "haloCH":1005, "vinylester":1006,
           "malei":1007, "conjdiene":1020, "vinylether":1030, "tertcatCH":1031, "cycCH":1050, "aliphCH":1052, }
```

2. The numerical values used as keys in ```mon_dic``` should be registered in the corresponding tuple within ```mon_vals```. Aziridine has been added to ```mon_vals[0]```.
```
# origine
mon_vals = ((1,2,3,4,5,6,7,8,9,10,11,12,13), (51,52,53,54,55,56,57,58), (200, 201, 202, 203, 204, 205, 206),
            (1001, 1002, 1003, 1004, 1005, 1006, 1007,1020, 1030, 1031, 1050, 1052))

# update
mon_vals = ((1,2,3,4,5,6,7,8,9,10,11,12,13, 14),
            (51,52,53,54,55,56,57,58), (200, 201, 202, 203, 204, 205, 206),
            (1001, 1002, 1003, 1004, 1005, 1006, 1007,1020, 1030, 1031, 1050, 1052))
```
3. Set up monL and exclL entries for aziridine.  
```
n=mon_dic['aziridine']
monL[n]=('[CX4H2]1[NX3H1][CX4H2]1', '[CX4H2]1[NX3]([CX3](=[OX1])[CX4])[CX4H2]1',
         '[CX4H2]1[NX3]([SX4](=[OX1])(=[OX1])[c])[CX4H2]1')
exclL[n]=('[N&X3;H2,H1;!$(N[C,S]=*);!R]', )
```

### Defining a New Polymerization Reaction  
Open ```2_Ps_rxnL.ipynb```.  
For polymerization involving a single monomer, use the value associated with the monomer class as the key in Ps_rxnL, taken from mon_dic. In cases where multiple monomers are involved or other polymerization reactions have already been defined, use an undefined integer as the key (see the relevant section in 2_Ps_rxnL.ipynb).  

```
n= mon_dic['aziridine']
aziridine_homo = '[CX4;H2;R:1]1[NX3;R:2][CX4;H2;R:3]1>>*-[CX4;H2:1][CX4;H2:3][NX3:2]-*'
Ps_rxnL[n] = AllChem.ReactionFromSmarts(aziridine_homo)
```

### Defining a New Polymer Class  
Open ```3_Ps_GenL.ipynb```.  
When the product of a new polymerization reaction does not fit into an existing polymer class, a suitable polymer class should be created.　　
For aziridine polymerization, update dictionary of polymerclass:```Ps_classL```. Generally, the value corresponding to a key in the polymer class should follow the Polymer Class No. in the [PoLyInfo database (https://polymer.nims.go.jp/PoLyInfo/guide/jp/term_polymer.html#chap06)](https://polymer.nims.go.jp/PoLyInfo/guide/jp/term_polymer.html#chap06).  

```
# origine
Ps_classL = {'polyolefin':11, 'polyester':6, 'polyether':12, 'polyamide':2, 'polyimide':8, 'polyurethane':19, 
            'polyoxazolidone':23, }
# update
Ps_classL = {'polyolefin':11, 'polyester':6, 'polyether':12, 'polyamide':2, 'polyimide':8, 'polyurethane':19, 
            'polyoxazolidone':23, 'polyamine':9}
```

Then define the combination of polymer class, applied monomer(s) and polymerization reaction.  
```
Ps_GenL['polyamine'] = (('aziridine', 'none', Ps_rxnL[mon_dic['aziridine']]), )
```
Where, 
'polyamine': newly defined polymer class  
'aziridine': the key of the monomer in mon_dic
'none': For polymerization involving a single monomer, set the second element to 'none'.
Ps_rxnL[mon_dic['aziridine']]): newly defined polymerization reaction. If an integer is specified as the key in "Defining a New Polymerization Reaction", use that value as the key for ```Ps_rxnL[]```.

If the resulting polymer belongs to an existing polymer class, append a tuple to that class consisting of:
- the keys of all monomers used in the newly defined polymerization reaction (from mon_dic), and
- the newly defined Ps_rxnL[n].  

After all modifications have been completed, run "Rule Change Operations".  

## TIPs
Points to note when using SMARTS notation, especially reaction SMARTS: 
Atom mapping should be applied to the following:
**Examples**
- Polymerization sites  
```
✔️ Correct: [CX3;H2:1]=[CX3;H1:2]  
❌ Incorrect: [CX3;H2]=[CX3;H]
```

- Atoms bonded to atoms outside the substructure  
```
✔️ Correct: [CX3;H2:1]=[CX3;H1:2][CX3](=[OX1])[OX2,SX2:3] 
❌ Incorrect: [CX3;H2:1]=[CX3;H1:2]CX3[OX2,SX2]
```

- Atoms for which multiple elements are possible candidates  
```
✔️ Correct: [CX3;H2:1]=[CX3;H0:2][F,Cl:3]
❌ Incorrect: [CX3;H2:1]=[CX3;H0:2][F,Cl]
```

