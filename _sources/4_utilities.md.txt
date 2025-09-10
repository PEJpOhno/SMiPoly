# 4. utilities; Rule Modification Utilities

## Overview
The .ipynb files in the utilities directory are tools for implementing rules related to monomer extraction and classification, polymerization reaction definitions, and polymer generation.
You can modify the rules using the notebooks in the ```utility``` directory: ```1_MonomerDefiner.ipynb```, ```2_Ps_rxnL.ipynb```, and ```3_Ps_GenL.ipynb```.  
These files require Jupyter Notebook to run. 

By using the files in './utilities' directory, one can modify or add the definition of monomers, the rules of polymerization reactions and polymer classes.  
To apply the new rule(s), replace the old './smipoly/rules' directory by the new one. The files must be run according to the number assigned the head of the each filename.  

  - 1_MonomerDefiner.ipynb: definitions of monomers  
  - 2_Ps_rxnL.ipynb: rules of polymerization reactions    
  - 3_Ps_GenL.ipynb: definitions of polymer classes with combinations of starting monomer(s) and polymerization reaction  

Each definations were written by [SMARTS (SMiles ARbitrary Target Specification](https://www.daylight.com/dayhtml_tutorials/languages/smarts/)) strings.  

The difined monomer classificatin, polymerization reaction and the polymer genration rules were located in src/smipoly/rulse as pickle or json file.  
These were genarated by using utilities/1_MonomerDefiner.ipynb, 2_Ps_rxnL.ipynb and 3_Ps_GenL.ipynb. 

|Generator .ipynb file name |Valiable name |Saved file name |Defined rule |
|---------------------------|--------------|----------------|-------------|
|1_MonomerDefiner.ipynb |mon_vals |mon_vals |Clategorization of monomer classes |
|                       |mon_dic, mon_dic_inv |mon_dic.json, mon_dic_inv.json |The defination of monomer |
|                       |monL |mon_lst.json |The definition of the polymerization site |
|                       |exclL |excl_lst.json |The definition of functional groups that should not coexist in a monomer molecule |
|2_Ps_rxnL.ipynb |Ps_rxnL  |ps_rxn.pkl    |polymerization reaction |
|3_Ps_GenL.ipynb |Ps_classL |ps_class.json |Polymer classes |
|               |Ps_GenL   |ps_gen.pkl    |Polymer classes with applied polymerization reactions and starting monomer(s) |



## Common Operations To Modifying Rules
### Pre-modification Operations  
1. Clone the [SMiPoly GitHub repository](https://github.com/PEJpOhno/SMiPoly.git) and note the path on your local machine (e.g., [path_to_cloned_repo]).
2. Delete the 8 files within the ```utilities/rules``` directory of the cloned repository, leaving an empty directory. 


### Post-modification Operations
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
- Exclude vinyl monomers with F at the terminal methylidene of the polymerization site.  　　
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

A detailed definition of olefin monomers (mon_dic[\<'olefinic monomer class'\>] >= 1001) is as follows:
- Position the C=C involved in polymerization as close to the beginning of the SMARTS string as possible.  
- Assigning atom mapping to the atoms involved in the reaction. This allows the polymerization reaction equation to be automatically generated by ```ole_rxnsmarts_gen(reactant)``` and ```ole_cru_gen()```.    
- In SMARTS notation, atoms are denoted by their element symbols, not by their atomic numbers.

### Defining a New Polymerization Reaction  
Open ```2_Ps_rxnL.ipynb```.  
For polymerization involving a single monomer, use the value associated with the monomer class as the key in Ps_rxnL, taken from mon_dic. In cases where multiple monomers are involved or other polymerization reactions have already been defined, use an undefined integer as the key (see the relevant section in 2_Ps_rxnL.ipynb).  

```
n= mon_dic['aziridine']
aziridine_homo = '[CX4;H2;R:1]1[NX3;R:2][CX4;H2;R:3]1>>*-[CX4;H2:1][CX4;H2:3][NX3:2]-*'
Ps_rxnL[n] = AllChem.ReactionFromSmarts(aziridine_homo)
```  
- Describe the skeleton directly involved in polymerization.  

Except in special cases, for polymerization via ```ole_copolym()```, this definition isn't required because the polymerization reaction equation is automatically generated by ```ole_cru_gen()```.


### Defining a New Polymer Class  
Open ```3_Ps_GenL.ipynb```.  
When the product of a new polymerization reaction does not fit into an existing polymer class, a suitable polymer class should be created.　　
For aziridine polymerization, update dictionary of polymerclass:```Ps_classL```. Generally, the value corresponding to a key in the polymer class should follow the Polymer Class No. in the [PoLyInfo database (https://polymer.nims.go.jp/PoLyInfo/guide/jp/term_polymer.html#chap06)](https://polymer.nims.go.jp/PoLyInfo/guide/jp/term_polymer.html#chap06).  
However, when you add a definition for a polymerization reaction using ```ole_copolym()``` for olefinic monomers that have been classified in detail by ```olecls()```, do not append to ```Ps_classL```.  

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

In the case of running polymerization via ole_copolym(), you may also add the recommended polymerization conditions for the corresponding olefinic monomer to:  
```Ps_GenL['rec:radi']``` , ```Ps_GenL['rec:cati'] ```, ```Ps_GenL['rec:ani'] ```, ```Ps_GenL['rec:coord']```


After all modifications have been completed, run "Post-modification Operations".  

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
