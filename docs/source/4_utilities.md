# 4. utilities; Overview

.ipynb files in directiry "utirities" are tools to imprementing the rules of monomer extraction-cllasification , polymerization reaction definition and polymer generation.  

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
