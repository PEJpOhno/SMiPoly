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
### Modifying an Existing Monomer Class  
```1_MonomerDefiner.ipynb```  


### Defining a New Monomer Class  

### Defining a New Polymerization Reaction  
```2_Ps_rxnL.ipynb```  


When the product of the new polymerization reaction is included in an existing polymer class  
```3_Ps_GenL.ipynb```

When the product of the new polymerization reaction is not included in an existing polymer class  
```3_Ps_GenL.ipynb```

