# 2. Installing SMiPoly  
SMiPoly can be installed with pip or conda. 
## 2-1. Install with pip  
Create new virtual environment and activate it.
To install this package, run as follows.

```sh
$pip install smipoly
```
## 2-2. Install with conda  
Add the channel "conda-forge" if it have not been enable.  

```sh
$conda config --add channels conda-forge
```

Create a new environment. 
```sh
$conda create -n "YOUR_NEW_ENVIRONMNT_NAME" python  
or 
$conda create -n "YOUR_NEW_ENVIRONMNT_NAME" python="required version (ex. 3.10)"
```
Then activate it. 
```sh
$conda activate "YOUR_NEW_ENVIRONMNT_NAME"
```
And install SMiPoly. 
```sh
$conda install smipoly
```

Or after create and activate a new environment, 
```sh
$conda install conda-forge::smipoly
```  
