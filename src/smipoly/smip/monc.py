#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# 07/27/2021, M. Ohno
# smilesで与えた化合物リストを、モノマー別に分類する
# monomer categolization system of the compound list in SMILES
#
# Refernce:
# https://future-chem.com/rdkit-chemical-rxn/
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
# https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

import os
from pathlib import Path
import numpy as np
import pandas as pd
import json
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
from .funclib import genmol, genc_smi, monomer_sel_mfg, monomer_sel_pfg

db_file = os.path.join(str(Path(__file__).resolve().parent.parent), 'rules')
with open(os.path.join(db_file, 'mon_vals.json'), 'r') as f:
    mon_vals = json.load(f) 
with open(os.path.join(db_file, 'mon_dic.json'), 'r') as f:
    mon_dic = json.load(f)
with open(os.path.join(db_file, 'mon_dic_inv.json'), 'r') as f:
    mon_dic_inv = json.load(f)
with open(os.path.join(db_file, 'mon_lst.json'), 'r') as f:
    monL = json.load(f)
with open(os.path.join(db_file, 'excl_lst.json'), 'r') as f:
    exclL = json.load(f)

monLg = {int(k): v for k, v in monL.items()}
exclLg = {int(k): v for k, v in exclL.items()}
mon_dic_inv = {int(k): v for k, v in mon_dic_inv.items()}


def moncls(df, smiColn, minFG=None, maxFG=None, dsp_rsl=None):
    # #The default number of the samle class of FG were limited 2 to 4 in the same molecule for poly functionalized monomer.
    # (dataframe, smiColn, minFG = 2, maxFG = 4)
    if minFG == None:
        minFG = 2
    if maxFG == None:
        maxFG = 4
    if dsp_rsl == None:
        dsp_rsl = False

    monL = {k: v for k, v in monLg.items() if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}
    exclL = {k: v for k, v in exclLg.items() if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}

    #read source file
    DF01 = df
    smiColn = smiColn

    #append CO and HCHO for carbonate
    DF02 = DF01.copy()
    if smiColn in DF02.columns.to_list():
        #DF02 = DF02[['CID', smiColn]] #remove this row if not required
        DFadd = pd.DataFrame([['[C-]#[O+]'], ['C=O'], ], columns=[smiColn])
        DF02 = pd.concat([DF02, DFadd], ignore_index=True)
    else:
        print("invalid SMILES column name")


    #drop NA of smiles, and add chemical structure
    DF02['ROMol'] = DF02[smiColn].apply(genmol)
    DF02['smip_cand_mons'] = DF02['ROMol'].apply(genc_smi)

    #classification for mono-functionalized monomer
    #count functional groupe, remove exclude compounds andjudge targetted monomer or not.
    for i in mon_vals[0]: 
        mons=()
        excls=()
        mons = monL[i]
        excls = list(exclL[i])
        DF02[mon_dic_inv[i]] = [e[0] for e in DF02['ROMol'].apply(monomer_sel_mfg, mons=mons, excls=excls)]
        if dsp_rsl==True:
            print(i)
            print(mon_dic_inv[i], ' = ', len(DF02[DF02[mon_dic_inv[i]]==True]), ' / ', len(DF02))
        else:
            pass


    #classification for poly-functionalized monomer
    for i in mon_vals[1]:
        mons=()
        excls=()
        mons=monL[i]
        excls=exclL[i]
        DF02[mon_dic_inv[i]] = [e[0] for e in DF02['ROMol'].apply(monomer_sel_pfg, mons=mons, excls=excls, minFG=minFG, maxFG=maxFG)]
        if dsp_rsl==True:
            print(i)
            print(mon_dic_inv[i], ' = ', len(DF02[DF02[mon_dic_inv[i]]==True]), ' / ', len(DF02))
        else:
            pass

    DF02 = DF02.drop('ROMol', axis=1) #2024/01 modified
    return DF02


    #classification for olefinic monomer **UNDER THE CONSTRUCTION**
    #count functional groupe, remove exclude compounds andjudge targetted monomer or not.
def olecls(df, smiColn, minFG=None, maxFG=None, dsp_rsl=None):
    if minFG == None:
        minFG = 1
    if maxFG == None:
        maxFG = 4
    if dsp_rsl == None:
        dsp_rsl = False

    monL = {k: v for k, v in monLg.items() if k in mon_vals[3]}
    exclL = {k: v for k, v in exclLg.items() if k in mon_vals[3]}

    #read source file
    DF02 = df
    smiColn = smiColn
    #drop NA of smiles, and add chemical structure
    DF02['ROMol'] = DF02[smiColn].apply(genmol)
    DF02['smip_cand_mons'] = DF02['ROMol'].apply(genc_smi)

    for i in mon_vals[3]: 
        mons=()
        excls=()
        mons = monL[i]
        excls = list(exclL[i])
        DF02[mon_dic_inv[i]] = [e[0] for e in DF02['ROMol'].apply(monomer_sel_pfg, mons=mons, excls=excls, minFG=minFG, maxFG=maxFG)]
        if dsp_rsl==True:
            print(mon_dic_inv[i], ' = ', len(DF02[DF02[mon_dic_inv[i]]==True]), ' / ', len(DF02))
        else:
            pass

    DF02 = DF02.drop('ROMol', axis=1)
    return DF02
#end
