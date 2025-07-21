#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# license that can be found in the LICENSE file.

# 07/27/2021, M. Ohno

"""
Monomer categolization system of the compound list in SMILES. 

Classifies monomers based on functional groups and other criteria.
This function processes a DataFrame containing SMILES strings, 
extracts and classifies  monomers, and appends the results 
to the DataFrame. It also supports optional 
display of classification results.

"""
import os
from pathlib import Path
import json
import pickle
import numpy as np
import pandas as pd
# from rdkit import rdBase, Chem
# from rdkit.Chem import AllChem, Draw
from .funclib import (
    genmol,
    genc_smi,
    monomer_sel_mfg,
    monomer_sel_pfg,
    ole_sel_cru,
    diene_14, update_nested_dict
)

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
with open(os.path.join(db_file, 'ps_rxn.pkl'), 'rb') as f:
    Ps_rxnL = pickle.load(f)

monLg = {int(k): v for k, v in monL.items()}
exclLg = {int(k): v for k, v in exclL.items()}
mon_dic_inv = {int(k): v for k, v in mon_dic_inv.items()}


def moncls(df, smiColn, minFG=None, maxFG=None, dsp_rsl=None):
    """
    Select monomers from given dataset of small molecule compounds and 
    categolize them into a monomer class. 

    Args:
        df (pd.DataFrame): Input DataFrame containing chemical data.
            smiColn (str): Column name in the DataFrame containing 
            SMILES strings.
        minFG (int, optional): Minimum number of functional groups 
            for poly-functionalized monomers. Defaults to 2.
        maxFG (int, optional): Maximum number of functional groups 
            for poly-functionalized monomers. Defaults to 4.
        dsp_rsl (bool, optional): Whether to display classification 
            results. Defaults to False.

    Returns:
        pd.DataFrame: A modified DataFrame with classification 
        results appended.

    Raises:
        ValueError: If the specified SMILES column name is invalid.

    Notes:
        - This function forcibly adds CO for carbonate synthesis 
          and HCHO for future works."

    """
    # #The default number of the samle class of FG were limited
    # 2 to 4 in the same molecule for poly functionalized monomer.
    # (dataframe, smiColn, minFG = 2, maxFG = 4)
    if minFG is None:
        minFG = 2
    if maxFG is None:
        maxFG = 4
    if dsp_rsl is None:
        dsp_rsl = False

    monL = {k: v for k, v in monLg.items(
    ) if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}
    exclL = {k: v for k, v in exclLg.items(
    ) if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}

    # read source file
    DF01 = df
    smiColn = smiColn

    # append CO and HCHO for carbonate
    DF02 = DF01.copy()
    if smiColn in DF02.columns.to_list():
        # DF02 = DF02[['CID', smiColn]] #remove this row if not required
        DFadd = pd.DataFrame([['[C-]#[O+]'], ['C=O'], ], columns=[smiColn])
        DF02 = pd.concat([DF02, DFadd], ignore_index=True)
    else:
        print("invalid SMILES column name")

    # drop NA of smiles, and add chemical structure
    DF02['ROMol'] = DF02[smiColn].apply(genmol)
    DF02['smip_cand_mons'] = DF02['ROMol'].apply(genc_smi)

    # classification for mono-functionalized monomer
    # count functional groupe, remove exclude compounds andjudge
    # targetted monomer or not.
    for i in mon_vals[0]:
        mons = ()
        excls = ()
        mons = monL[i]
        excls = list(exclL[i])
        DF02[mon_dic_inv[i]] = [e[0]
                                for e in DF02['ROMol'].apply(
                                    monomer_sel_mfg, mons=mons,
                                    excls=excls)]
        if dsp_rsl:
            print(i)
            print(mon_dic_inv[i], ' = ', len(
                DF02[DF02[mon_dic_inv[i]] == True]), ' / ', len(DF02))
        else:
            pass

    # classification for poly-functionalized monomer
    for i in mon_vals[1]:
        mons = ()
        excls = ()
        mons = monL[i]
        excls = exclL[i]
        DF02[mon_dic_inv[i]] = [
            e[0] for e in DF02['ROMol'].apply(
                monomer_sel_pfg, mons=mons, excls=excls,
                minFG=minFG, maxFG=maxFG
            )
        ]
        if dsp_rsl:
            print(i)
            print(mon_dic_inv[i], ' = ', len(
                DF02[DF02[mon_dic_inv[i]] == True]), ' / ', len(DF02))
        else:
            pass

    DF02 = DF02.drop('ROMol', axis=1)  # 2024/01 modified
    return DF02


# classification for olefinic monomer
# count functional groupe, remove exclude compounds andjudge targetted
# monomer or not.


def olecls(df, smiColn, minFG=None, maxFG=None, dsp_rsl=None):
    """
    Select olefinic monomers from given dataset of small molecule 
    compounds and categolize them into a olefinic monomer class. 

    Args:
        df (pd.DataFrame): The input DataFrame containing chemical data.
            Must include the structure of a compound written in SMILES.
        smiColn (str): The column name in the DataFrame containing 
            SMILES strings.
        minFG (int, optional): Minimum number of functional groups 
            to consider. Defaults to 1.
        maxFG (int, optional): Maximum number of functional groups 
            to consider. Defaults to 4.
        dsp_rsl (bool, optional): Whether to display results 
            during processing. Defaults to False.

    Returns:
        pd.DataFrame: The updated DataFrame with olefin classification 
        results.

    Notes:
        - The function assumes the existence of several global
          variables such as `monLg`, `exclLg`,
          `mon_vals`, `mon_dic_inv`, and `Ps_rxnL`.
        - The function modifies the input DataFrame by adding
          new columns for olefin classification.
        - The `genmol`, `genc_smi`, `ole_sel_cru`,
          `update_nested_dict` and `diene_14` functions
          are defined in 'funclib.py'.
        - The `ole_cls` column is refined for conjugated
          diene classification using a specific reaction.

    """
    if minFG is None:
        minFG = 1
    if maxFG is None:
        maxFG = 4
    if dsp_rsl is None:
        dsp_rsl = False

    monL = {k: v for k, v in monLg.items() if k in mon_vals[3]}
    exclL = {k: v for k, v in exclLg.items() if k in mon_vals[3]}

    template_ole_keys = [mon_dic_inv[i] for i in mon_vals[3]]
    print(template_ole_keys)

    # read source file
    DF02 = df
    smiColn = smiColn
    # drop NA of smiles, and add chemical structure
    DF02['ROMol'] = DF02[smiColn].apply(genmol)
    DF02['smip_cand_mons'] = DF02['ROMol'].apply(genc_smi)
    # create null column for olefin classification
    DF02['ole_cls'] = [
        {
            k: v for k, v in zip(
                template_ole_keys, [
                    np.nan for x in range(len(template_ole_keys))]
            )
        }
        for y in range(len(DF02))
    ]

    for i in mon_vals[3]:
        mons = ()
        excls = ()
        mons = monL[i]
        excls = list(exclL[i])
        DF02['temp'] = ['' for e in range(len(DF02))]
        DF02['temp'] = DF02['ROMol'].apply(
            ole_sel_cru, mons=mons, excls=excls, minFG=minFG, maxFG=maxFG)
        DF02 = DF02.apply(update_nested_dict, axis=1, args=(
            'ole_cls', 'temp', mon_dic_inv[i]))

        if dsp_rsl:
            print(mon_dic_inv[i], ' = ',
                  list(DF02['ole_cls'].apply(
                      lambda x: x[mon_dic_inv[i]][0] == True)).count(True),
                  ' / ', len(DF02))

        else:
            pass
    DF02['ole_cls'] = df['ole_cls'].apply(
        diene_14, rxn=Ps_rxnL[209])  # refine conjugated diene CRU
    DF02 = DF02.drop('ROMol', axis=1)
    DF02 = DF02.drop('temp', axis=1)
    return DF02

# end
