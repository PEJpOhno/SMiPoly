#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# polymer generator from classfied monomers.

import os
from pathlib import Path
import numpy as np
import pandas as pd
import json
import pickle
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
from .funclib import monomer_sel_MFG, monomer_sel_PFG, seq_chain, seq_successive, genmol, gencSMI

def biplym(df, targ = None, Pmode = None, dsp_rsl = None):
    if targ == None:
        targ = ['all', ]
    if Pmode == None:
        Pmode = 'a'
    if dsp_rsl == None:
        dsp_rsl = False


    db_file = os.path.join(str(Path(__file__).resolve().parent.parent), 'rules')
    with open(os.path.join(db_file, 'mon_dic.json'), 'r') as f:
        mon_dic = json.load(f)
    with open(os.path.join(db_file, 'mon_lst.json'), 'r') as f:
        monL = json.load(f)
    with open(os.path.join(db_file, 'excl_lst.json'), 'r') as f:
        exclL = json.load(f)
    with open(os.path.join(db_file, 'ps_rxn.pkl'), 'rb') as f:
        Ps_rxnL = pickle.load(f)
    with open(os.path.join(db_file, 'ps_class.json'), 'r') as f:
        Ps_classL = json.load(f)
    with open(os.path.join(db_file, 'ps_gen.pkl'), 'rb') as f:
        Ps_GenL = pickle.load(f)


    monL = {int(k): v for k, v in monL.items()}
    exclL = {int(k): v for k, v in exclL.items()}


    #set the generated polymer class
    targL = []
    if targ==['all', ]:
        targL = Ps_classL.keys()
    else:
        targL.extend(targ)
    for x in targL:
        if x not in Ps_classL.keys():
            print('oops! no such polymer class!')
            return
        else:
            pass

    #export fixed form DataFrame as pickle for next step
    DF =df.drop('ROMol', axis=1).dropna(subset=['smip_cand_mons'])
    DF_L = ['smip_cand_mons', ]
    for mon_class in mon_dic.keys():
        DF_L.append(mon_class)
    DF = DF[DF_L]
    DF_L = DF_L[1:]
    for col_nam in DF_L:
        DF[col_nam] = DF[col_nam].replace('False', '')
        DF[col_nam] = DF[col_nam].astype('bool')


    #select mode
    if Pmode == 'r':
        from .funclib import bipolymR, homopolymR
    elif Pmode == 'a':
        from .funclib import bipolymA, homopolymA
    else:
        raise Exception('invalid mode!')

    DF_Pgen = pd.DataFrame(columns=['mon1', 'mon2', 'polym', 'polymer_class'])

    #generate polymer
    for P_class in targL:
        for P_set in Ps_GenL[str(P_class)]:
            targ_mon1 = ''
            targ_mon2 = ''
            targ_mon1 = P_set[0]
            targ_mon2 = P_set[1]
            temp1 = []
            temp2 = []
            DF10 = DF[DF[targ_mon1]]
            temp1 = list(DF10['smip_cand_mons'])
            if len(temp1) != 0:
                if targ_mon2 != 'none':
                    DF20 = DF[DF[targ_mon2]]
                    temp2 = list(DF20['smip_cand_mons'])
                    del DF10, DF20
                    if len(temp2) != 0:
                        combs=[[m1, m2] for m1 in temp1 for m2 in temp2]
                        temp11=[]
                        temp21=[]
                        for comb in combs:
                            if comb[0]!=comb[1]:
                                temp11.append(comb[0])
                                temp21.append(comb[1])
                        DF_temp = pd.DataFrame()
                        DF_temp = pd.DataFrame(data={'mon1':temp11, 'mon2':temp21}, columns=['mon1', 'mon2'])
                        DF_temp['polymer_class'] = str(P_class)
                        targ_rxn=P_set[2]
                        if Pmode == 'r':
                            DF_temp['polym'] = DF_temp.apply(lambda x: [genmol(x['mon1']), genmol(x['mon2'])], axis=1).apply(bipolymR, targ_rxn=targ_rxn, monL=monL, Ps_rxnL=Ps_rxnL, P_class=P_class)
                            DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True, copy=False)
                        elif Pmode == 'a':
                            DF_temp['polym'] = DF_temp.apply(lambda x: [genmol(x['mon1']), genmol(x['mon2'])], axis=1).apply(bipolymA, targ_rxn=targ_rxn, monL=monL, Ps_rxnL=Ps_rxnL, P_class=P_class)
                            DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True, copy=False)
                else:
                    temp2 = ['' for i in range(len(temp1))]
                    DF_temp = pd.DataFrame()
                    DF_temp = pd.DataFrame(data={'mon1':temp1, 'mon2':temp2}, columns=['mon1', 'mon2'])
                    DF_temp['polymer_class'] = str(P_class)
                    mons=monL[mon_dic[targ_mon1]]
                    excls=exclL[mon_dic[targ_mon1]]
                    if Pmode == 'r':
                        DF_temp['polym'] = DF_temp.apply(lambda x: genmol(x['mon1']), axis=1).apply(homopolymR, mons=mons, excls=excls, targ_mon1=targ_mon1, Ps_rxnL=Ps_rxnL, mon_dic=mon_dic, monL=monL)
                        DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True, copy=False)
                    elif Pmode == 'a':
                        DF_temp['polym'] = DF_temp.apply(lambda x: genmol(x['mon1']), axis=1).apply(homopolymA, mons=mons, excls=excls, targ_mon1=targ_mon1, Ps_rxnL=Ps_rxnL, mon_dic=mon_dic, monL=monL)
                        DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True, copy=False)

    num_polym_react = len(DF_Pgen)
    DF_gendP = DF_Pgen.explode('polym')
    DF_gendP = DF_gendP.reset_index(drop=True)
    DF_gendP = DF_gendP.dropna(subset=['polym'])

    #adjust DataFrame
    DF_gendP['polym'].replace('', np.nan, inplace=True)

    #drpo duplicated polymerization reaction
    DF_gendP=DF_gendP.dropna(subset=['polym'])
    DF_gendP=DF_gendP[DF_gendP['mon1']!=DF_gendP['mon2']]
    DF_gendP['reactset']=np.sort(DF_gendP.loc[:,['mon1', 'mon2']].values).tolist()
    DF_gendP['reactset']=DF_gendP['reactset'].apply(set).apply(tuple)
    DF_gendP = DF_gendP.drop_duplicates(subset=['reactset', 'polym'])
    DF_gendP = DF_gendP.reset_index(drop=True)
    if dsp_rsl == True:
        if Pmode == 'a':
            print('run at advanced mode')
        elif Pmode == 'r':
            print('run at rapid mode')
        else:
            print('invalid mode')
        print('number of polymerization reactions = ', num_polym_react)
        print('number of generated polymers = ', len(DF_gendP))
    else:
        pass
    return DF_gendP
# #end
