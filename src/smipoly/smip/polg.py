#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# polymer generator from classfied monomers.

import os
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import json
import pickle
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem
#from .funclib import monomer_sel_mfg, monomer_sel_pfg, seq_chain, seq_successive
from .funclib import genmol, coord_polym
import warnings #for warning 


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
with open(os.path.join(db_file, 'ps_class.json'), 'r') as f:
    Ps_classL = json.load(f)
with open(os.path.join(db_file, 'ps_gen.pkl'), 'rb') as f:
    Ps_GenL = pickle.load(f)

monLg = {int(k): v for k, v in monL.items()}
exclLg = {int(k): v for k, v in exclL.items()}
mon_dic_inv = {int(k): v for k, v in mon_dic_inv.items()}

def biplym(df, targ = None, Pmode = None, dsp_rsl = None):
    if targ == None:
        targ = ['all', ]
    if Pmode == None:
        Pmode = 'a'
    if dsp_rsl == None:
        dsp_rsl = False

    #FOR FUTURE WORKS!! temporary reduced the dictionaly on 04/21/2004  
    monL = {k: v for k, v in monLg.items() if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}
    exclL = {k: v for k, v in exclLg.items() if k in mon_vals[0]+mon_vals[1]+mon_vals[2]}

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

    #treat source DataFrame
    
    if 'ROMol' in df.columns: #2024/01 modified
        DF =df.drop('ROMol', axis=1).dropna(subset=['smip_cand_mons'])
    else:
        DF =df.dropna(subset=['smip_cand_mons'])
   
    DF_L = ['smip_cand_mons', ]
    for col_nam in DF.columns.values:
        if col_nam in list(mon_dic):
            DF_L.append(col_nam)
        else:
            pass
    DF = DF[DF_L]
    DF_L = DF_L[1:]
    for col_nam in DF_L:
        DF[col_nam] = DF[col_nam].replace('False', '')
        DF[col_nam] = DF[col_nam].astype('bool')

    #select mode
    if Pmode == 'r':
        warnings.warn("Mode 'r' will be deprecated and merged into mode 'a'. Use mode 'a' instead. ")
        from .funclib import bipolymR, homopolymR
    elif Pmode == 'a':
        from .funclib import bipolymA, homopolymA
    else:
        raise Exception('invalid mode!')

    DF_Pgen = pd.DataFrame(columns=['mon1', 'mon2', 'polym', 'polymer_class', 'Ps_rxnL']) #20240826added

    #generate polymer
    for P_class in targL:
        for P_set in Ps_GenL[str(P_class)]:
            targ_mon1 = ''
            targ_mon2 = ''
            targ_mon1 = P_set[0]
            targ_mon2 = P_set[1]
            temp1 = []
            temp2 = []

            #20240826 addded
            Ps_rxnL_key = []
            Ps_rxnL_key = [k for k,v in Ps_rxnL.items() if AllChem.ReactionToSmarts(v)==AllChem.ReactionToSmarts(P_set[2])]

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
                        DF_temp['Ps_rxnL'] = int(Ps_rxnL_key[0]) #20240826added
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
                    DF_temp['Ps_rxnL'] = int(Ps_rxnL_key[0]) #20240826added
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
    DF_gendP.replace({'polym':{'':np.nan}}, inplace=True)

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

#set the olefin class(es) of the copolymer
def ole_copolym(df, targ = None, ncomp = None, dsp_rsl = None, drop_dupl = None):
    if targ == None or len(targ)==0:
        print('Plz define the olefin class(es)')
        return
    if ncomp == None:
        ncomp = 1
    if dsp_rsl == None:
        dsp_rsl = False
    if drop_dupl == None:
        drop_dupl = True

    #explanation
    template_ole_keys = [mon_dic_inv[i] for i in mon_vals[3]]
    print('\n', 'valid arguments \n', template_ole_keys , 
    ',\n', Ps_GenL['rec:coord'], '\n')

    #confirm the list of copolymerization unit(s)
    if type(targ) != list:
        print('This arg. should be given as list')
        return
    for x in targ:
        if x not in [mon_dic_inv[k] for k in mon_vals[3]]+list(Ps_GenL['rec:coord']):
            print('oops! no such olefin class!')
            return
        else:
            pass
    if len(targ)>ncomp:
        print('reconfirm ncomp')
        return
    if any(e in targ for e in Ps_GenL['rec:coord']):
        if len(targ) > 1:
            print('ROMP, ROMPH, and COC should each be used solely. ')
            return
        if ncomp > 2:
            print('Reccomend: nocmp=1 for ROMP, ROMPH and 2 for COC')

    #reconsruct the list of  cllasified olefin monomers for co-polymerization
    cand_cru = []
    comb_cru = []
    copoly_cru = []

    ole_clsL = df['ole_cls'].to_list()
    cand_monsL = df['smip_cand_mons'].to_list()
    ole_targL = [e for e in list(zip(ole_clsL, cand_monsL)) if True in [d[0] for d in e[0].values()]] #[({オレフィン種:[該非, 官能基数, CRU], }, 出発モノマー), ()]
    ole_targL2 = []
    if all([item in template_ole_keys for item in targ]):
        for ole in [mon_dic_inv[k] for k in mon_vals[3]]:
            for m in ole_targL:
                if m[0][ole][0]==True:
                    ole_targL2.append([ole, m[0][ole][2], m[1]])

        #Extract CRUs corresponding to defined components
        cand_cru = list(itertools.chain.from_iterable([[m for m in ole_targL2 if t==m[0]] for t in targ]))
    
    elif targ==['ROMP',] or targ==['ROMPH',]:
        ole = 'cycCH'
        for m in ole_targL:
            if m[0][ole][0]==True:
                ole_targL2.append([ole, m[0][ole][2], m[1]])
        for e in ole_targL2:
            e[1] = coord_polym(e[2], Ps_rxnL[1050])
        #explode the list of generated CRU
        cand_cru = [[e[0], sub_e, e[2]] for e in ole_targL2 for sub_e in e[1]]
        if targ==['ROMPH']:
            for e in cand_cru:
                e[1] = e[1].replace("=", "")

    elif targ==['COC',]:
        targ.append('') #dummy element to count components
        if ncomp<2:
            print('reconfirm ncomp')
            return
        else:
            coc_mons = ['cycCH', 'aliphCH']
            ole_targL2cyc = []
            ole_targL2chain = []
            for ole in coc_mons:
                if ole=='cycCH':
                    for m in ole_targL:
                        if m[0][ole][0]==True:
                            if ole=='cycCH':
                                ole_targL2cyc.append([ole, m[0][ole][2], m[1]])
                                for e in ole_targL2cyc:
                                    e[1] = coord_polym(e[2], Ps_rxnL[1051])
                elif ole=='aliphCH':
                    for m in ole_targL:
                        if m[0][ole][0]==True:
                            if ole=='cycCH':
                                ole_targL2chain.append([ole, m[0][ole][2], m[1]])
                                for e in ole_targL2:
                                    e[1] = coord_polym(e[2], Ps_rxnL[1052])
                            ole_targL2chain.append([ole, m[0][ole][2], m[1]])
                            for e in ole_targL2chain:
                                e[1] = coord_polym(e[2], Ps_rxnL[1052])
            ole_targL2 = list(itertools.chain(ole_targL2cyc, ole_targL2chain))
            #explode the list of generated CRU
            cand_cru = [[e[0], sub_e, e[2]] for e in ole_targL2 for sub_e in e[1]]

    else:
        print('reconfirm olefin class(es)')
        return


    #generate copolymers
    comb_cru = [e for e in itertools.combinations(cand_cru, ncomp)]

    #Exclude if it does not contain all defined olefin class as the component
    for e in comb_cru:
        if len(set([l[0] for l in e]))>=len(targ):
            copoly_cru.append(([l[0] for l in e], [l[1] for l in e], [l[2] for l in e]))
        else:
            pass

    #Export as Pandas DataFrame
    DF_Pgen = pd.DataFrame(columns=['mon1', 'mon2', 'polym', 'polymer_class', 'Ps_rxnL', 'reactset'])
    DF_Pgen[['polymer_class', 'polym', 'reactset']] = pd.DataFrame(copoly_cru)
    DF_Pgen = DF_Pgen.fillna({'mon1':'', 'mon2':''})

    #Type of initiator  
    l_initiator = ['rec:radi', 'rec:cati', 'rec:ani']
    rec_initiator = []
    if any(e in targ for e in Ps_GenL['rec:coord']):
        rec_initiator = targ
    for k in l_initiator:
        if all([e in Ps_GenL[k] for e in targ]):
            rec_initiator.append(k)
        else:
            pass
    if len(rec_initiator)==0:
        rec_initiator=np.nan
    DF_Pgen['Ps_rxnL'] = DF_Pgen.apply(lambda x: rec_initiator, axis=1)

    #Drop duplicated copolymer. It takes long time when the DataFrame (DF_Pgen) was large.
    if drop_dupl == True:
        DF_Pgen['temp_dro_dupl'] = DF_Pgen.apply(lambda x:"".join(sorted(x['polym'])), axis=1)
        DF_gendP = DF_Pgen.drop_duplicates(subset='temp_dro_dupl')
        DF_gendP = DF_gendP.drop('temp_dro_dupl', axis=1)
        DF_gendP = DF_gendP.reset_index(drop=True)
    else:
        DF_gendP = DF_Pgen

    if dsp_rsl==True:
        print('Number of generated (co)polymer, ',  ncomp, ' component(s) system : ', format(len(DF_gendP), ','))
    else:
        pass

    return DF_gendP

# #end
