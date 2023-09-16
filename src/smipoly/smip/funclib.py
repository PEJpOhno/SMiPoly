#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# 08/02/2021, M. Ohno
# functions for MonomerClassifier and PolymerGenerator.

import numpy as np
import pandas as pd
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem

def genmol(s):
    try:
        m = Chem.MolFromSmiles(s)
    except:
        m = np.nan
    return m


def gencSMI(m):
    try:
        cS = Chem.MolToSmiles(m)
    except:
        cS = np.nan
    return cS


#classify candidate compounds for mono-FG monomer
def monomer_sel_MFG(m, mons, excls):
    if pd.notna(m):
        chk = []
        if len(mons)!=0:
            for mon in mons:
                patt = Chem.MolFromSmarts(mon)
                if m.HasSubstructMatch(patt):
                    chk_excl=[]
                    for excl in excls:
                        excl_patt=Chem.MolFromSmarts(excl)
                        if m.HasSubstructMatch(excl_patt):
                            chk_excl.append(False)
                        else:
                            chk_excl.append(True)
                    if False in chk_excl:
                        chk.append(False)
                    else:
                        chk.append(True)
                else:
                    chk.append(False)
            if True in chk:
                fchk = True
            else:
                fchk = False
        else:
            fchk = False
    else:
        fchk = False
    return fchk


#classify candidate compounds for poly-FG monomer
#count objective FGs
def monomer_sel_PFG(m, mons, excls, minFG, maxFG):
    if pd.notna(m):
        chk_c = 0
        fchk_c = 0
        if len(mons)!=0:
            for mon in mons:
                patt = Chem.MolFromSmarts(mon)
                chk_c = len(m.GetSubstructMatches(patt))
                fchk_c = fchk_c + chk_c
            if minFG <= fchk_c <= maxFG:
                chk=[]
                for excl in excls:
                    excl_patt=Chem.MolFromSmarts(excl)
                    if m.HasSubstructMatch(excl_patt):
                        chk.append(False)
                    else:
                        chk.append(True)
                if False in chk:
                    fchk = False
                else:
                    fchk=True
            else:
                fchk = False
        else:
            fchk = (False)
    else:
        fchk = (False)
    return fchk


#define sequential polymerization for chain polymerization except polyolefine
def seq_chain(prod_P, targ_mon1, Ps_rxnL, mon_dic, monL):
    if Chem.MolToSmiles(prod_P) != '':
        if targ_mon1 not in ['vinyl', 'cOle']:
            seqFG2=Chem.MolFromSmarts(monL[[202][0]])
            seqFG3=Chem.MolFromSmarts(monL[[203][0]])
            seqFG4=Chem.MolFromSmarts(monL[[204][0]])
            while prod_P.HasSubstructMatch(seqFG2):
                prods = Ps_rxnL[202].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG3):
                prods = Ps_rxnL[203].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG4):
                prods = Ps_rxnL[204].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
        else:
            prod_P=prod_P
    return prod_P


#define sequential polymerization for successive polymerization
def seq_successive(prod_P, targ_rxn, monL, Ps_rxnL, P_class):
    if Chem.MolToSmiles(prod_P) != '':
        seqFG0=Chem.MolFromSmarts(monL[[200][0]])
        seqFG1=Chem.MolFromSmarts(monL[[201][0]])
        seqFG2=Chem.MolFromSmarts(monL[[202][0]])
        seqFG3=Chem.MolFromSmarts(monL[[203][0]])
        seqFG4=Chem.MolFromSmarts(monL[[204][0]])
        seqFG5=Chem.MolFromSmarts(monL[[205][0]])
        seqFG6=Chem.MolFromSmarts(monL[[206][0]])
        if P_class not in ['polyolefin', 'polyoxazolidone', ]:
            while prod_P.HasSubstructMatch(seqFG1):
                prods = Ps_rxnL[201].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG2):
                prods = Ps_rxnL[202].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG3):
                prods = Ps_rxnL[203].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG4):
                prods = Ps_rxnL[204].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG5):
                prods = Ps_rxnL[205].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG6):
                if P_class =='polyimide':
                    prods = Ps_rxnL[207].RunReactants([prod_P])
                elif P_class =='polyester':
                    prods = Ps_rxnL[206].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
        elif P_class in ['polyoxazolidone', ]:
            while prod_P.HasSubstructMatch(seqFG1):
                prods = Ps_rxnL[201].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
            while prod_P.HasSubstructMatch(seqFG5):
                prods = Ps_rxnL[208].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
        elif P_class in ['polyolefin', ]:
            while prod_P.HasSubstructMatch(seqFG0):
                prods = Ps_rxnL[200].RunReactants([prod_P])
                prod_P = prods[0][0]
                Chem.SanitizeMol(prod_P)
        else:
            prod_P=prod_P
    return prod_P


#homopolymerization
def homopolymR(mon1,mons,excls, targ_mon1, Ps_rxnL, mon_dic, monL):
    prod_P=mon1
    while monomer_sel_MFG(prod_P, mons, excls)== True: #生成したポリマーがさらに重合可能な場合、再度反応
        prods = Ps_rxnL[mon_dic[targ_mon1]].RunReactants([prod_P])
        try:
            prod_P = prods[0][0]
            Chem.SanitizeMol(prod_P)
            prod_P = seq_chain(prod_P, targ_mon1=targ_mon1, Ps_rxnL=Ps_rxnL, mon_dic=mon_dic, monL=monL)
        except:
            pass
    return [gencSMI(prod_P)] #20230904 revised returned Molobject to SMILES


#binarypolymerization
def bipolymR(reactant, targ_rxn, monL, Ps_rxnL, P_class):
    prod_P = Chem.MolFromSmiles('')
    prods = targ_rxn.RunReactants(reactant)
    try:
        prod_P = prods[0][0]
        Chem.SanitizeMol(prod_P)
        prod_P = seq_successive(prod_P, targ_rxn=targ_rxn, monL=monL, Ps_rxnL=Ps_rxnL, P_class=P_class)
    except:
        pass
    return [gencSMI(prod_P)] #20230904 revised returned Molobject to SMILES


#homopolymerization
def homopolymA(mon1,mons,excls, targ_mon1, Ps_rxnL, mon_dic, monL):
    prod_P=mon1
    while monomer_sel_MFG(prod_P, mons, excls)== True: #生成したポリマーがさらに重合可能な場合、再度反応
        prods = Ps_rxnL[mon_dic[targ_mon1]].RunReactants([prod_P])
        prod_Ps = []
        for prod_P in prods:
            try:
                Chem.SanitizeMol(prod_P[0])
                prod_P = prod_P[0]
                prod_P = seq_chain(prod_P, targ_mon1=targ_mon1, Ps_rxnL=Ps_rxnL, mon_dic=mon_dic, monL=monL)
                prod_Ps.append(prod_P)
            except:
                pass
    return [gencSMI(m) for m in prod_Ps]


#binarypolymerization
def bipolymA(reactant, targ_rxn, monL, Ps_rxnL, P_class):
    prod_P = Chem.MolFromSmiles('')
    prods = targ_rxn.RunReactants(reactant)
    prod_Ps = []
    for prod_P in prods:
        try:
            Chem.SanitizeMol(prod_P[0])
            prod_P = prod_P[0]
            prod_P = seq_successive(prod_P, targ_rxn=targ_rxn, monL=monL, Ps_rxnL=Ps_rxnL, P_class=P_class)
            prod_Ps.append(prod_P)
        except:
            pass
    return [gencSMI(m) for m in prod_Ps]
# end
