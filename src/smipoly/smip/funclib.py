#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2021 Mitsuru Ohno

# Use of this source code is governed by a BSD-3-style
# license that can be found in the LICENSE file.

# 08/02/2021, M. Ohno

"""
functions for MonomerClassifier (monc.py) and 
PolymerGenerator (polyg.py).

"""
import numpy as np
import pandas as pd
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem


def genmol(s):
    """
    Generates a molecular object from a SMILES string.

    Args:
        s (str): A SMILES (Simplified Molecular Input Line Entry System) string
            representing the molecular structure.

    Returns:
        rdkit.Chem.Mol or numpy.nan: A molecular object if the SMILES string 
        is valid, otherwise returns numpy.nan.

    """
    try:
        m = Chem.MolFromSmiles(s)
    except:
        m = np.nan
    return m


def genc_smi(m):
    """
    Generates a RDkit canonical SMILES string from a molecule object.

    Args:
        m (rdkit.Chem.Mol): A molecule object, from the RDKit library.

    Returns:
        str or np.nan: The SMILES string representation of the molecule if successful, 
        otherwise returns np.nan.

    """
    try:
        cS = Chem.MolToSmiles(m)
    except:
        cS = np.nan
    return cS

# count the number of the targetted functional group


def count_fg(m, patt):
    """
    Counts the number of functional groups (FG) in a molecule 
    based on a given pattern. 

    Args:
        m (rdkit.Chem.Mol): The molecule object to search
            for substructure matches.
        patt (rdkit.Chem.Mol): The pattern molecule used
            to identify substructure matches.

    Returns:
        int: The number of functional groups identified in the molecule.

    """
    numFG = 0
    matchs = m.GetSubstructMatches(patt)
    if len(matchs) >= 2:
        not_match = []
        for i in range(0, len(matchs)-1):
            if len(set(matchs[i]) ^ set(matchs[i+1])) != 2:
                pass
            else:
                not_match.append(i+1)
            numFG = len([matchs[i]
                        for i in range(0, len(matchs)) if i not in not_match])
    else:
        numFG = len(matchs)
    return numFG

# classify candidate compounds for mono-FG monomer


def monomer_sel_mfg(m, mons, excls):
    """
    Determining whether the given small molecule compound 
    qualifies as a self-polymerizable monomer and 
    categolize it into a monomer class. 

    Args:
        m (rdkit.Chem.Mol): The molecule to be analyzed. 
            If None or NaN, the function returns default values.
        mons (list of str): A list of SMARTS strings representing 
            monomer patterns to match against the molecule.
        excls (list of str): A list of SMARTS strings representing 
            exclusion patterns to check against the molecule.

    Returns:
        list: A list containing:  

            - fchk (bool): True if the molecule matches any monomer pattern 
            and does not match any exclusion pattern, otherwise False.  

            - fchk_c (int): The total count of substructure matches 
            for all monomer patterns.

    """
    if pd.notna(m):
        chk_c = 0
        fchk_c = 0
        chk = []
        if len(mons) != 0:
            for mon in mons:
                patt = Chem.MolFromSmarts(mon)
                if m.HasSubstructMatch(patt):
                    chk_c = len(m.GetSubstructMatches(patt))
                    fchk_c = fchk_c+chk_c
                    chk_excl = []
                    for excl in excls:
                        excl_patt = Chem.MolFromSmarts(excl)
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
    return [fchk, fchk_c]


# classify candidate compounds for poly-FG monomer
# count objective FGs
def monomer_sel_pfg(m, mons, excls, minFG, maxFG):
    """
    Determining whether the given small molecule compound qualifies as 
    a monomer or not. If so, count a number of polymerizeble 
    functional group and categolize it into a monomer class. 

    Args:
        m (rdkit.Chem.Mol): The monomer molecule to evaluate.
        mons (list of str): A list of SMARTS patterns representing the 
            functional groups to count in the monomer.
        excls (list of str): A list of SMARTS patterns representing the 
            exclusion patterns to check against the monomer.
        minFG (int): The minimum number of functional groups required.
        maxFG (int): The maximum number of functional groups allowed.

    Returns:
        list: A list containing:  
        
            - fchk (bool): True if the monomer satisfies the conditions,
            False otherwise.  
            
            - fchk_c (int): The total count of functional groups found in
            the monomer.

    """
    if pd.notna(m):
        chk_c = 0
        fchk_c = 0
        if len(mons) != 0:
            for mon in mons:
                patt = Chem.MolFromSmarts(mon)
                chk_c = count_fg(m, patt)
                fchk_c = fchk_c + chk_c
            if minFG <= fchk_c <= maxFG:
                chk = []
                for excl in excls:
                    excl_patt = Chem.MolFromSmarts(excl)
                    if m.HasSubstructMatch(excl_patt):
                        chk.append(False)
                    else:
                        chk.append(True)
                if False in chk:
                    fchk = False
                else:
                    fchk = True
            else:
                fchk = False
        else:
            fchk = (False)
    else:
        fchk = (False)
    return [fchk, fchk_c]


# define sequential polymerization for chain polymerization except polyolefine
def seq_chain(prod_P, targ_mon1, Ps_rxnL, mon_dic, monL):
    """
    This function applied to multifunctional monomers 
    for chain polymerization except polyolefine.
    Processes a molecular structure by applying a sequential reactions 
    based on specific substructure matches.

    Args:
        prod_P (rdkit.Chem.Mol): The input molecule to be processed.
        targ_mon1 (str): Target monomer type, used to
            determine processing logic.
        Ps_rxnL (dict): A dictionary of polymerization reaction objects
            indexed by integers.
        mon_dic (dict): A dictionary containing monomer class
            (not used in this function).
        monL (list): A list of monomer SMARTS patterns indexed by integers.

    Returns:
        rdkit.Chem.Mol: The processed molecule after applying the reactions.

    """
    if Chem.MolToSmiles(prod_P) != '':
        if targ_mon1 not in ['vinyl', 'cOle']:
            seqFG2 = Chem.MolFromSmarts(monL[[202][0]])
            seqFG3 = Chem.MolFromSmarts(monL[[203][0]])
            seqFG4 = Chem.MolFromSmarts(monL[[204][0]])
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
            prod_P = prod_P
    return prod_P


# define sequential polymerization for successive polymerization
def seq_successive(prod_P, targ_rxn, monL, Ps_rxnL, P_class):
    """
    This function applied to multifunctional monomers for 
    successive polymerization. 
    Processes a molecular structure by applying a sequential reactions 
    based on specific substructure matches. 

    Args:
        prod_P (rdkit.Chem.Mol): The product molecule to be processed.
        targ_rxn (Any): Target reaction
            (not used in the current implementation).
        monL (list): A list containing SMARTS patterns for functional groups.
        Ps_rxnL (list): A list of polymerization reaction objects
            to be applied to the product molecule.
        P_class (str): The polymer class of the product molecule,
            which determines the reaction sequence.

    Returns:
        rdkit.Chem.Mol: The processed product molecule 
        after applying the reaction sequence.

    Notes:
        - The function uses substructure matching to determine 
          which reactions to apply.
        - The behavior of the function depends on the `P_class` 
          of the molecule.
        - Specific reaction sequences are applied for classes 
          such as 'polyolefin', 'polyoxazolidone', 'polyimide',
          and 'polyester'.
        - If the `P_class` is not recognized,
          the product molecule is returned unchanged.

    """
    if Chem.MolToSmiles(prod_P) != '':
        seqFG0 = Chem.MolFromSmarts(monL[[200][0]])
        seqFG1 = Chem.MolFromSmarts(monL[[201][0]])
        seqFG2 = Chem.MolFromSmarts(monL[[202][0]])
        seqFG3 = Chem.MolFromSmarts(monL[[203][0]])
        seqFG4 = Chem.MolFromSmarts(monL[[204][0]])
        seqFG5 = Chem.MolFromSmarts(monL[[205][0]])
        seqFG6 = Chem.MolFromSmarts(monL[[206][0]])
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
                if P_class == 'polyimide':
                    prods = Ps_rxnL[207].RunReactants([prod_P])
                elif P_class == 'polyester':
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
            prod_P = prod_P
    return prod_P


# homopolymerization
def homopolymA(mon1, mons, excls, targ_mon1, Ps_rxnL, mon_dic, monL):
    """
    Generates a polymer CRU formed from a single monomer by
    iteratively reacting a monomer until no further reactions are possible.

    Args:
        mon1 (rdkit.Chem.Mol): The initial monomer
            to start the polymerization process.
        mons (list): A list of SMARTS strings representing monomer patterns to
            match against the molecule.
        excls (list): A list of SMARTS strings representing exclusion patterns
            to check against the molecule.
        targ_mon1 (object): The target monomer class for
            the polymerization process.
        Ps_rxnL (list): A dictionary of polymerization reaction objects
            indexed by integers.
        mon_dic (dict):  A dictionary containing monomer class. 
        monL (list): A list of monomer SMARTS patterns indexed by integers.

    Returns:
        list: A list of SMILES strings representing the generated homopolymers.

    """
    prod_P = mon1
    # 生成したポリマーがさらに重合可能な場合、再度反応
    while monomer_sel_mfg(prod_P, mons, excls)[0] == True:
        prods = Ps_rxnL[mon_dic[targ_mon1]].RunReactants([prod_P])
        prod_Ps = []
        for prod_P in prods:
            try:
                Chem.SanitizeMol(prod_P[0])
                prod_P = prod_P[0]
                prod_P = seq_chain(prod_P, targ_mon1=targ_mon1,
                                   Ps_rxnL=Ps_rxnL, mon_dic=mon_dic, monL=monL)
                prod_Ps.append(prod_P)
            except:
                pass
    return [genc_smi(m) for m in prod_Ps]


# binarypolymerization
def bipolymA(reactant, targ_rxn, monL, Ps_rxnL, P_class):
    """
    Generates a polymer CRU formed from two monomers by iteratively reacting 
    a monomer until no further reactions are possible. 

    Args:
        reactant (tuple): A tuple of reactant molecules
            to be used in the reaction.
        targ_rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            The target chemical reaction to apply.
        monL (list):  A list of monomer SMARTS patterns indexed by integers.
        Ps_rxnL (dict): A list of monomer SMARTS patterns indexed by integers.
        P_class (type): A class type used for polymer processing.

    Returns:
        list: A list of SMILES strings representing 
        the generated polymer products.

    """
    prod_P = Chem.MolFromSmiles('')
    prods = targ_rxn.RunReactants(reactant)
    prod_Ps = []
    for prod_P in prods:
        try:
            Chem.SanitizeMol(prod_P[0])
            prod_P = prod_P[0]
            prod_P = seq_successive(
                prod_P,
                targ_rxn=targ_rxn,
                monL=monL,
                Ps_rxnL=Ps_rxnL,
                P_class=P_class
            )
            prod_Ps.append(prod_P)
        except:
            pass
    return [genc_smi(m) for m in prod_Ps]


# Copyright (c) 2024 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# (license that can be found in the LICENSE file. )

# 08/15/2024, M. Ohno
# generate CRU of olefinic polymers
# NEED def ole_rxnsmarts_gen(reactant)

def ole_rxnsmarts_gen(reactant):
    """
    Generates a polymerization reaction SMARTS string 
    for a given olefinic monomer. Place this function
    right before def ole_cru_gen() so that it can be used 
    within the function ole_cru_gen.

    Args:
        reactant (str): The input reactant string in SMARTS format.

    Returns:
        str: The reaction SMARTS string representing the transformation
        from the reactant to the product.

    """
    prod = ''
    prod1 = ''
    prod2 = ''
    prod3 = ''
    prod4 = ''
    inv_reactant = reactant[::-1]
    C1_i = inv_reactant.find(':1]'[::-1])
    C1_j = inv_reactant.find('[CX3'[::-1], C1_i)
    C2_i = inv_reactant.find(':2]'[::-1])
    C2_j = inv_reactant.find('=[CX3'[::-1], C2_i)
    prod1 = inv_reactant[:C2_i]+'(-*)'[::-1]
    prod2 = inv_reactant[C2_i:C2_j]+'-[CX4'[::-1]
    prod3 = inv_reactant[C2_j+5:C1_i]+'(-*)'[::-1]
    prod4 = inv_reactant[C1_i:C1_j]+'[CX4'[::-1]+inv_reactant[C1_j+4:]
    prod = prod4[::-1] + prod3[::-1] + prod2[::-1] + prod1[::-1]
    rxn_smarts = reactant + '>>' + prod
    return rxn_smarts


def ole_cru_gen(m, mon):
    """
    Generates a CRU from olefinic monomer by applying a reaction
    iteratively until no further reactions are possible.

    Args:
        m (rdkit.Chem.Mol): The input molecule to which
            the reaction will be applied.
        mon (str): A SMARTS string representing the monomer pattern.
    Returns:
        list: A list containing:
            - rdkit.Chem.Mol: The final CRU after all reactions.
            - list of str: A list of SMILES strings for CRUs.
    Raises:
        Exception: If there is an issue with sanitizing the molecule 
        during reaction processing.

    """
    reactant = [m, ]
    patt = Chem.MolFromSmarts(mon)
    targ_rxn = AllChem.ReactionFromSmarts(ole_rxnsmarts_gen(mon))
    targ_rxn.Initialize()  # need initialization
    while m.HasSubstructMatch(patt):
        prod_P = Chem.MolFromSmiles('')
        prods = targ_rxn.RunReactants(reactant)
        prod_Ps = []
        for prod_P in prods:
            try:
                Chem.SanitizeMol(prod_P[0])
                prod_P = prod_P[0]
                prod_Ps.append(prod_P)
            except:
                pass
        m = prod_Ps[0]
        reactant = [m, ]
    return [m, [genc_smi(m) for m in prod_Ps]]


def diene_12to14(smi, rxn):
    """
    Convert the structure of the 1,2-adducted CRU to a 1,4-adduct.
    Place this function right before def diene_14() so that
    it can be used within the function diene_14.

    Args:
        smi (str): The input SMILES string containing
            asterisks (*) as placeholders.
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            Ps_rxnL[209] was applied.

    Returns:
        str: The resulting SMILES string after the reaction, 
        with placeholders 
        replaced back to asterisks (*).

    Raises:
        rdkit.Chem.rdchem.KekulizeException: 
        If the molecule sanitization fails.
        IndexError: If the reaction does not produce any products.

    """
    # rxn == Ps_rxnL[209]
    # Replace all asterisks with [3H] in the SMILES string:
    repl_smi = smi.replace("*", "[3H]")
    new_m = Chem.MolFromSmiles(repl_smi)
    Chem.SanitizeMol(new_m)
    prods = rxn.RunReactants([new_m])
    for m in prods:
        Chem.SanitizeMol(m[0])
    new_smi = Chem.MolToSmiles(prods[0][0]).replace("[3H]", "*")
    return new_smi


def diene_14(x, rxn):
    """
    Generate 1,4-addition CRU from  a conjugated diene monomer. 

    Args:
        x (dict): The results of olefin classification and
            the chemical structure of these CRU generated by ole_sel_cru.
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            Ps_rxnL[209] was applied.

    Returns:
        dict: The modified dictionary `x` with the transformed SMILES string in 
        `x['conjdiene'][2]`, if applicable. If `'conjdiene'` is not present or 
        empty, the dictionary is returned unchanged.

    """
    if 'conjdiene' in x and x['conjdiene'][0]:
        x['conjdiene'][2] = diene_12to14(x['conjdiene'][2], rxn)
    return x

 

# Copyright (c) 2024 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# (license that can be found in the LICENSE file. )

# 08/15/2024, M. Ohno
# classify olefinic monomers and generate CRU


def ole_sel_cru(m, mons, excls, minFG, maxFG):
    """
    Selects and processes a molecule based on specific criteria and 
    generates a SMILES representation.

    Args:
        m (rdkit.Chem.Mol): The molecule to be processed.
        mons (list of str): A list of SMARTS patterns representing
            the functional groups to count in the monomer.
        excls (list of str): A list of SMARTS patterns representing
            the exclusion patterns to check against the monomer.
        minFG (int): The minimum number of olefinic polymerizable site
            required.
        maxFG (int): The maximum number of olefinic polymerizable site allowed.

    Returns:
        list: A list containing:
            - The result of the `monomer_sel_pfg` function
              (list of bool and other values).
            - The SMILES representation of the processed molecule (str).

    """
    judge = monomer_sel_pfg(m, mons, excls, minFG, maxFG)
    if judge[0] == True:
        for mon in mons:
            patt = Chem.MolFromSmarts(mon)
            if m.HasSubstructMatch(patt):
                CRU = ole_cru_gen(m, mon)
                m = CRU[0]
            else:
                pass
    else:
        m = np.nan
    smi_p = genc_smi(m)
    judge.append(smi_p)
    return judge


def update_nested_dict(row, dict_col, new_val, updated_k):
    """
    Used in the function 'olecls'. If the classification result for 
    the olefin class is True, writes the number of functional groups and 
    the SMILES notation of the CRU.

    Args:
        row (dict): The dictionary representing a row of data.
        dict_col (str): The key in the row that contains
            the nested dictionary to be updated.
        new_val (str): The key in the row whose value will be assigned to
            the nested dictionary.
        updated_k (str): The key in the nested dictionary to be updated.

    Returns:
        dict: The updated row with the modified nested dictionary.

    """
    row[dict_col][updated_k] = row[new_val]
    return row


def coord_polym(smi, targ_rxn):
    """
    Generate a list of CRUs for olefin copolymer polymers from  
    the input SMILES string with a target reaction.

    Args:
        smi (str): The SMILES string of the input molecule.
        targ_rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            The target reaction
            to apply to the input molecule.

    Returns:
        list: A list of unique SMILES strings representing 
        the products of the reaction.

    """
    prods = targ_rxn.RunReactants([genmol(smi),])
    prod_Ps = []
    for prod_P in prods:
        try:
            Chem.SanitizeMol(prod_P[0])
            prod_P = prod_P[0]
            prod_Ps.append(prod_P)
        except:
            pass
    rsl = set([genc_smi(m) for m in prod_Ps])
    return list(rsl)

# end
