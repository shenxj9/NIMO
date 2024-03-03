import re
from collections import Counter
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import HeavyAtomCount
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem.BRICS import FindBRICSBonds,BreakBRICSBonds
import ast
import re

def f_index(lst,id):
    return [i for i, x in enumerate(lst) if x[0] == id]

def swap(s1, s2, l):
    seg1=l[:s1.start]+l[s2]
    seg2=l[s1.stop : s2.start]
    seg3=l[s1]+l[s2.stop:]
    return seg1+seg2+seg3

def Recursive_Fragment_recombination1(n, frags, map_id_list,Standardized):
    mols1 = [Chem.MolFromSmiles(smi) for smi in frags]
    rxn = AllChem.ReactionFromSmarts("[$(*[#2:1]):1][#2].[$(*[#2:1]):2][#2]>>[*:1][*:2]")
    if n == 0:
        return frags[0], Standardized
    else:
        regex = re.compile(r"(\d+)\*")
        map_id_n = regex.findall(frags[n])
        map_id_n_1 = regex.findall(frags[n - 1])
        if len(map_id_n) == 0 or len(map_id_n_1) == 0:
            return 'CC', 'CC'
        if map_id_n[0] != map_id_n_1[0]:
            swap_id = f_index(map_id_list, map_id_n[0])[0]
            if swap_id != n:
                frags[n - 1], frags[swap_id] = frags[swap_id], frags[n - 1]
                map_id_list[n - 1], map_id_list[swap_id] = map_id_list[swap_id], map_id_list[n - 1]
            else:
                return None, None
        Standardized.append(frags[n - 1])
        frags[n] = re.sub(r'\[\d+\*\]', "[#2]", frags[n], count=1)
        frags[n - 1] = re.sub(r'\[\d+\*\]', "[#2]", frags[n - 1], count=1)
        mols1[n] = Chem.MolFromSmiles(frags[n])
        mols1[n - 1] = Chem.MolFromSmiles(frags[n - 1])
        mols1[n - 1] = rxn.RunReactants((mols1[n], mols1[n - 1]))
        if len(mols1[n - 1]) is 0:
            mols1[n - 1] = Chem.MolFromSmiles(frags[n - 1])
            rxn = AllChem.ReactionFromSmarts("([#2]=[C:1]).([#2]=[C:2])>> [C:1]=[C:2]")
            mols1[n - 1] = rxn.RunReactants((mols1[n], mols1[n - 1]))
        frags[n - 1] = Chem.MolToSmiles(mols1[n - 1][0][0])
        return Recursive_Fragment_recombination1((n - 1), frags, map_id_list,Standardized)

def Recursive_Fragment_recombination2(n, frags, Standardized,R_list,len_dict):
    mols1 = [Chem.MolFromSmiles(smi) for smi in frags]
    rxn = AllChem.ReactionFromSmarts("[$(*[#2:1]):1][#2].[$(*[#2:1]):2][#2]>>[*:1][*:2]")
    if n == 0:
        return frags[0], Standardized
    else:
        regex = re.compile(r"\*\:(\d+)")
        true_R_list = regex.findall(frags[n])
        front = true_R_list[0]
        true_next = regex.findall(frags[n-1])
        if len(true_R_list) >= 2 and (front == true_next[0] == true_R_list[1]):
                start = len_dict[front].start
                stop = len_dict[front].stop - 1
                len_dict[front] = slice(start, stop, None)
        elif front == true_next[0]:
             pass
        else:
            back = true_next[0]
            frags = swap(len_dict[front], len_dict[back], frags)
            front_index = R_list.index(front)
            back_index = R_list.index(back)
            if true_R_list.count(front)==1 & true_R_list.count(back)==1:
                R_list[front_index], R_list[back_index] = R_list[back_index], R_list[front_index]
                len_dict[back], len_dict[front] = len_dict[front], len_dict[back]
            else:
                R_list_R = R_list[::-1]
                R_list_R = swap(len_dict[front], len_dict[back], R_list_R)
                R_list = R_list_R[::-1]
                regex = re.compile(r"\*\:(\d+)")
                true_R_list = regex.findall(''.join(frags[:n - 1]))
                R_list1 = list(set(true_R_list))
                R_list1.sort(key=true_R_list.index)
                R_count = dict(Counter(true_R_list))
                len_dict = {}
                start = 0
                for r in R_list1:
                    end = start + R_count[r]
                    len_dict[r] = slice(start, end)
                    start = end
        Standardized.append(frags[n - 1])
        frags[n] = re.sub(r'\[\*\:\d+\]', "[#2]", frags[n], count=1)
        frags[n - 1] = re.sub(r'\[\*\:\d+\]', "[#2]", frags[n - 1], count=1)

    mols1[n] = Chem.MolFromSmiles(frags[n])
    mols1[n - 1] = Chem.MolFromSmiles(frags[n - 1])
    if mols1[n] == None or mols1[n - 1] == None:
        return 'CC', 'CC'
    mols1[n - 1] = rxn.RunReactants((mols1[n], mols1[n - 1]))
    if len(mols1[n - 1]) is 0:
        return 'CC', 'CC'
    frags[n - 1] = Chem.MolToSmiles(mols1[n - 1][0][0])
    return Recursive_Fragment_recombination2((n - 1), frags, Standardized,R_list,len_dict)

def Fragment_recombination1(true_smiles,frags,n):

    i = 1
    pre_smiles = None
    while pre_smiles == None:
        initial_frags = frags[:]
        Standardized = []
        Standardized.append(frags[0])
        regex = re.compile(r"(\d+)\*")
        map_id_list = [regex.findall(frag) for frag in frags]
        map_id_list.reverse()
        frags.reverse()
        pre_smiles,Standardized_frags = Recursive_Fragment_recombination1((n-1),frags,map_id_list,Standardized)
        initial_frags = initial_frags[1:] + [initial_frags[0]]
        frags = initial_frags[:]
        i = i+1
        if i == n+1:
            print("Error arised , special standardization required")
            return None,False,False
    pre_smiles = re.sub(r'\[CH?]', 'C', pre_smiles)
    true_smiles = re.sub(r'\[CH?]', 'C', true_smiles)
    true_smiles_unchiral = Chem.MolToSmiles(Chem.MolFromSmiles(true_smiles), isomericSmiles=False)
    pre_smiles_unchiral = Chem.MolToSmiles(Chem.MolFromSmiles(pre_smiles), isomericSmiles=False)
    return Standardized_frags,pre_smiles == true_smiles, true_smiles_unchiral == pre_smiles_unchiral

def Fragment_recombination2(true_smiles,frags,n,R_list):

    Standardized= []
    Standardized.append(frags[0])
    frags.reverse()
    R_list_R = R_list[::-1]
    R_list1 = list(set(R_list_R))
    R_list1.sort(key=R_list_R.index)
    R_count = dict(Counter(R_list_R))
    len_dict = {}
    start = 0
    for r in R_list1:
        end = start + R_count[r]
        len_dict[r] = slice(start, end)
        start = end
    pre_smiles,Standardized = Recursive_Fragment_recombination2((n-1),frags,Standardized,R_list,len_dict)
    true_smiles_unchiral = Chem.MolToSmiles(Chem.MolFromSmiles(true_smiles),isomericSmiles=False)
    pre_smiles_unchiral = Chem.MolToSmiles(Chem.MolFromSmiles(pre_smiles), isomericSmiles=False)
    return Standardized,pre_smiles == true_smiles, true_smiles_unchiral == pre_smiles_unchiral

def Motif_Fragment(mol,smi):
    # print("*******************")
    # print('smi:', smi)
    mol = Chem.MolFromSmiles(smi)
    AtomsPairs1 = list(mol.GetSubstructMatches(Chem.MolFromSmarts('[R][R]')))
    AtomsPairs2 = list(mol.GetSubstructMatches(Chem.MolFromSmarts('[R:1]@[R:2]')))
    AtomsPairs3 = [x for x in AtomsPairs1 if x not in AtomsPairs2]
    AtomsPairs4 = list(mol.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]')))
    BRICSBonds = list(FindBRICSBonds(mol))
    AtomsPairs5 = [BRICSBonds[i][0] for i in range(len(BRICSBonds))]
    AtomsPairs = AtomsPairs3 + AtomsPairs4 + AtomsPairs5
    AtomsPairs  = [tuple(sorted(x)) for x in AtomsPairs ]
    AtomsPairs = list(set(AtomsPairs))
    Bonds_id = [mol.GetBondBetweenAtoms(i[0], i[1]).GetIdx() for i in AtomsPairs]
    str_Bonds_id = [(str(Bonds_id[i]), str(Bonds_id[i])) for i in range(len(Bonds_id))]
    new_BRICSBonds = [(AtomsPairs[i], str_Bonds_id[i]) for i in range(len(AtomsPairs))]
    mols = BreakBRICSBonds(mol, new_BRICSBonds)
    frags = Chem.MolToSmiles(mols).split('.', -1)
    frags = [re.sub(r'(?<!\d)\*', "[0*]", frag) for frag in frags]
    n = len(frags)
    if n ==1:
        return smi,smi,n,True,True
    Standardized_frags, chiral_equal, unchiral_equal = Fragment_recombination1(smi, frags, n)
    return smi,Standardized_frags,n,chiral_equal, unchiral_equal

def Scaffold_Fragment(mol, smi):
#     print("*******************")
#     print('smi:',smi)
    core = MurckoScaffold.GetScaffoldForMol(mol)
    res, unmatched = rdRGD.RGroupDecompose([core], [mol], asSmiles=True)
    if len(res) == 0:
        str = [smi]
    else:
        str = res[0].values()
    frags = [s.split('.', -1) for s in str]
    frags = [b for a in frags for b in a]
    for smile in frags:
        if HeavyAtomCount(Chem.MolFromSmarts(smile)) == 0:
            frags.remove(smile)
    regex = re.compile(r"\:(\d+)")
    old_list = regex.findall(''.join(frags[1:]))
    old_list = list(map(int, old_list))
    old_list_counter= dict(Counter(old_list))
    dup_old_list = [key for key, value in old_list_counter.items() if value > 1]
    rule = regex.findall(''.join(frags[0]))
    rule = list(map(int, rule))
    key = [i for i in range(len(rule))]
    rule = dict(zip(rule, key))
    new_list = [rule[i] for i in old_list]
    new_list_count = dict(Counter(new_list))
    dup_num = [key for key, value in new_list_count.items() if value > 1]
    dup_num.sort(reverse=True)
    dup_key = [new_list.index(i) for i in dup_num]
    for j in dup_num:
        new_list = [i + 1 if i >= j else i for i in new_list]
    for i in dup_key:
        new_list[i] = new_list[i] - 1
    d = dict(zip(new_list, frags[1:]))
    Standardized_frags = list(dict(sorted(d.items(), key=lambda d: d[0], reverse=False)).values())
    Standardized_frags.insert(0, frags[0])
    if dup_old_list is not None:
        sub0 = [f'\(\[\*\:{i}\]\)' for i in dup_old_list]
        sub1 = [f'\[\*\:{i}\]' for i in dup_old_list]
        sub2 = [f'([*:{i}])' * 2 for i in dup_old_list]
        for i in range(len(dup_old_list)):
            if re.search(sub0[i], Standardized_frags[0]) is not None:
                Standardized_frags[0] = re.sub(sub0[i], sub2[i], Standardized_frags[0])
            else:
                Standardized_frags[0] = re.sub(sub1[i], sub2[i], Standardized_frags[0])
    regex = re.compile(r"\*\:(\d+)")
    R_list =regex.findall(Standardized_frags[0])
    R_list.insert(0,'0')
    n = len(Standardized_frags)
    Standardized_frags, chiral_equal, unchiral_equal = Fragment_recombination2(smi, Standardized_frags, n, R_list)
    return smi,Standardized_frags,n,chiral_equal, unchiral_equal

def del_mapid(x):
    x = re.sub(r'\[\*\:\d+]', "[*]", f'{x}')
    x = re.sub(r'\[\d+\*\]', "[*]", f'{x}')
    x = ast.literal_eval(x)
    x = ' '.join(x)
    return x
