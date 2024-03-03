import re
import argparse
from rdkit.Chem import AllChem
from rdkit import Chem

def Fragment_recombination(n, frags):
    mols1 = [Chem.MolFromSmiles(smi) for smi in frags]
    rxn = AllChem.ReactionFromSmarts("[$(*[#2:1]):1][#2].[$(*[#2:1]):2][#2]>>[*:1][*:2]")
    if n == 0:
        return frags[0]
    else:
        frags[n] = re.sub(r'\*', "He", frags[n], count=1)
        frags[n - 1] = re.sub(r'\*', "He", frags[n - 1], count=1)
        frags[n] = re.sub(r'(?<!\[)He(?<!\])', "[He]", frags[n])
        frags[n - 1] = re.sub(r'(?<!\[)He(?<!\])', "[He]", frags[n - 1])
        mols1[n] = Chem.MolFromSmiles(frags[n])
        mols1[n - 1] = Chem.MolFromSmiles(frags[n - 1])
        if mols1[n] == None or mols1[n - 1] == None:
            print('erro')
            return 'CC', 'CC'
        mols1[n - 1] = rxn.RunReactants((mols1[n], mols1[n - 1]))
        if len(mols1[n - 1]) is 0:
            mols1[n - 1] = Chem.MolFromSmiles(frags[n - 1])
            rxn = AllChem.ReactionFromSmarts("([#2]=[C:1]).([#2]=[C:2])>> [C:1]=[C:2]")
            mols1[n - 1] = rxn.RunReactants((mols1[n], mols1[n - 1]))
        if len(mols1[n - 1]) == 0:
            return None
        frags[n - 1] = Chem.MolToSmiles(mols1[n - 1][0][0])
        return Fragment_recombination((n - 1), frags)

def remove_ele(src,ele):
    n=len(src)
    bek=[]
    for i in range(n):
        if ele in src[i]:
            bek.append(i)

    src = [src[i] for i in range(n) if (i not in bek)]
    return src

def recombination(input,output,scaffold=None ):
    F = T = num = 0
    moleculars = []
    id = []
    s = 0
    for line in input:
        s += 1
        line = line.strip('\n')
        line = line.split(' ')
        # line = line.split(' ')[1:]
        Standardized_frags = remove_ele(line,'dummy')
        if scaffold != None:
            Standardized_frags.insert(0, scaffold)
        Standardized_frags.reverse()
        generation = Fragment_recombination(len(Standardized_frags)-1,Standardized_frags)
        num = num+1
        if generation is not None and '#2' not in generation and 'He' not in generation and '*' not in generation :
            T += 1
            moleculars.append(generation)
            id.append(s)
        else:
            F += 1
    print('Valid_percent: {:.2%}'.format(T/(T+F)))
    output.write('id'+ ',' + 'smiles' + '\n')
    for (molecular, s) in zip(moleculars, id):
        if isinstance(molecular,str):
            output.write(str(s)+','+molecular+'\n')
            output.flush()












