
import pandas as pd
from joblib import Parallel, delayed

from rdkit.Chem import Crippen, QED
from rdkit.Chem import  Descriptors, rdMolDescriptors
from rdkit import Chem
from onmt.sascorer.sascorer import calculateScore
import pandas as pd
def condition_convert(con_df):

    con_df['logP'][con_df['logP'] < 5] = 1
    con_df['logP'][con_df['logP'] >= 5] = 0
    con_df['Qed'][con_df['Qed'] >= 0.5] = 1
    con_df['Qed'][con_df['Qed'] < 0.5] = 0
    con_df['SAS'][con_df['SAS'] <= 5] = 1
    con_df['SAS'][con_df['SAS'] > 5] = 0

    con_df['logP'][con_df['logP'] == 1] = 'good_logP'
    con_df['logP'][con_df['logP'] == 0] = 'bad_logP'
    con_df['Qed'][con_df['Qed'] == 1] = 'high_QED'
    con_df['Qed'][con_df['Qed'] == 0] = 'low_QED'
    con_df['SAS'][con_df['SAS'] == 1] = 'good_SA'
    con_df['SAS'][con_df['SAS'] == 0] = 'bad_SA'

    return con_df

def canonicalize(smi, clear_stereo=False):
    mol = Chem.MolFromSmiles(smi)
    if clear_stereo:
        Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def mol_from_smiles(smi):
    smi = canonicalize(smi)
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol)
    return mol

def mols_from_smiles(mols):
    return [mol_from_smiles(m) for m in mols]


def logp(mol):
    return Crippen.MolLogP(mol) if mol else None


def qed(mol):
    return QED.qed(mol) if mol else None


def sas(mol):
    return calculateScore(mol) if mol else None


def mr(mol):
    return Crippen.MolMR(mol) if mol else None


def Mw(mol):
    return Descriptors.MolWt(mol) if mol else None

def HBD(mol):
    return rdMolDescriptors.CalcNumLipinskiHBD(mol) if mol else None

def HBA(mol):
    return rdMolDescriptors.CalcNumLipinskiHBA(mol) if mol else None

def RB(mol):
    return rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else None

def AromaticRings(mol):
    return rdMolDescriptors.CalcNumAromaticRings(mol) if mol else None

def FractionCSP3(mol):
    return rdMolDescriptors.CalcFractionCSP3(mol) if mol else None

def TPSA(mol):
    return Descriptors.TPSA(mol) if mol else None



def add_property(dataset, name, n_jobs):

    # fndict = {"qed": qed, "SAS": sas, "logP": logp, "Mr": mr,"Mw":Mw,"HBD":HBD,"HBA":HBA,"RB":RB,"TPSA":TPSA,"AromaticRings":AromaticRings,"FractionCSP3":FractionCSP3}
    fndict = {"Qed": qed, "SAS": sas, "logP": logp}
    fn = fndict[name]
    smiles = dataset.smiles.tolist()
    mols = mols_from_smiles(smiles)
    pjob = Parallel(n_jobs=n_jobs, verbose=0)
    prop = pjob(delayed(fn)(mol) for mol in mols)
    new_data = pd.DataFrame(prop, columns=[name])
    return pd.concat([dataset, new_data], axis=1, sort=False)



