import argparse
import pandas as pd
import json
import csv
import os
from onmt.utils.properties import add_property,mols_from_smiles,condition_convert
from onmt.utils.motif_extration import Motif_Fragment, Scaffold_Fragment, del_mapid
from sklearn.model_selection import train_test_split

from itertools import chain

dummy_list = dict()
for i in range(30):
    dummy_list[i + 1] = (f'dummy_{i + 1}')

def add_dummy(line):

    line = line.replace('\t', ' ')
    line = line.split(' ')
    dummy = '[*]'
    number = [i.count(dummy) for i in line[3:]]
    dummy_number = [dummy_list[i] for i in number]
    dummy_line = line[:3]+list(chain.from_iterable(zip(dummy_number, line[3:])))
    return ' '.join(dummy_line)

def build_specific_fragment(specific_fragment):
    if specific_fragment  == "Nimom":
        return Motif_Fragment
    elif specific_fragment  == "Nimos":
        return Scaffold_Fragment
    else:
        raise ValueError(f"No Fragmentation mode defined for {specific_fragment}")

def add_fragments(dataset, specific_fragment):
    smiles = dataset.smiles.tolist()
    mols = mols_from_smiles(smiles)
    results = []
    fun = build_specific_fragment(specific_fragment)
    for m, s in zip(mols, smiles):
        results.append(fun(m, s))
    smiles,fragments, lengths,chiral_equal, unchiral_equal = zip(*results)
    dataset["smiles"] = smiles
    dataset["fragments"] = fragments
    dataset["n_fragments"] = lengths
    dataset["chiral_equal"] = chiral_equal
    dataset["unchiral_equal"] = unchiral_equal
    return dataset


def save_dataset(dataset, save_path, info):
    testset = dataset[dataset.fragments.notnull()]
    trainset = testset[testset.n_fragments >= info['min_length']]
    trainset = trainset[trainset.n_fragments <= info['max_length']]
    trainset = trainset[trainset['unchiral_equal'] == True]

    trainset['fragments'] = trainset['fragments'].apply(lambda x: del_mapid(x))
    trainset.to_csv(save_path+'\\dataset_.smi', index=False,sep = '\t',quoting=csv.QUOTE_NONE,escapechar='\t',columns = ['logP','Qed','SAS','fragments'],header = None)

    train_data, test_data = train_test_split(trainset, test_size=0.1, random_state=1000)
    train_data.to_csv(save_path+'\\train_.smi', index=False,sep = '\t',quoting=csv.QUOTE_NONE,escapechar='\t',columns = ['logP','Qed','SAS','fragments'],header = None)
    test_data.to_csv(save_path +'\\valid_.smi', index=False,sep = '\t',quoting=csv.QUOTE_NONE,escapechar='\t',columns = ['logP','Qed','SAS','fragments'],header = None)

    filename1 =[save_path+'\\dataset_.smi',save_path+'\\train_.smi',save_path +'\\valid_.smi' ]
    filename2 = [save_path + '\\dataset.smi', save_path + '\\train.smi', save_path + '\\valid.smi']
    for i, j in zip(filename1, filename2):
        f1 = open(i, 'r')
        f2 = open(j, 'w+')
        for line in f1:
            f2.write(add_dummy(line))
        f1.close()
        f2.close()
        os.remove(i)

def preprocess_dataset(raw_data,save_path,specific_fragment,n_jobs):
    json_path = 'onmt/utils/properties.json'
    info =  json.load(open(json_path, "r"))
    dataset = pd.read_csv(raw_data, encoding="gbk", sep='\t', header=0, usecols=["smiles"],nrows=300) #,nrows=100
    for prop in info['properties']:
        if prop not in dataset.columns:
            dataset = add_property(dataset, prop,n_jobs)
    dataset = condition_convert(dataset)
    dataset = add_fragments(dataset, specific_fragment)
    save_dataset(dataset, save_path,info)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Main script for running the model")
    parser.add_argument('--raw-data', action='store', dest='raw_data')
    parser.add_argument('--save-path', action='store', dest='save_path')
    parser.add_argument('--specific_fragment', action='store', dest='specific_fragment',default='Nimom', choices=['Nimom','Nimos'])
    parser.add_argument('--n_jobs', default=1)
    arg_dict = vars(parser.parse_args())
    print(arg_dict)
    preprocess_dataset(**arg_dict)




