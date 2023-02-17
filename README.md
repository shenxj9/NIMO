# NIMO
![img_1.png](img_1.png)

This is the code for the "Nimo: a Nature Inspired Molecular Generative Model Based on Fragments of Natural Products" paper.

The package is based on [OpenNMT-py 2.0](http://opennmt.net/OpenNMT-py/).

## Environment
- python = 3.8

## Installation Instructions

1. Clone this repository
2. Install the following python packages:
* pip install pandas
* pip install torch==1.8.1
* pip install torchtext==0.6.0
* pip install joblib
* pip install rdkit
* pip install scikit-learn
* pip install configargparse
* pip install pyyaml
* pip install pyonmttok
* pip install tensorboard

## On MacOS

* pip install setuptools==59.5.0

## Notes
- Nimom is a generic model for de novo generation.
Nimos is a scaffold-based model for lead optimization by specifying an extra scaffold. They are distinguished by two different motif extraction methods.


- The default task of our code is multi-constraint molecular generation. 
The Wildman–Crippen partition coefficient (logP), drug-likeness (QED) and synthetic accessibility score (SAS) are selected as constraints to train model. 


- Users can customize their own tasks by modifying the code or providing own data.


## Pre-processing 

The pre-processed datasets can be found on the `data/` folder. For raw input file, run: `preprocess.py`
as follows:

```bash
python ./preprocess.py --raw-data data/test/raw.csv --save-path data/test --specific_fragment Nimom
```

**Notes**:
- `-specific_fragment` is required here to specify motif extraction methods. `"Nimom"` or `"Nimos"` is optional.


## building vocabulary

From this configuration, we can build the vocabulary before training the model by run:`build_vocab.py`:
```bash
python build_vocab.py -train_data data/test/dataset.smi -src_vocab data/test/run/test.vocab.src --n_sample -1
```


## Training

Model training can be started by running the `training.py` script:
```bash
python train.py -train_steps 10000 -train_data data/test/train.smi -valid_data data/test/valid.smi -src_vocab data/test/run/test.vocab.src -save_model data/test/run/models/model_lm -tensorboard_log_dir data/test/run/tensorboard
```


## Sampling 

Model sampling use the `generation.py` script:
```bash
python generation.py -model data/coconut_M/run/models/{your model} -src data/coconut_M/lm_input.txt -output data/coconut_M/lm_pred.txt -n_best 5 -beam_size 10
```

