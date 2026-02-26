# The ESKAPE Model Standalone

This repository provides a standalone application to the web version of ESKAPE Model at eskape.mcmaster.ca. The ESKAPE Model is a machine learning-based online resource to facilitate discovery of novel antibiotics against the ESKAPE pathogens, a group of multidrug-resistant bacteria that are responsible for the majority of hospital-acquired infections. 

The ESKAPE Model predicts the antibacterial activity of inputted molecules against each of the following ESKAPE pathogens:

- EF - **Enterococcus faecium**
- SA - **Staphylococcus aureus**
- KP - **Klebsiella pneumoniae**
- AB - **Acinetobacter baumannii**
- PA - **Pseudomonas aeruginosa**
- BW - **Escherichia coli** (wildtype)
- DKO - **Escherichia coli** (hyperpermeable and efflux deficient)

Models were trained on in-house growth inhibition screening datasets against common laboratory strains of each pathogen. A total of 21 models were trained - three model architectures for each pathogen:

- Random forest using Morgan fingerprints
- Chemprop graph neural network
- Chemprop with RDKit features

# How to Use

### Input:

Molecules are inputted as a CSV file containing SMILES (one per row) with the column heading "smiles". An example csv file with two SMILES (eskape_test_input.csv) is available on this repository.

### Output:

Results are outputted as a TSV file containing the following:
- Prediction scores from each of the 21 models are computed for each molecule. A prediction score is a value between 0 and 1 that denotes how confident the model is that a molecule is antibacterial. Predicted antibacterial molecules will have prediction scores closer to 1, while predicted non-antibacterial molecules will have prediction scores closer to 0.
- For any input compounds that were tested against the ESKAPE pathogens during training data acquisition, this tool will additionally output the experimental optical density (OD) values in the "validated" row. OD is a measure of bacterial cell growth, where a high OD means the bacteria grew in the presence of the compound, and a low OD means the compound was able to inhibit the growth of the bacteria. For reference, an OD less than 0.06 denotes full growth inhibition. All OD values were normalized by plate based on the interquartile mean.
- Several metrics are also calculated for each compound:
  - Sum of PS: Sum of prediction scores from all pathogen models for one compound. This metric can be used to prioritize broad-spectrum antibacterial compounds.
  - PPF: The ratio of the highest prediction score for a compound (PS1) to the second highest (PS2). This metric can be used to prioritize pathogen-prioritized antibacterial compounds.
  - Molecular weight: Size of the molecule in g/mol
  - clogP: Calculated octanol-water partition coefficient, where high clogP values mean the compound is more lipophilic. clogP is an important metric for solubility and bioavailability.
  - TNN: The TNN similarity measures the structural similarity (value between 0-1) of an input molecule to the most similar molecule (nearest neighbour) from the training set. TNN similarity closer to 1 indicates the molecules are more similar (TNN similarity = 1 means the molecules are equal). Predictions on compounds that are more similar to the training set are likely to be more accurate. Nearest neighbor SMILES from the training set are included in the TSV.

### Interpretation:

While all models were trained on the same datasets using the same training scheme, the three model types differ in terms of architecture and molecular representation. Prediction scores for the same molecule and pathogen will therefore vary based on the model type. Note that prediction scores do not correlate directly with likelihood of activity or potency, but rather represent model confidence.

###Runtime:

Note: Predictions on 1 molecule takes ~2 minutes. Predictions on 100 molecules takes ~3.5 minutes.

# Installation

The tool requires Python 3.10. Python versions more recent than 3.10 have been tested and do not work. Installation takes ~5 minutes.

# Create a virtual environment

```
python3 -m venv eskape_env
source eskape_env/bin/activate
```

# Install eskape_model using pip

The latest release can be installed directly from pip or this repository which will also install the dependencies `chemprop` and `chemfunc`.

```
pip install eskape_model
```

Or

# Install eskape_model using tarball

Install the `eskape_model` application within the created `eskape_model` python environment using a tarball.


```
(eskape_env) amos@Amogelangs-MacBook-Pro % python3 -m pip install /path/to/eskape_model-1.0.0.tar.gz
```

# Dependencies

The following are required dependencies (listed below):

- chemprop version 1.6.1 - https://github.com/chemprop/chemprop.git
- chemfunc version 1.0.10 - https://github.com/swansonk14/chemfunc.git

# Install dependencies

## install chemprop v1.6.1

```
wget https://github.com/chemprop/chemprop/archive/refs/tags/v1.6.1.tar.gz
python3 -m pip install v1.6.1.tar.gz
```

## install chemfunc v_1.0.10

```
wget https://github.com/swansonk14/chemfunc/archive/refs/tags/v_1.0.10.tar.gz
python3 -m pip install v_1.0.10.tar.gz
```

## install specific scikit-learn and numpy

```
(eskape_env) amos@Amogelangs-MacBook-Pro % pip install scikit-learn==1.3.2
(eskape_env) amos@Amogelangs-MacBook-Pro % pip install numpy==1.26.4
```

## test functions

```
(eskape_env) amos@Amogelangs-MacBook-Pro % chemprop_predict -h
(eskape_env) amos@Amogelangs-MacBook-Pro % sklearn_predict -h
(eskape_env) amos@Amogelangs-MacBook-Pro % chemfunc -h
(eskape_env) amos@Amogelangs-MacBook-Pro % eskape_model -h
```

# Download ESKAPE model models from eskape.mcmaster.ca or GitHub

Please download the models and training data at [GitHub](https://github.com/arnolp3/The-ESKAPE-Model).

Create a directory `db` with two sub-directories `canonical_data` and `models`. From the downloaded models data, add `training_data_canonical.csv` to `db/canonical_data/` directory. Add all models to directory `db/models/all/`.

The tree structure of db should look like so:

```
(eskape_env) amos@Amogelangs-MacBook-Pro db % tree -L 3
.
├── canonical_data
│   └── training_data_canonical.csv
└── models
    └── all
        ├── AB_chemprop
        ├── AB_rdkit
        ├── AB_rf
        ├── BW_chemprop
        ├── BW_rdkit
        ├── BW_rf
        ├── DKO_chemprop
        ├── DKO_rdkit
        ├── DKO_rf
        ├── EF_chemprop
        ├── EF_rdkit
        ├── EF_rf
        ├── KP_chemprop
        ├── KP_rdkit
        ├── KP_rf
        ├── PA_chemprop
        ├── PA_rdkit
        ├── PA_rf
        ├── SA_chemprop
        ├── SA_rdkit
        └── SA_rf
```


# run eskape_model

```
(eskape_env) amos@Amogelangs-MacBook-Pro % eskape_model \
--input_file input.txt \
--output_directory output \
--models_directory db \
--debug > run.log 2>&1 &
```
