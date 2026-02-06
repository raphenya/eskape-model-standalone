# The ESKAPE Model Standalone

This repository provides a stanalone application to the web version of ESKAPE Model at eskape.mcmaster.ca. The ESKAPE Model is a machine learning-based online resource to facilitate discovery of novel antibiotics against the ESKAPE pathogens, a group of multidrug-resistant bacteria that are responsible for the majority of hospital-acquired infections. 

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

# Installation

The tool requires Python >= 3.10 and conda >= 4.12.0. The latest release can be installed directly from pip or this repository.

# Create a virtual environment

```
python3 -m venv test_env
source test_env/bin/activate
```

# Install eskape_model using tarball

Install the `eskape_model` application within the created `eskape_model` python environment using a tarball.


```
(test_env) amos@Amogelangs-MacBook-Pro % python3 -m pip install /path/to/eskape_model-1.0.0.tar.gz
```

# Dependencies

The following are required dependencies (listed below):

- chemprop version 1.6.1 - https://github.com/chemprop/chemprop.git
- chemfunc version 1.0.10 - https://github.com/swansonk14/chemfunc.git

# Install dependencies

## install chemprop v1.6.1

```
wget https://github.com/chemprop/chemprop/archive/refs/tags/v1.6.1.tar.gz
tar xvf v1.6.1.tar.gz
cd chemprop-1.6.1
pip install -e .
```

## install chemfunc v_1.0.10

```
wget https://github.com/swansonk14/chemfunc/archive/refs/tags/v_1.0.10.tar.gz
tar xvf v_1.0.10.tar.gz
cd chemfunc-v_1.0.10
pip install -e .
```

## cleanup chemprop and chemfunc

```
(test_env) amos@Amogelangs-MacBook-Pro % rm -r chemfunc-v_1.0.10/
(test_env) amos@Amogelangs-MacBook-Pro % rm -r chemprop-1.6.1/
```

## install specific scikit-learn and numpy

```
(test_env) amos@Amogelangs-MacBook-Pro % pip install scikit-learn==1.3.2
(test_env) amos@Amogelangs-MacBook-Pro % pip install numpy==1.26.4
```

## test functions

```
(test_env) amos@Amogelangs-MacBook-Pro % chemprop_predict -h
(test_env) amos@Amogelangs-MacBook-Pro % sklearn_predict -h
(test_env) amos@Amogelangs-MacBook-Pro % chemfunc -h
(test_env) amos@Amogelangs-MacBook-Pro % eskape_model -h
```

# Download ESKAPE model models from eskape.mcmaster.ca or GitHub

```
wget /path/to/models
tar xvf models.1.0.0.tar.gz
mv models db
```

The tree structure of db should look like so:

```
(test_env) amos@Amogelangs-MacBook-Pro db % tree -L 3
.
├── cannonical_data
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
        ├── KP_rf_bak
        ├── PA_chemprop
        ├── PA_rdkit
        ├── PA_rf
        ├── SA_chemprop
        ├── SA_rdkit
        └── SA_rf
```


# run eskape_model

```
(test_env) amos@Amogelangs-MacBook-Pro % eskape_model \
--input_file input.txt \
--output_directory output \
--models_directory db \
--debug > run.log 2>&1 &
```