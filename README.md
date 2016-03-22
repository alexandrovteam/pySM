# README #
This repository accompanies the article: *Palmer et al., FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry*.

## What is this repository for? ##
This is a reference implementation of our pipeline for false-discovery-rate controlled annotation of high-resolution imaging mass spectrometry data. The pipeline is developed and implemented by the [Alexandrov Team](http://www.embl.de/research/units/scb/alexandrov/) at EMBL Heidelberg.

The code was tested on Ubuntu 14.04 and Mac OS X 10.11.3.

## Requirements

* Python 2.7
* Linux or Mac OS X operating system
* git
* wget
* unzip
* 16 GB RAM

## Installation ##

We recommend installing inside a virtual environment

Create a convenient directory and clone the repository
```bash
mkdir spatial_metabolomics
cd spatial_metabolomics
git clone https://github.com/alexandrovteam/pySM
```

If you have Anaconda installation of Python, use [conda](#conda) installation instructions.
Otherwise, use [virtualenv](#virtualenv).

### virtualenv

Setup and activate a new virtual environment
```bash
pip install virtualenv
virtualenv venv
source venv/bin/activate
```

Install `pySM` and dependencies with `pip`
```bash
cd pySM
pip install . -r requirements.txt
```

### conda

Initialize and activate `pySM` environment with all the dependencies
```bash
cd pySM
conda env create
source activate pySM
```

Install `pySM` package with `pip`
```bash
pip install .
```

## Reproducing the results from the paper ##

Download and unzip MALDI imaging MS datasets from the EBI MetaboLights repository.

```bash
cd pySM/example/datasets
wget -O _RB_a1s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a1s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
wget -O _RB_a2s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
wget -O _RB_a2s2_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s2_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
unzip '_RB_a*.zip'
unzip 'RB_a*.zip'
```

Download precomputed isotope patterns (optional; if not provided will be computed, it might take 10 to 20 hours).

```bash
git clone https://github.com/alexandrovteam/precomputed_isotope_patterns pySM/example/precomputed_isotope_patterns
```
Run script to produce molecular annotations for each dataset at the desired FDR equal 0.1:
```bash
cd pySM/example
python run_example.py
```

The annotations will be printed to the terminal.

## Processing a dataset ##
To process a dataset three things are needed: a high-resolution imaging MS dataset; a metabolite database; and a configuration file

### 1. Dataset ###
* Data should be in the .imzML format
    * The pipeline is currently only designed for centroid data. Using centroided data is **highly** recommended for run time and annotation performance.

### 2. Database ###
The database is a csv with columns for id, name, exact_mass, formula

### 3. Configuration file ###
A complete example configuration can be found [here](https://github.com/alexandrovteam/pySM/blob/master/pySM/example/example_config.json).
The following parameters should be set individually for every dataset, other parameters can generally be left at their default values

```json
   "name":"dataset_short_name",
   "image_generation":{
      "ppm":4.0,
      },
   "file_inputs":{
     "data_file":"/path/to/imaging_ms_dataset.imzML",
     "database_load_folder":"/path/to/tmp_folder_for_storing_isotope_patterns",
     "results_folder":"/path/to/folder_for_storing_results",
     "database_file":"/path/to/database.csv"
   },
   "fdr":{
      "pl_adducts":[
         {"adduct":"+H"},
         {"adduct":"+Na"},
         {"adduct":"+K"}
         ],
   },
   "isotope_generation":{
      "charge":[
        {"polarity":"+", "n_charges":1}
        ],
      "isocalc_sig":0.01,
```
* name: a short name for the dataset, if ```"name":"" ``` the imzml filename will be used
* ppm: the m/z window for ion images
* file_inputs: path for loading data/storing results
* fdr: false discovery rate properties
  * pl_adducts: real adducts to search for
* isotope_generation:
  * charge: polarity and charge state to search for (the pipeline currently only supports one charge state at a time). e.g. for negative mode singly charged use ```"charge":[
        {"polarity":"-", "n_charges":1}
        ],```
  * isocalc_sig: peaks are predicted with a gaussian shape. This parameter is the sigma parameter. sigma = FWHM/2.3548.
  * isocalc_resoultion is *not* mass spectral resolution, it is the digitisation rate of the isotope patterns


### Who do I talk to? ###
If you are having difficulty running the pipeline, please get in touch with [Andrew Palmer](andrew.palmer@embl.de).
