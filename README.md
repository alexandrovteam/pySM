# README #
This repository accompanies the article: *Palmer et al., FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry* and provides a reference implementation of our pipeline for False Discovery Rate-controlled metabolite annotation of high-resolution imaging mass spectrometry data. 

The pipeline is developed and implemented by the [Alexandrov Team](http://www.embl.de/research/units/scb/alexandrov/) at EMBL Heidelberg.

For those wishing to reproduce the results from the paper, please follow the instructions from [Reproducing the results from the paper](#reproducing-the-results-from-the-paper).

For more general use, please consider instructions in [General usage](#general-usage).

## Requirements

* Linux or Mac OS X operating system (was tested on Ubuntu 14.04 and Mac OS X 10.11.3)
* Python 2.7
* `git`
* `wget` (Ubuntu) or `curl` (MacOS) 
* `unzip` (normally provided on both Ubuntu and MacOS)
* 16GB RAM

## Reproducing the results from the paper ##

### Simple installation inside a virtual environment ###
We recommend installing `pySM` and its dependencies inside a virtual environment as follows.

Create a convenient directory, for example `spatial_metabolomics` and clone the repository into there:
```bash
mkdir spatial_metabolomics
cd spatial_metabolomics
git clone https://github.com/alexandrovteam/pySM
```

Next, if you have Anaconda installation of Python, follow the installation instructions [Setting up a virtual environment using conda](#setting-up-a-virtual-environment-using-conda). Otherwise, follow the instructions [Setting up a virtual environment using virtualenv](#setting-up-a-virtual-environment-using-virtualenv).

#### Setting up a virtual environment using `conda`

Initialize and activate an 'pySM' environment with all the dependencies:
```bash
cd pySM
conda env create
source activate pySM
```

Install `pySM` package with `pip`:
```bash
pip install .
```

#### Setting up a virtual environment using `virtualenv`

Setup and activate a new virtual environment:
```bash
pip install virtualenv
virtualenv venv
source venv/bin/activate
```

Install `pySM` and dependencies with `pip`:
```bash
cd pySM
pip install . -r requirements.txt
```


### Evoking FDR-controlled molecular annotation

This section explains how to run FDR-controlled molecular annotation against HMDB for the three MALDI-FTICR-imaging MS datasets from the coronal rat brain sections considered in the paper.

#### Download data
Download and unzip MALDI-imaging MS datasets from the EBI MetaboLights repository as following [for Ubuntu](#for-ubuntu) or [for MacOS](#for-macos):

##### for Ubuntu
```bash
cd pySM/example/datasets
wget -O _RB_a1s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a1s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
wget -O _RB_a2s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
wget -O _RB_a2s2_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s2_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
unzip '_RB_a*.zip'
unzip 'RB_a*.zip'
cd ..
```

##### for MacOS 
```bash
cd pySM/example/datasets
curl -o _RB_a1s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a1s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
curl -o _RB_a2s1_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s1_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
curl -o _RB_a2s2_data.zip http://www.ebi.ac.uk/metabolights/MTBLS313/files/RB_a2s2_data.zip?token=11e11f8d-789d-47e4-83ed-01f4eb768cb6
unzip '_RB_a*.zip'
unzip 'RB_a*.zip'
cd ..
```

#### Download precomputed isotope patterns ####
Download precomputed isotope patterns for sum formulas from HMDB for the considered target and decoy adducts (optional; if not provided will be computed, this computation might take 10 to 20 hours):

```bash
git clone https://github.com/alexandrovteam/precomputed_isotope_patterns 
```

#### Perform annotation ####

Perform molecular annotations for each dataset at the desired FDR level of 0.1 as follows (N.B.: the approximate time needed is 10 hours for the RB_a1s1 dataset, and 5 hours for each of the RB_a2s1 and RB_a2s2 datasets):
```bash
python run_example.py
```

The annotations will be printed to the terminal.





## General usage ##
To process a dataset three things are needed: a high-resolution imaging MS dataset; a metabolite database; and a configuration file

### Dataset ###
* Data should be in the .imzML format. The pipeline is designed for and was tested on **centroided** data. 

### Database ###
The database is a CSV with columns for id, name, exact_mass, formula

### Configuration file ###
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


### Contact
* email: [Andrew Palmer](andrew.palmer@embl.de)
* twitter:  [@alexandrovteam](https://twitter.com/alexandrovteam "alexandrovteam on twitter")
 
### Licence
The source code in this repository is distributed under [the Apache 2.0 licence](https://github.com/alexandrovteam/pySM/blob/master/licence).
