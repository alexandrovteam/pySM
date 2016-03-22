# README #
This repository accompanies the article: *Palmer et al., FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry*. 


## What is this repository for? ##
This is a reference implementation of our pipeline for false-discovery-rate controlled annotation of high-resolution imaging mass spectrometry data. The pipeline is developed and implemented by the [Alexandrov Team](http://www.embl.de/research/units/scb/alexandrov/) at EMBL Heidelberg.

## Installation ##

* Tested on Ubuntu 14, Mac OS-X 10.11.3
* Python 2.7
* 16 GB RAM (depending on dataset size)
We recommend installing inside a virtual environment 

Create a convenient directory clone the repository

```
mkdir spatial_metabolomics
cd spatial_metabolomics
git clone https://github.com/alexandrovteam/pySM
```

### virtualenv

Setup and activate a new virtual environment
```
virtualenv venv
source venv/bin/activate
```

Install `pySM` and dependencies with `pip`
``` 
cd pySM
pip install . -r requirements.txt
```

### conda

Initialize and activate `pySM` environment with all the dependencies
```
cd pySM
conda env create
source activate pySM
```

Install `pySM` package with `pip`
```
pip install .
```

## Reproducing the Results from the Paper ##
1.  Download files
  * Data: 
    * Go to the EBI MetaboLights repository ([MTBLS313](http://www.ebi.ac.uk/metabolights/MTBLS313) (under embargo during review, for early access contact [Andrew Palmer](andrew.palmer@embl.de)))
    * From 'study files': download RB_a1s1_data.zip, RB_a2s1_data.zip, RB_a2s2_data.zip
    * Unzip so all RB_x.imzml and RB_x.ibd files are in the ```/pySM/example/datasets``` folder
  * Isotope patterns (optional)
    * generating isotope patterns can be time-consuming, we provide pre-computed patterns [here](https://github.com/alexandrovteam/precompiled_isotope_patterns/archive/master.zip). Unzip into ```/pySM/example/precomputed_patterns```
2. Run script to produce molecular annotations at 10% FDR: 
  * run script from example directory
```
cd pySM/example
python run_example.py
```
  * The annotations will be printed to the command line

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
