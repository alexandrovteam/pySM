# README #
This repository accompanies the article: *Palmer et al., FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry*. 


## What is this repository for? ##
This is a reference implementation of our pipeline for false-discovery-rate controlled annotation of high-resolution imaging mass spectrometry data. The pipeline is developed and implemented by the [Alexandrov Team](http://www.embl.de/research/units/scb/alexandrov/) at EMBL Heidelberg.

## Installation ##

* Tested on Ubuntu 14, Mac OS-X
* Python 2.7
* 16 GB RAM (depending on dataset size)

Install with pip:
``` sudo pip install git+https://github.com/alexandrovteam/pySM ```

## Reproducing the Results from the Paper ##
The annotations reported in the paper can be reproduced using the scips in the [example folder](https://github.com/alexandrovteam/pySM/blob/master/pySM/example/) of this repository 
* [Clone](https://github.com/alexandrovteam/pySM.git) or [download and unizip](https://github.com/alexandrovteam/pySM/archive/master.zip) this repository 
* Download the imagingMS datasets from the EBI MetaboLights repositories ([MTBLS313](http://www.ebi.ac.uk/metabolights/MTBLS313), [MTBLS317](http://www.ebi.ac.uk/metabolights/MTBLS317) (under embargo during review, for early access please contact [Andrew Palmer](andrew.palmer@embl.de)))
* Unzip the files in config.zip files within the [example folder](https://github.com/alexandrovteam/pySM/blob/master/pySM/example/config.zip)
* Edit the following paths within'file_inputs'section of the config:
  ```
  "data_file":"/path/to/imaging_ms_dataset.imzML"
  "database_load_folder":"/path/to/tmp_folder_for_storing_isotope_patterns
  "results_folder":"/path/to/folder_for_storing_results
  "database_file":"/path/to/database.csv" 
  ```
* From the command line, run the python script for reproducing the annotations reported in the manuscript [link](https://github.com/alexandrovteam/pySM/blob/master/pySM/example/run_example.py) 
```
>> python run_example.py
```


## Processing a dataset ##
To process a dataset three things are needed:
1. a high-resolution imaging MS dataset
2. a metabolite database
3. a configuration file 

### Dataset ###
* Data should be in the .imzML format
    * The pipeline is currently only designed for centroid data. Using centroided data is **highly** recommended for run time and annotation performance.

### Database ###
The database is a csv with columns for database_id, name, molecular_formula

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
  * isocalc_sig: peaks are predicted with a gaussian shape. This parameter is the sigma parameter. sigma = FWHM/2.3548 or 


### Who do I talk to? ###
If you are having difficulty running the pipeline, please get in touch with [Andrew Palmer](andrew.palmer@embl.de).
