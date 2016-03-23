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

# Reproducing the results from the paper #

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

Perform molecular annotations for each dataset at the desired FDR level of 0.1 as follows (N.B.: run time is approximately 5 hours for the RB_a1s1 dataset, and 1 hour for each of the RB_a2s1 and RB_a2s2 datasets on a MacBook Pro 3GHz i7 with 16GB DDR3):
```bash
python run_example.py
```

After all datasets have been processed a table of molecular formulas annotated for each dataset at an FDR of 0.1 will be displayed on the command line.  Results are displayed as comma separated values, showing the molecular formula-adduct, the values for each of the individual measures considered in the manuscript (&rho; chaos &rho; spatial, &rho; spectral) and the composite MSM score.   
e.g.

```
Pass formula for RBa2s2
sum formula,adduct,p_chaos,p_spatial,p_spectral,msm
C24H38O4,+K,0.998696103579,0.962699782944,0.965102356221,0.927892373696
C24H50NO7P,+K,0.999294484149,0.791215125324,0.980649030323,0.775356932611
C26H42O4,+K,0.998390741419,0.954013608755,0.961014346342,0.915345362936
C26H54NO7P,+K,0.999508321881,0.778993545242,0.966106492883,0.752220689582
C27H46O,+K,0.999367695796,0.922418016749,0.976702607076,0.900358421159
...
...
...
```

The corresponding table showing the scores for every molecule/adduct from the database (including decoys) `pySM\example\RB_x\RB_x_spatial_all_adducts_full_results.txt`



# General usage #
The pipeline has two core parts: 
1. Calculation of Metabolite Signal Match (MSM) scores for every molecular formula in a metabolite database
2. Reporting of molecular formula at a specified FDR 

## Installation ##
See [Simple installation inside a virtual environment](#simple-installation-inside-a-virtual-environment)

## Annotating a dataset ##
### Inputs ###
To process a dataset three things are needed: 
1. a high-resolution imaging MS dataset; 
2. a metabolite database
3. a configuration file

#### Dataset ####
Data should be in the .imzML format. The pipeline is designed for and was tested on **centroided** data. 

#### Database ####
The database is a CSV with columns for id, name, exact_mass, formula

#### Configuration file ####
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

### Calculating MSM Scores ###
The `spatial_metabolomics` module runs the pipeline for calculating MSM scores. To calculate MSM scores for a whole dataset and database simply pass the configuration file to the `run_pipeline` function:
```python
from pySM import spatial_metabolomics
json_filename = '/path/to/config.json'
spatial_metabolomics.run_pipeline(json_filename)
```
This will then write the MSM score for every combination of molecular formula and target adduct found in the metabolite database to a text file in the `"results_folder"` specified in the config file. Additionally a randomly selected set of decoy adducts will be chosen for , and their MSM scores calculated. The number of decoy adducts is controlled by the config parameter `fdr\n_im`. 

### Reporting FDR ###
The main use of FDR control is to report which molecular formulas are annotated at a fixed FDR. This uses the results file generated by `spatial_metabolomics.run_pipeline` and the target and decoy adducts specified in the configuration file.
```python
from  pySM import spatial_metabolomics, fdr_measures
json_filename = '/path/to/config.json'
results_fname = spatial_metabolomics.generate_output_filename(spatial_metabolomics.get_variables(json_filename),[],'spatial_all_adducts')
target_adducts,decoy_adducts = fdr_measures.get_adducts(json_filename)
fdr = fdr_measures.decoy_adducts(results_fname,target_adducts,decoy_adducts)
fdr_target=0.1
print fdr.decoy_adducts_get_pass_list(fdr_target,n_reps=20,col='msm')
```
This will print a list of molecular-formula for each target adduct that have an MSM score which results in an FDR of less than `fdr_target`. 
## Contact
* email: [Andrew Palmer](andrew.palmer@embl.de)
* twitter:  [@alexandrovteam](https://twitter.com/alexandrovteam "alexandrovteam on twitter")
 
## Licence
The source code in this repository is distributed under [the Apache 2.0 licence](https://github.com/alexandrovteam/pySM/blob/master/licence).
