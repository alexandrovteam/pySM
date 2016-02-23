# README #

This repository accompanies the article: 
FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry
*Andrew Palmer, Michael Becker, Ilya Chernyavsky, Dominik Fay, Artem Tarasov, Vitaly Kovalev, Jens Fuchser, Sergey Nikolenko, Theodore Alexandrov*


### What is this repository for? ###

This is a reference implementation of our pipeline for false-discovery-rate controlled annotation of high-resolution imaging mass spectrometry data. 


### How do I get set up? ###

* The following python dependancies can be installed with pip
    * pyimzml
    * numpy
    * scipy
* The following dependancies must be cloned from their github repositories
    * https://github.com/alexandrovteam/pySM
    * https://github.com/alexandrovteam/pyMS
    * https://github.com/alexandrovteam/pyIMS
* remember to add the clone directories to the python path
    '''python
    import sys
    sys.path.append(../clone/directory/)
    '''


### Processing a dataset ###
* Data should be in .imzml format
    * The pipeline is currently only designed for centroid data. Using centroided data is **highly** recommended for run time and annotation performance
* Instructions for running the pipeline on the data accompanying the paper can be found in the example folder

### Who do I talk to? ###
If you're having difficulty running the pipeline get in touch with palmer@embl.de
http://www.embl.de/research/units/scb/alexandrov/
