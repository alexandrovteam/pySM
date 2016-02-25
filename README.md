# README #

This repository accompanies the article: 
FDR-controlled metabolite annotation for high-resolution imaging mass spectrometry
*Andrew Palmer, Michael Becker, Ilya Chernyavsky, Dominik Fay, Artem Tarasov, Vitaly Kovalev, Jens Fuchser, Sergey Nikolenko, Theodore Alexandrov*


### What is this repository for? ###

This is a reference implementation of our pipeline for false-discovery-rate controlled annotation of high-resolution imaging mass spectrometry data. 


### How do I get set up? ###

* The following python dependancies can be installed from pypi with pip
    *  ```
       sudo pip install numpy
      ```
    * ```
       sudo pip install scipy
      ```
   *  ```
       sudo pip install pyimzml
      ```
* The following libraries can be installed directly from git
    * ```
       sudo pip install git+https://github.com/alexandrovteam/pySM
      ```
    * ```
       sudo pip install git+https://github.com/alexandrovteam/pyMS@adp
      ```
    * ```
       sudo pip install git+https://github.com/alexandrovteam/pyIMS@adp
      ```
* remember to add the clone directories to the python path
    '''python
    import sys
    sys.path.append(../clone/directory/)
    '''


### Processing a dataset ###
* Data should be in .imzml format
    * The pipeline is currently only designed for centroid data. Using centroided data is **highly** recommended for run time and annotation performance
* An example for running the pipeline on the data accompanying the paper can be found [here](https://github.com/alexandrovteam/pySM/tree/master/pySM/example)

### Who do I talk to? ###
If you're having difficulty running the pipeline get in touch with palmer@embl.de
http://www.embl.de/research/units/scb/alexandrov/
