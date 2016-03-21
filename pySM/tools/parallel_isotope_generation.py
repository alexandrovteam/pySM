__author__ = 'palmer'
from pySM import spatial_metabolomics
from pySM.parse_databases import parse_databases
from pyMSpec.pyisocalc import pyisocalc
import numpy as np
import os
from joblib import Parallel, delayed

def calculate_isotope_pattern( sum_formula, adduct,isocalc_sig, isocalc_resolution, charge, verbose):
    sf_str ="{}{}".format(sum_formula,adduct)
    if verbose:
        print sf_str
    try:
        sf = pyisocalc.parseSumFormula(sf_str)
    except pyisocalc.ParseError as e:
        print "error->", str(e), sum_formula, adduct
        return ["foo",[]]
    except KeyError as e:
        print  str(e)
        return ["foo",[]]
    except pyisocalc.InvalidFormulaError as e:
        print  str(e)
        return ["foo",[]]
    if verbose:
        print sf
    try:
        isotope_ms = pyisocalc.complete_isodist(sf,sigma=isocalc_sig,charge=charge, pts_per_mz=isocalc_resolution)
    except MemoryError as e:
        print sf, str(e)
        return ["foo",[]]
    except KeyError as e:
        print sf, str(e)
        return ["foo",[]]
    #mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')
    return [sum_formula,isotope_ms.get_spectrum(source='centroids')]

def calculate_isotope_patterns(sum_formulae, adduct='', isocalc_sig=0.01, isocalc_resolution=200000.,
                                   isocalc_do_centroid=True, charge=1,verbose=True):
    ### Generate a mz list of peak centroids for each sum formula with the given adduct
    mz_list = {}
    n_parts =5
    #for ii in range(int(n_parts)):
     #   ix_start = int(float(ii)*len(sum_formulae)/n_parts)
      #  ix_end = int(min((float(ii+1)*len(sum_formulae)/n_parts, len(sum_formulae))))
       # sum_formulae_short = sum_formulae.keys()[ix_start:ix_end]
    mz_list_tmp = Parallel(n_jobs=6)(delayed(calculate_isotope_pattern)(sf, adduct,isocalc_sig, isocalc_resolution, charge, verbose
                                      ) for sf in sum_formulae)
        #mz_list_tmp = map(lambda sf:calculate_isotope_pattern(sf, adduct,isocalc_sig, isocalc_resolution, charge, verbose), sum_formulae)
        #print "mz_list", len(mz_list_tmp), len(sum_formulae)
    for row in mz_list_tmp:
        mz_list[row[0]]={}
        mz_list[row[0]][adduct]=row[1]
    return mz_list

def generate_isotope_pattern(sf_list,adduct, config):
    # Extract variables from config dict
    db_dump_folder = config['file_inputs']['database_load_folder']
    isocalc_sig = float(config['isotope_generation']['isocalc_sig'])
    isocalc_resolution = float(config['isotope_generation']['isocalc_resolution'])
    if len(config['isotope_generation']['charge']) > 1:
        print 'Warning: only first charge state currently accepted'
    charge = int('{}{}'.format(config['isotope_generation']['charge'][0]['polarity'],
                               config['isotope_generation']['charge'][0][
                                   'n_charges']))  # currently only supports first charge!!
    db_name = os.path.splitext(os.path.basename(db_filename))[0]
    load_file = '{}/{}_{}_{}_{}.dbasedump'.format(db_dump_folder, db_name, adduct, isocalc_sig, isocalc_resolution)
    print "{} -> generating".format(load_file)
    mz_list_tmp = calculate_isotope_patterns(sum_formulae, adduct=adduct, isocalc_sig=isocalc_sig,
                                                    isocalc_resolution=isocalc_resolution, charge=charge, verbose=False)
    spatial_metabolomics.save_pattern(load_file, db_dump_folder, mz_list_tmp)
    return True


if __name__== '__main__':
    json_filename = "/home/palmer/Documents/tmp_data/parallel_chebi/RBa2s1.json"
    config = spatial_metabolomics.get_variables(json_filename)
    db_filename = config['file_inputs']['database_file']
    sum_formulae = parse_databases.read_generic_csv(db_filename)
    adducts = [a['adduct'] for a in config['isotope_generation']['adducts']]
    if '' in sum_formulae:
        del sum_formulae['']
    for a in adducts:
        generate_isotope_pattern(sum_formulae,a,config)
