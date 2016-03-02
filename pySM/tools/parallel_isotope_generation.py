__author__ = 'palmer'
from pySM import spatial_metabolomics
from pySM.parse_databases import parse_databases

import os
from joblib import Parallel, delayed

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
    mz_list_tmp = spatial_metabolomics.calculate_isotope_patterns(sum_formulae, adduct=adduct, isocalc_sig=isocalc_sig,
                                                    isocalc_resolution=isocalc_resolution, charge=charge, verbose=True)
    spatial_metabolomics.save_pattern(load_file, db_dump_folder, mz_list_tmp)
    return True


if __name__== '__main__':
    json_filename = "/Users/palmer/Documents/Projects/2016/SM_Apps/RatBrain/ChEBI/RBa2s1.json"
    config = spatial_metabolomics.get_variables(json_filename)
    db_filename = config['file_inputs']['database_file']
    sum_formulae = parse_databases.read_generic_csv(db_filename)
    adducts = [a['adduct'] for a in config['isotope_generation']['adducts']]
    if '' in sum_formulae:
        del sum_formulae['']
    Parallel(n_jobs=4)(delayed(generate_isotope_pattern)(sum_formulae,a,config) for a in adducts)
