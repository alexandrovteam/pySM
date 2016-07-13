__author__ = 'palmer'
from pySM import spatial_metabolomics
from pySM.parse_databases import parse_databases
from pyMSpec.pyisocalc import pyisocalc
from cpyMSpec.legacy_interface import complete_isodist
import numpy as np
import os
from joblib import Parallel, delayed
def calculate_isotope_pattern( sum_formula, adduct,isocalc_sig, isocalc_resolution, charge, verbose, instrument):
    sf_str ="{}{}".format(sum_formula,adduct)
    if verbose:
        print sf_str
    try:
        #sf = pyisocalc.parseSumFormula(sf_str)
        isotope_ms = instrument.get_isotope_pattern(sum_formula+adduct, charge=charge)
    except pyisocalc.ParseError as e:
        print "error->", str(e), sum_formula, adduct
        return ["foo",[]]
    except KeyError as e:
        print  str(e)
        return ["foo",[]]
    except pyisocalc.InvalidFormulaError as e:
        print  str(e)
        return ["foo",[]]
    #try:
        #isotope_ms = complete_isodist(sf,sigma=isocalc_sig,charge=charge, pts_per_mz=isocalc_resolution)
    except MemoryError as e:
        print sum_formula, adduct, str(e)
        return ["foo",[]]
    except KeyError as e:
        print sum_formula, adduct, str(e)
        return ["foo",[]]
    #mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')
    return [sum_formula,isotope_ms.get_spectrum(source='centroids')]

def calculate_isotope_patterns(sum_formulae, adduct='', isocalc_sig=0.01, isocalc_resolution=200000.,
                                   isocalc_do_centroid=True, charge=1,verbose=True, instrument=[]):
    ### Generate a mz list of peak centroids for each sum formula with the given adduct
    mz_list = {}
    n_parts =5
    #for ii in range(int(n_parts)):
     #   ix_start = int(float(ii)*len(sum_formulae)/n_parts)
      #  ix_end = int(min((float(ii+1)*len(sum_formulae)/n_parts, len(sum_formulae))))
       # sum_formulae_short = sum_formulae.keys()[ix_start:ix_end]
    mz_list_tmp = Parallel(n_jobs=6)(delayed(calculate_isotope_pattern)(sf, adduct,isocalc_sig, isocalc_resolution, charge, verbose, instrument
                                      ) for sf in sum_formulae)
        #mz_list_tmp = map(lambda sf:calculate_isotope_pattern(sf, adduct,isocalc_sig, isocalc_resolution, charge, verbose), sum_formulae)
        #print "mz_list", len(mz_list_tmp), len(sum_formulae)
    for row in mz_list_tmp:
        mz_list[row[0]]={}
        mz_list[row[0]][adduct]=row[1]
    return mz_list

def generate_isotope_pattern(sum_formulae, adduct, config, db_filename):
    # Extract variables from config dict
    from pySM.spatial_metabolomics import get_instrument
    instrument = get_instrument(config["isotope_generation"]["instrument"], config["isotope_generation"]["resolving_power"], config["isotope_generation"]["at_mz"])
    isocalc_sig = np.round(instrument.sigma_at_mz(200), decimals=4)
    isocalc_resolution = instrument.points_per_mz(isocalc_sig)
    db_dump_folder = config['file_inputs']['database_load_folder']
    if len(config['isotope_generation']['charge']) > 1:
        print 'Warning: only first charge state currently accepted'
    charge = int('{}{}'.format(config['isotope_generation']['charge'][0]['polarity'],
                               config['isotope_generation']['charge'][0][
                                   'n_charges']))  # currently only supports first charge!!
    db_name = os.path.splitext(os.path.basename(db_filename))[0]
    load_file = '{}/{}_{}_{}_{}_{}.dbasedump'.format(db_dump_folder, db_name, instrument.__class__.__name__, adduct, isocalc_sig, isocalc_resolution)
    print "{} -> generating".format(load_file)

    mz_list_tmp = calculate_isotope_patterns(sum_formulae, adduct=adduct, isocalc_sig=isocalc_sig,
                                                    isocalc_resolution=isocalc_resolution, charge=charge, verbose=False, instrument=instrument)
    spatial_metabolomics.save_pattern(load_file, db_dump_folder, mz_list_tmp)
    return True


if __name__== '__main__':
    json_filename = "/home/palmer/Documents/tmp_data/MB/URenn/16135s1-3_mousebrain_dhbsub.json"
    config = spatial_metabolomics.get_variables(json_filename)
    adducts = [a['adduct'] for a in config['isotope_generation']['adducts']]
    db_filenames = config['file_inputs']['database_file']
    if isinstance(db_filenames, basestring):
        db_filenames = [db_filenames,]
    for db_filename in db_filenames:
        sum_formulae = parse_databases.read_generic_csv(db_filename,header=1,idcol=0,namecol=1,sfcol=2, sep='\t')
        if '' in sum_formulae:
            del sum_formulae['']
        for a in adducts:
            generate_isotope_pattern(sum_formulae,a,config, db_filename)
