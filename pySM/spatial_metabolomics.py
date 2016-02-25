import numpy as np
import sys
import os
sys.path.append('/Users/palmer/Documents/python_codebase/')


def d_update(d, u):
    import collections
    #http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = d_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


def get_variables(json_filename):
    import json
    config = json.loads(open(json_filename).read())
    # maintain compatibility with previous versions
    # defaults are be how everything *should* be -> update makes sure that whatever is loaded conforms to this
    compatability_defaults = {'image_generation':{'smooth': True, 'smooth_params':{} } }
    config = d_update(compatability_defaults,config)
    return compatability_defaults


### We simulate a mass spectrum for each sum formula/adduct combination. This generates a set of isotope patterns (see http://www.mi.fu-berlin.de/wiki/pub/ABI/QuantProtP4/isotope-distribution.pdf) which can provide additional informaiton on the molecule detected. This gives us a list of m/z centres for the molecule
def calculate_isotope_patterns(sum_formulae, adduct='', isocalc_sig=0.01, isocalc_resolution=200000.,
                                   isocalc_do_centroid=True, charge=1,verbose=True):
    from pyMS.pyisocalc import pyisocalc
    ### Generate a mz list of peak centroids for each sum formula with the given adduct
    mz_list = {}
    for n, sum_formula in enumerate(sum_formulae):
        try:
            if verbose:
                print sum_formula, adduct
            sf = pyisocalc.complex_to_simple(sum_formula+adduct)
        except KeyError as e:
            if str(e).startswith("KeyError: "):
                print str(e)
                continue
        except ValueError as e:
            if str(e).startswith("Element not recognised"):
                print str(e)
                continue
            else:
                print sum_formula, adduct
                raise
        if sf == None: #negative atoms as a result of simplification
            print 'negative adduct for {} : {}'.format(sum_formula,adduct)
            continue
        try:
            isotope_ms = pyisocalc.isodist(sf, plot=False, sigma=isocalc_sig, charges=charge,
                                       resolution=isocalc_resolution, do_centroid=isocalc_do_centroid)
        except KeyError as e:
            if str(e).startswith("KeyError: "):
                print str(e)
                continue

        if not sum_formula in mz_list:
            mz_list[sum_formula] = {}
        mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')

    return mz_list


def generate_isotope_patterns(config,verbose=True):
    from pySpatialMetabolomics.parse_databases import parse_databases

    import pickle
    # Extract variables from config dict
    db_filename = config['file_inputs']['database_file']
    db_dump_folder = config['file_inputs']['database_load_folder']
    isocalc_sig = float(config['isotope_generation']['isocalc_sig'])
    isocalc_resolution = float(config['isotope_generation']['isocalc_resolution'])
    if len(config['isotope_generation']['charge']) > 1:
        print 'Warning: only first charge state currently accepted'
    charge = int('{}{}'.format(config['isotope_generation']['charge'][0]['polarity'],
                               config['isotope_generation']['charge'][0][
                                   'n_charges']))  # currently only supports first charge!!
    adducts = [a['adduct'] for a in config['isotope_generation']['adducts']]

    # Read in molecules
    sum_formulae = parse_databases.read_generic_csv(db_filename)
    if '' in sum_formulae:
        if verbose:
            print 'empty sf removed from list'
        del sum_formulae['']
    # Check if already genrated and load if possible, otherwise calculate fresh   
    db_name = os.path.splitext(os.path.basename(db_filename))[0]
    mz_list = {}
    for adduct in adducts:
        load_file = '{}/{}_{}_{}_{}.dbasedump'.format(db_dump_folder, db_name, adduct, isocalc_sig, isocalc_resolution)
        if os.path.isfile(load_file):
            if verbose:
                print  "{} -> loading".format(load_file)
            mz_list_tmp = pickle.load(open(load_file, 'r'))
        else:
            if verbose:
                print "{} -> generating".format(load_file)
            mz_list_tmp = calculate_isotope_patterns(sum_formulae, adduct=adduct, isocalc_sig=isocalc_sig,
                                                     isocalc_resolution=isocalc_resolution, charge=charge)
            if db_dump_folder != "":
                pickle.dump(mz_list_tmp, open(load_file, 'w'))
        # add patterns to total list
        for sum_formula in sum_formulae:
            if sum_formula not in mz_list_tmp:# could be missing if [M-a] would have negative atoms
                continue
            if sum_formula not in mz_list:
                mz_list[sum_formula]={}
            ## this limit of 4 is hardcoded to reduce the number of calculations
            n = np.min([4,len(mz_list_tmp[sum_formula][adduct][0])])
            mz_list[sum_formula][adduct] = [mz_list_tmp[sum_formula][adduct][0][0:n],mz_list_tmp[sum_formula][adduct][1][0:n]]
    if verbose:
        print  'all isotope patterns generated and loaded'
    return sum_formulae, adducts, mz_list


def hot_spot_removal(xics, q):
    print 'moved to pyIMS.smoothing - should be called'
    for xic in xics:
        xic_q = np.percentile(xic, q)
        xic[xic > xic_q] = xic_q
    return xics

def apply_image_processing(config, ion_datacube):
    """
    Function to apply pre-defined image processing methods to ion_datacube
    #todo: expose parameters in config
    :param ion_datacube:
        object from pyIMS.ion_datacube already containing images
    :return:
        ion_datacube is updated in place.
        None returned
    """
    from pyIMS import smoothing
    #todo hot_spot_removal shouldn't be separately coded - should be within smooth_methods of config and iterated over
    # every method in smoothing should accept (im,**args)
    q = config['image_generation']['q']
    if q > 0:
        for xic in ion_datacube.xic:
            smoothing.hot_spot_removal(xic, q) #updated in place
    smooth_method = config['image_generation']['smooth']
    smooth_params = config['image_generation']['smooth_params']
    if smooth_method != '':
        for ii in range(ion_datacube.n_im):
            im = ion_datacube.xic_to_image(ii)
            #todo: for method in smoothing_methods:
            methodToCall = getattr(smoothing,smooth_method)
            im_s = methodToCall(im,**smooth_params)
            ion_datacube.xic[ii] = ion_datacube.image_to_xic(im_s)
    return None

def run_search(config, IMS_dataset, sum_formulae, adducts, mz_list):
    from pyIMS.image_measures import level_sets_measure, isotope_image_correlation, isotope_pattern_match
    import time
    ### Runs the main pipeline
    # Get sum formula and predicted m/z peaks for molecules in database
    ppm = config['image_generation']['ppm']  # parts per million -  a measure of how accuracte the mass spectrometer is
    nlevels = config['image_generation']['nlevels']  # parameter for measure of chaos
    do_preprocessing = config['image_generation']['do_preprocessing']
    interp = config['image_generation']['smooth']
    measure_value_score = {}
    iso_correlation_score = {}
    iso_ratio_score = {}
    t0 = time.time()
    t_el = 0
    for adduct in adducts:
        print 'searching -> {}'.format(adduct)
        for ii,sum_formula in enumerate(sum_formulae):
            if adduct not in mz_list[sum_formula]:#adduct may not be present if it would make an impossible formula, is there a better way to handle this?
                # print '{} adduct not found for {}'.format(adduct, mz_list[sum_formula])
                continue
            if time.time() - t_el > 10.:
                t_el = time.time()
                print '{:3.2f} done in {:3.0f} seconds'.format(float(ii)/len(sum_formulae),time.time()-t0)
            # Allocate dicts if required
            if not sum_formula in measure_value_score:
                    measure_value_score[sum_formula] = {}
            if not sum_formula in iso_correlation_score:
                    iso_correlation_score[sum_formula] = {}
            if not sum_formula in iso_ratio_score:
                    iso_ratio_score[sum_formula] = {}
            try:
                # 1. Generate ion images
                ion_datacube = IMS_dataset.get_ion_image(mz_list[sum_formula][adduct][0],
                                                         ppm)  # for each spectrum, sum the intensity of all peaks within tol of mz_list
                if do_preprocessing:
                    apply_image_processing(config,ion_datacube) #currently just supports hot-spot removal
                # 2. Spatial Chaos
                measure_value_score[sum_formula][adduct] = level_sets_measure.measure_of_chaos(
                    ion_datacube.xic_to_image(0), nlevels, interp=None, clean_im=False)[0]
                if measure_value_score[sum_formula][adduct] == 1:
                    measure_value_score[sum_formula][adduct] = 0
                # 3. Score correlation with monoiso
                if len(mz_list[sum_formula][adduct][1]) > 1:
                    iso_correlation_score[sum_formula][adduct] = isotope_image_correlation.isotope_image_correlation(
                        ion_datacube.xic, weights=mz_list[sum_formula][adduct][1][1:])
                else:  # only one isotope peak, so correlation doesn't make sense
                    iso_correlation_score[sum_formula][adduct] = 1
                # 4. Score isotope ratio
                iso_ratio_score[sum_formula][adduct] = isotope_pattern_match.isotope_pattern_match(ion_datacube.xic,
                                                                                                   mz_list[sum_formula][
                                                                                                       adduct][1])
            except KeyError as e:
                print str(e)
                print "bad key in: \"{}\" \"{}\" ".format(sum_formula, adduct)
        output_results(config, measure_value_score, iso_correlation_score, iso_ratio_score, sum_formulae, [adduct], mz_list)
    return measure_value_score, iso_correlation_score, iso_ratio_score


def check_pass(pass_thresh, pass_val):
    tf = []
    for v, t in zip(pass_val, pass_thresh):
        tf.append(v > t)
    if all(tf):
        return True
    else:
        return False


def score_results(config, measure_value_score, iso_correlation_score, iso_ratio_score):
    measure_tol = config['results_thresholds']['measure_tol']
    iso_corr_tol = config['results_thresholds']['iso_corr_tol']
    iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']
    sum_formulae, adducts, mz_list = generate_isotope_patterns(config)
    pass_formula = []
    for sum_formula in sum_formulae:
        for adduct in adducts:
            if check_pass((measure_tol, iso_corr_tol, iso_ratio_tol), (
            measure_value_score[sum_formula][adduct], iso_correlation_score[sum_formula][adduct],
            iso_ratio_score[sum_formula][adduct])):
                pass_formula.append('{} {}'.format(sum_formula, adduct))
    return pass_formula


def output_results(config, measure_value_score, iso_correlation_score, iso_ratio_score, sum_formulae, adducts, mz_list, fname=''):
    import os
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    measure_tol = config['results_thresholds']['measure_tol']
    iso_corr_tol = config['results_thresholds']['iso_corr_tol']
    iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']
    # sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    # Save the processing results
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    filename_out = generate_output_filename(config,adducts,fname=fname)
    with open(filename_out, 'w') as f_out:
        f_out.write('sf,adduct,mz,moc,spat,spec,pass\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                if adduct not in mz_list[sum_formula]:
                    continue
                p_vals = (
                    measure_value_score[sum_formula][adduct],
                    iso_correlation_score[sum_formula][adduct],
                    iso_ratio_score[sum_formula][adduct])
                moc_pass = check_pass((measure_tol, iso_corr_tol, iso_ratio_tol), p_vals)
                str_out = '{},{},{},{},{},{},{}\n'.format(
                    sum_formula,
                    adduct,
                    mz_list[sum_formula][adduct][0][0],
                    measure_value_score[sum_formula][adduct],
                    iso_correlation_score[sum_formula][adduct],
                    iso_ratio_score[sum_formula][adduct],
                    moc_pass)
                str_out.replace('[',"\"")
                str_out.replace(']',"\"")
                f_out.write(str_out)

def generate_output_filename(config,adducts,fname=''):
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    if fname == '':
        for adduct in adducts:
            fname='{}_{}'.format(fname,adduct)
    filename_out = '{}{}{}_{}_full_results.txt'.format(output_dir, os.sep,
                                                    os.path.splitext(os.path.basename(filename_in))[0],fname)
    return filename_out

def output_results_exactMass(config, ppm_value_score, sum_formulae, adducts, mz_list, fname=''):
    import os
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    # sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    # Save the processing results
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    if fname == '':
        for adduct in adducts:
            fname='{}_{}'.format(fname,adduct)

    filename_out = '{}{}{}_{}_exactMass_full_results.txt'.format(output_dir, os.sep,
                                                    os.path.splitext(os.path.basename(filename_in))[0],fname)
    with open(filename_out, 'w') as f_out:
        f_out.write('sf,adduct,mz,ppm\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                if adduct not in mz_list[sum_formula]:
                    continue
                str_out = '{},{},{},{}\n'.format(
                    sum_formula,
                    adduct,
                    mz_list[sum_formula][adduct][0][0],
                    ppm_value_score[sum_formula][adduct])
                str_out.replace('[',"\"")
                str_out.replace(']',"\"")
                f_out.write(str_out)

def output_results_frequencyFilter(config, ppm_value_score, sum_formulae, adducts, mz_list, fname=''):
    import os
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    # sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    # Save the processing results
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    if fname == '':
        for adduct in adducts:
            fname='{}_{}'.format(fname,adduct)

    filename_out = '{}{}{}_{}_frequencyFilter_full_results.txt'.format(output_dir, os.sep,
                                                    os.path.splitext(os.path.basename(filename_in))[0],fname)
    with open(filename_out, 'w') as f_out:
        f_out.write('sf,adduct,mz,fraction\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                if adduct not in mz_list[sum_formula]:
                    continue
                str_out = '{},{},{},{}\n'.format(
                    sum_formula,
                    adduct,
                    mz_list[sum_formula][adduct][0][0],
                    ppm_value_score[sum_formula][adduct])
                str_out.replace('[',"\"")
                str_out.replace(']',"\"")
                f_out.write(str_out)

def output_pass_results(config, measure_value_score, iso_correlation_score, iso_ratio_score, sum_formulae, adducts, mz_list, fname=''):
    import os
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    measure_tol = config['results_thresholds']['measure_tol']
    iso_corr_tol = config['results_thresholds']['iso_corr_tol']
    iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']
    # sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    # Save the processing results
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    if fname == '':
        for adduct in adducts:
            fname='{}_{}'.format(fname,adduct)
    filename_out = '{}{}{}_{}_pass_results.txt'.format(output_dir, os.sep,
                                                    os.path.splitext(os.path.basename(filename_in))[0],fname)
    with open(filename_out, 'w') as f_out:
        f_out.write('ID,sf,adduct,mz,moc,spec,spat\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                if adduct not in mz_list[sum_formula]:
                    continue
                if check_pass((measure_tol, iso_corr_tol, iso_ratio_tol), (
                        measure_value_score[sum_formula][adduct], iso_correlation_score[sum_formula][adduct],
                        iso_ratio_score[sum_formula][adduct])):
                    f_out.write('{},{},{},{},{},{},{}\n'.format(
                        sum_formulae[sum_formula]['db_id'],
                        sum_formula, adduct,
                        mz_list[sum_formula][adduct][0][0],
                        measure_value_score[sum_formula][adduct],
                        iso_correlation_score[sum_formula][adduct],
                        iso_ratio_score[sum_formula][adduct]))


def load_data(config):
    # Parse dataset
    from pyIMS.inMemoryIMS import inMemoryIMS
    IMS_dataset = inMemoryIMS(config['file_inputs']['data_file'])
    return IMS_dataset


def takeClosest(myList, myNumber):
        import bisect
        """
        Assumes myList is sorted. Returns closest value to myNumber.
        If two numbers are equally close, return the smallest number.
        """
        pos = bisect.bisect_left(myList, myNumber)
        if pos == 0:
            return (myList[0],pos)
        if pos == len(myList):
            return (myList[-1],pos)
        before = abs(myList[pos - 1]-myNumber)
        after = abs(myList[pos]-myNumber)
        if after <  before:
           return (after,pos)
        else:
           return (before,pos-1)


def run_exact_mass_search(config, mzs,counts, sum_formulae, adducts, mz_list):
    ### Runs the main pipeline
    # Get sum formula and predicted m/z peaks for molecules in database
    ppm_value_score = {}
    for sum_formula in sum_formulae:
        ppm_value_score[sum_formula]={}
    for adduct in adducts:
        for ii,sum_formula in enumerate(sorted(sum_formulae.keys())):
            if adduct not in mz_list[sum_formula]:#adduct may not be present if it would make an impossible formula, is there a better way to handle this?
                continue
            target_mz = mz_list[sum_formula][adduct][0][0]
            mz_nearest,pos = takeClosest(mzs, target_mz)
            ppm_value_score[sum_formula][adduct] = 1e6*mz_nearest/target_mz
        output_results_exactMass(config, ppm_value_score, sum_formulae, [adduct], mz_list)
    return ppm_value_score


def run_frequency_mass_search(config, IMS_dataset, sum_formulae, adducts, mz_list):
    ### Runs the main pipeline
    # Get sum formula and predicted m/z peaks for molecules in database
    freq_value_score = {}
    for sum_formula in sum_formulae:
        freq_value_score[sum_formula]={}
    for adduct in adducts:
        for ii,sum_formula in enumerate(sorted(sum_formulae.keys())):
            if adduct not in mz_list[sum_formula]:#adduct may not be present if it would make an impossible formula, is there a better way to handle this?
                continue
            target_mz = mz_list[sum_formula][adduct][0][0]
            ion_image = IMS_dataset.get_ion_image(np.asarray([target_mz,]),np.asarray([config['image_generation']['ppm'],]))
            freq_value_score[sum_formula][adduct] = np.sum(np.asarray(ion_image.xic)>0)/float(len(ion_image.xic[0]))
        output_results_exactMass(config, freq_value_score, sum_formulae, [adduct], mz_list)
    return freq_value_score


def fdr_selection(mz_list,pl_adducts, n_im):
    # produces a random subset of the adducts loaded in mz_list to actually calculate with
    pl_adducts = set(pl_adducts)
    for sf in mz_list:
        adduct_list = set(mz_list[sf].keys()) - pl_adducts # can be different for each molecule (e.g if adduct loss would be imposisble)
        rep=False
        if len(adduct_list)<n_im:
            rep=True
        keep_adducts = set(np.random.choice(list(adduct_list),n_im,replace=rep))|pl_adducts
        for a in adduct_list- keep_adducts:
            del mz_list[sf][a]
    return mz_list


def run_pipeline(JSON_config_file):
    config = get_variables(JSON_config_file)
    sum_formulae, adducts, mz_list = generate_isotope_patterns(config)
    if 'fdr' in config:
        mz_list = fdr_selection(mz_list,[str(a["adduct"]) for a in config['fdr']["pl_adducts"]], config['fdr']['n_im'])
    IMS_dataset = load_data(config)
    measure_value_score, iso_correlation_score, iso_ratio_score = run_search(config, IMS_dataset, sum_formulae, adducts,mz_list)
    # pass_list = score_results(config,measure_value_score, iso_correlation_score, iso_ratio_score)
    output_results(config, measure_value_score, iso_correlation_score, iso_ratio_score, sum_formulae, adducts, mz_list,fname='spatial_all_adducts')


def exact_mass(JSON_config_file):
    config = get_variables(JSON_config_file)
    sum_formulae, adducts, mz_list = generate_isotope_patterns(config)
    IMS_dataset = load_data(config)
    spec_axis,mean_spec =IMS_dataset.generate_summary_spectrum(summary_type='mean',ppm=config['image_generation']['ppm']/2)
    from pyMS.centroid_detection import gradient
    import numpy as np
    mzs,counts,idx_list = gradient(np.asarray(spec_axis),np.asarray(mean_spec),weighted_bins=2)
    ppm_value_score = run_exact_mass_search(config,  mzs,counts, sum_formulae, adducts, mz_list)
    output_results_exactMass(config, ppm_value_score, sum_formulae, adducts, mz_list,fname='exactMass_all_adducts')


def frequency_filter(JSON_config_file):
    config = get_variables(JSON_config_file)
    sum_formulae, adducts, mz_list = generate_isotope_patterns(config)
    IMS_dataset = load_data(config)
    #spec_axis,mean_spec =IMS_dataset.generate_summary_spectrum(summary_type='hist',ppm=config['image_generation']['ppm']/2)
    #from pyMS.centroid_detection import gradient
    #import numpy as np
    #mzs,counts,idx_list = gradient(np.asarray(spec_axis),np.asarray(mean_spec),weighted_bins=2)
    ppm_value_score = run_frequency_mass_search(config,  IMS_dataset, sum_formulae, adducts, mz_list)
    output_results_frequencyFilter(config, ppm_value_score, sum_formulae, adducts, mz_list,fname='frequencyFilter_all_adducts')
