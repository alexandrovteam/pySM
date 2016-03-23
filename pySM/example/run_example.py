# Before running this:
# - add pyMSpec,pyIMS and pySpatialMetabolomics to the sys.path
# - grab the example data from metabolights
#   * http://www.ebi.ac.uk/metabolights/MTBLS317
#   * http://www.ebi.ac.uk/metabolights/MTBLS313
# - edit the  .json config file 
#!!! These are mandatory edits for every dataset
### "file_inputs":{  
###      "data_file":"PATH/TO/sm_examples/simulated_data/decoy_dataset_chemnoise_centroids_IMS.hdf5",
###      "database_load_folder":"/ADD_PATH/",

###      "results_folder":"PATH/TO/sm_examples/simulated_data",
###      "database_file":"PATH/TO/sm_examples/database/hmdb_organic_database.csv"
#       * "data_file" - path to imzml file
#       * "database_load_folder" - can be '' or a folder. if a folder then isotope patterns are saved and loaded from this location othewise they are generated each time
#       * "results_folder" - root folder for results, results will be saved in "results_folder"/name
#       * "database_file" - csv file with first three columns: <id, formula, name>
###    "name":"hmdb_decoy_simulated_dataset",
#       * "name" - short name for the dataset
# - other fields can be edited as required (e.g. change pl_adducts for negative mode)
# Run Pipeline
from  pySM import spatial_metabolomics, fdr_measures
from pySM.tools import results_tools
json_filenames = ["RB_a1s1.json",
                  "RB_a2s1.json",
                  "RB_a2s2.json"
                 ]
for json_filename in json_filenames:
    spatial_metabolomics.run_pipeline(json_filename)

# View results
for json_filename in json_filenames:
    target_adducts,decoy_adducts = fdr_measures.get_adducts(json_filename)
    config = spatial_metabolomics.get_variables(json_filename)
    results_fname = spatial_metabolomics.generate_output_filename(config,[],'spatial_all_adducts')
    fdr = fdr_measures.decoy_adducts(results_fname,target_adducts,decoy_adducts)
    fdr_target=0.1
    n_reps=20
    pass_sf = fdr.decoy_adducts_get_pass_list(fdr_target,n_reps,col='msm')
    (measure_value_score, iso_correlation_score,iso_ratio_score, moc_pass) = results_tools.load_results(results_fname)
    print 'Formula annotated at {} from {}:'.format(fdr_target, config['name'])
    print "{},{},{},{},{},{}".format('sum formula','adduct','p_chaos','p_spatial','p_spectral','msm')
    for adduct in pass_sf:
        for sf in pass_sf[adduct]:
            moc = measure_value_score[sf][adduct]
            spat = iso_correlation_score[sf][adduct]
            spec = iso_ratio_score[sf][adduct]
            msm = moc*spat*spec
            print "{},{},{},{},{},{}".format(sf,adduct,moc,spat,spec,msm)