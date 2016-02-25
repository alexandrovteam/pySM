# Before running this:
# - add pyMS,pyIMS and pySpatialMetabolomics to the sys.path
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
import sys
sys.path.append('/path/to/python_codebase/pyMS')
sys.path.append('/path/to/python_codebase/pyIMS')
sys.path.append('/path/to/python_codebase/pySM')
# Run Pipeline
from pySpatialMetabolomics import spatial_metabolomics
json_filename = "./sm_paper_data/sim001_squares/sim001_standalone_config.json"
spatial_metabolomics.run_pipeline(json_filename)
# View results
from  pySpatialMetabolomics import spatial_metabolomics, fdr_measures
target_adducts = ["H","Na","K"]   
decoy_adducts = ["He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Ir", "Th", "Pt", "Pu", "Os", "Yb", "Lu", "Bi", "Pb", "Re", "Tl", "Tm", "U", "W", "Au", "Er", "Hf", "Hg", "Ta"]   
config = spatial_metabolomics.get_variables(json_filename)
results_fname = spatial_metabolomics.generate_output_filename(config,[],'spatial_all_adducts')
fdr = fdr_measures.decoy_adducts(results_fname,target_adducts,decoy_adducts)
fdr_target=0.1
n_reps=10
fdr.decoy_adducts_get_pass_list(fdr_target,n_reps,col='msm')
