import h5py
import numpy as np
def get_coords(hdf,index_list):
	# retrieve coordinates from an _IMS.hdf5 document - each spectrum has a coordinate attached
        coords = np.zeros((len(index_list),3))
        for k in index_list:
            coords[k,:] = hdf['/spectral_data/'+str(k)+'/coordinates/']
        return coords

def dump_coord(data_id=1,hdf_filename,f_out):
	# write coords to file as formatted for sm server
    import h5py
    import numpy as np
    hdf = h5py.File(hdf_filename,'r')   #Readonly, file must exist
    index_list = map(int,hdf['/spectral_data'].keys())
    def get_coords(hdf,index_list):
        coords = np.zeros((len(index_list),3))
        for k in index_list:
            coords[k,:] = hdf['/spectral_data/'+str(k)+'/coordinates/']
        return coords
    coords=get_coords(hdf,index_list)
	#subtract minimum value from each coordinate set to pull to a useful frame of reference
    min_0 = min([c[0] for c in coords])
    min_1 = min([c[1] for c in coords])
    min_2 = min([c[2] for c in coords])
    # save to file
    with open(f_out,'w') as f_out:
        for ii,coord in enumerate(coords):
            f_out.write('{},{},{},{},{},\n'.format(data_id,ii,(coord[0]-min_0),(coord[1]-min_1),coord[2]-min_2))
        
f_in = '/Volumes/alexandr/shared/ASMS_Server_Test_Data/Ctrl3s2_SpheroidsCtrl_DHBSub_centroids_IMS.hdf5' #using a temporary hdf5 based format
f_out = '{}/ctrl3s2_coords.txt'.format('/Users/palmer/Documents/Projects/2015/SM_development/coordinates')
data_id = 1
dump_coord(data_id,f_in,f_out)
