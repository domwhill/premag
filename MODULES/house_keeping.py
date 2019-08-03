'''
	Module containing paths to run data and saving data
'''
class directory_paths(object):
    def __init__(self):
        self.src_dir = '/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/'
        self.data_dir_2D = '/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/PREMAG/2D_RUNS/'
        self.data_dir_1D = '/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/PREMAG/2D_RUNS/'
        self.log_dir = '' # path to normalisation file
        self.norm_dir = '/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/'
#-----------------------------------------------------------------------
class load_names(object):
    '''
    '''
    def __init__(self):
        path_in = directory_paths().data_dir_2D
        self.single_B = '%sr5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E' % (path_in)
        self.single_noB = '%sr5_v40_Z_FEOS_MODNE_5y_matchedf0_in_0T_E' % (path_in)
        
        self.speckle_B = '%sr5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E' % (path_in)
        self.speckle_noB = '%sr5_v40_Z_FEOS_MODNE_5y_matchedf0_in_1T_E' % (path_in)
        
        self.save_path = './pics/r_5v38/'
        self.time_glb = '13'
        return
