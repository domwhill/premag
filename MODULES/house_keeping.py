'''
	Module containing paths to run data and saving data
'''
import numpy as np
import os
import json, pdb


class directory_paths(object):

    def __init__(self):
        dirname, filename = os.path.split(os.path.abspath(__file__))
        data_path = os.path.abspath(dirname + '../../../TF/')

        self.src_dir = os.path.abspath(dirname + '/../') + '/'
        self.data_dir_2D = os.path.abspath(data_path + '/2D_RUNS/') + '/'
        self.data_dir_1D = os.path.abspath(data_path + '/1D_RUNS/') + '/'
        self.save_dir = '%sOUTPUT/' % self.src_dir    # path to normalisation file
        self.norm_dir = self.src_dir

    #-------------------------------------------------------------------
    def get_fprefix(self, dim, scale_length, lambda_p='5', bz_in='50', pert_amp='1p'):

        if pert_amp == '1p':
            amp_str = ''
        else:
            amp_str = '_' + pert_amp

        n2s = lambda num: str(int(num))

        filename = self.src_dir + 'MODULES/paths.json'
        with open(filename, "r") as f:
            data = json.loads(f.read())
        if dim == '1D' or str(dim) == '1':
            return data['1D'][n2s(bz_in) + amp_str]
        else:
            return data['2D'][n2s(scale_length)][n2s(lambda_p)][n2s(bz_in) + amp_str]

    #-------------------------------------------------------------------

    def get_path(self, scale_length, Bz, lambda_p=5, dim='2D', pert_amp='1p'):
        """
            scale_length = str
            Bz = integer (can also be string eg. '-1')
                        'bz_100p' = 100% perturbation
                        '-1' =  Bz switched off
                        '0' = 0 T (no applied field)
                        50 = 50 T applied field
                        100 = 100T applied field
                        400 - 400T applied field.
                
            dim = '1D' or '2D'
            lambda_p = integer: Perturbation wavelength in lambda_mfp,ref
             5, 10 (not applicable for 1D runs)
            
        """
        print('=== ', dim, scale_length, lambda_p, Bz, pert_amp)
        path = 0
        fprefix = self.get_fprefix(dim, scale_length, lambda_p, Bz, pert_amp)
        if dim == '1D' or dim == '1' or dim == 1:
            if scale_length != 1:
                print('-- only scale length LT1 runs available in 1D currently --')

            dir_path = directory_paths().data_dir_1D
            path = dir_path + fprefix

        elif dim == '2D':
            dir_path = self.data_dir_2D
            scale_length = int(scale_length)
            if scale_length == 1:
                dir_path = dir_path + 'LT1/'

            elif scale_length == 2:
                dir_path = dir_path + 'LT2/'

            elif scale_length == 3:
                dir_path = dir_path + 'LT3/'
            else:
                print('scale length must be integer between 1 and 3')
                print('chosen = ', scale_length)
            if path == 0:
                print('no path available for this selection: ', lambda_p, Bz)
        path = dir_path + fprefix

        return path
