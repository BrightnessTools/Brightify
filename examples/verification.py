#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
from pathlib import Path

this_file_path = Path(__file__)
examples_folder_path = this_file_path.parent
project_folder_path = examples_folder_path.parent

#%%
path = './examples/uniform/output.mcpl.gz'
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/mcpl_different_normal/dump_final.mcpl.gz'

path = str(project_folder_path) + path[1:] 
print(path)
from os import path as path_module
print(path_module.exists(path))
#%%
from brightify import Flat
#%%
# import brightify
# print(brightify.main.ENERGY_NAME)
#%%
flatModel = Flat(
    inputFile=path, 
    primary_protons = 4.7e7, # number
    pCurrent = 1.0,
    pos_size_x = 1.0, # in cm
    pos_size_y = 1.0, # in cm
    dir_size = 6e-3 # in sr (steradians)
    )
#%%
flatModel.apply_filter(
    particle = 'neutron',
    energy = (0.0, 2.0),
    )

flatModel.calculate()

flatModel.plot_brightness_map()

flatModel.plot_error_map()
