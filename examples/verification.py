#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
try:
    import IPython
    ipython = IPython.get_ipython()
    if ipython is not None:
        ipython.run_line_magic("load_ext", "autoreload")
        ipython.run_line_magic("autoreload", "2")
except (ImportError, AttributeError):
    pass

#%%
from pathlib import Path

this_file_path = Path(__file__)
examples_folder_path = this_file_path.parent # provided that this file is in the folder "examples"
Brightify_folder_path = examples_folder_path.parent

#%%

# enter paths relative to the top folder "Brightify" that contains the examples folder

path = './examples/uniform/output.mcpl.gz'
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/mcpl_different_normal/dump_final.mcpl.gz'


path = Path(Brightify_folder_path / path.lstrip('./'))
print(path, "was found")
print("valid path: "+str(path) if path.exists() else "unvalid path: "+str(path))


#%%
import brightify as bfy

print(bfy.main.ENERGY_NAME)

#%%
flatModel = bfy.Flat(
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
#%%
flatModel.calculate()
#%%
flatModel.plot_brightness_map()
#%%
flatModel.plot_error_map()

# %%
