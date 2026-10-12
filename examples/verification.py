#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# This script can be run from console like "python verification.py".
# Or interactively (with #%% syntax) where I don't know what the working directory is.
# To ensure both works, a few helpers are needed, in particular for file paths.
helpers=True
if helpers: # helpers
    try:
        import IPython
        ipython = IPython.get_ipython()
        if ipython is not None:
            ipython.run_line_magic("load_ext", "autoreload")
            ipython.run_line_magic("autoreload", "2")
    except (ImportError, AttributeError):
        pass

    from pathlib import Path

    def from_top_path(top_path):
        """
        resolves a given path that is relative to the top folder "Brightify"
        """
        this_file_path = Path(__file__)
        examples_folder_path = this_file_path.parent # provided that this file is in the folder "examples"
        Brightify_folder_path = examples_folder_path.parent # provided "examples" is in "Brightify"
        path = Path(Brightify_folder_path / top_path.lstrip('./')).resolve()
        print("valid path: "+str(path) if path.exists() else "unvalid path: "+str(path))
        return str(path) # to show that the package can take simple strings, convert again to string.

#%%
# Choose data (path defined relative to top directory "Brightify"):

path = './examples/uniform/output.mcpl.gz'
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/Brightify Fig 6/Fig6a/output.mcpl.gz' # The number of primary particles is 4.7e7 and the resolution is 1 cm2 and 0.006 sr, based on the paper 
# path = './examples/mcpl_different_normal/dump_final.mcpl.gz'

if helpers: path = from_top_path(top_path=path)

#%%
# Now we can bring in the brightify package:
import brightify as bfy
# it contains these model classes:
Flat = bfy.Flat
Spherical = bfy.Spherical
# brightify also has a method used for reading data (e.g. mcpl) into a dataframe object:
read_data = bfy.read_data
# and a method for restoring a previously saved model (saving via model.save_model)
load_model = bfy.load_model

#%%
# construct a Flat object, reading in data on the way:
flatModel = Flat(
    inputFile=path, 
    primary_protons = 4.7e7, # number
    pCurrent = 1.0,
    pos_size_x = 1.0, # in cm
    pos_size_y = 1.0, # in cm
    dir_size = 6e-3 # in sr (steradians)
    )
#%%
# apply a filter to the data:
flatModel.apply_filter(
    particle = 'neutron',
    energy = (0.0, 2.0),
    )
# (Note that you could also apply the particle filter first, and then a second one for the energy.)
# (That works because the apply_filter method is cumulative (a "monotonic filter"). To reset, use "reset_filter()")
#%%
# make a calculation and save the result
flatModel.calculate(method="mean")
flatModel.save_model(dumpFile="test", overwrite=True, include_filter=False, include_data=False)
# note that files are saved with the ".bfy" extension (but you can leave it out in the argument)
#%%
# look at the results:
flatModel.plot_brightness_map()
flatModel.plot_error_map()
#%%
#now load again from file
flatModel_from_file = load_model(dumpFile="test") 
# Note we left out the file extension here too, causing a file with the extension added to be searched.
#%%
# and try to plot something:
flatModel_from_file.plot_brightness_map()
# this does plot something!
#%%
# But, you might want to have a different filter this time.
# In our case you need to separately load and insert the data first (because we used include_data=False when saving)
flatModel_from_file.data = read_data(inputFile=path)
#%%
# Then apply the different filter
flatModel_from_file.apply_filter(
    particle = 'neutron',
    energy = (0.0, 1.0), # this is different than above
    )
# (Note that our disk save didn't include the filter (via include_filter=False).
# (If the old filter would exist (in the "filter" attribute), then we 
# could should run "reset_filter()" first.
#%%
flatModel_from_file.calculate(method="normal") # different too
flatModel_from_file.plot_brightness_map() # looks different
# %%
