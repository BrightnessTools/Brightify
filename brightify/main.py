#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
import pandas as pd
import numpy as np
from pathlib import Path
from io import StringIO
import pickle as pkl
from tqdm import tqdm
import matplotlib.pyplot as plt
import mcpl

from brightify import fast_computation as fc


FILE_ENDING = ".bfy" # File ending (including dot) used for brightify model saving

#%% Define particle types with their corresponding PDG codes
PDG_CODES = {
    'proton': 2212,  # Proton PDG code
    'neutron': 2112,  # Neutron PDG code
    'photon': 22,  # Photon PDG code
    'electron': 11  # Electron PDG code
}

#%% List of column names for different types of data (position, velocity, etc.)
# These will be used to organize and access particle data
PDG_CODE_NAME = 'Kf'
ENERGY_NAME = 'e'
WEIGHT_NAME = 'Wt'
POSITION_NAME_X, POSITION_NAME_Y, POSITION_NAME_Z = 'x', 'y', 'z'
POSITION_NAMES_XYZ = [POSITION_NAME_X, POSITION_NAME_Y, POSITION_NAME_Z] # Position components (x, y, z)
VELOCITY_NAME_X, VELOCITY_NAME_Y, VELOCITY_NAME_Z = 'Vx', 'Vy', 'Vz'
VELOCITY_NAMES_XYZ = [VELOCITY_NAME_X, VELOCITY_NAME_Y, VELOCITY_NAME_Z] # Velocity components (Vx, Vy, Vz)

# POSITION_NAMES = ['x', 'y', 'z']  # Position components (x, y, z)
# VELOCITY_NAMES = ['Vx', 'Vy', 'Vz']  # Velocity components (Vx, Vy, Vz)
# DIRECTION_NAMES = ['dirX', 'dirY', 'dirZ']  # Direction components (cosines in 3D space)
# IMAGE_DIRECTION_NAMES = ['imX', 'imY', 'imZ']  # Image plane direction components
# WEIGHTED_POSITION_NAMES = ['Xw', 'Yw', 'Zw']  # Weighted position components
# WEIGHTED_VELOCITY_NAMES = ['Vxw', 'Vyw', 'Vzw']  # Weighted velocity components
# MEAN_POSITION_NAMES = ['Xm', 'Ym', 'Zm']  # Mean position components
# MEAN_VELOCITY_NAMES = ['Vxm', 'Vym', 'Vzm']  # Mean velocity components

def read_data(inputFile, input_is_ascii=False):
    """
    Reads data from mcpl / mcpl.gz file and returns it as pandas dataframe.
    Optionally from ascii file.
    """
    inputFile=Path(inputFile).resolve()
    data = None
    if input_is_ascii:
        with open(str(inputFile), 'r') as f:
                data = f.read().replace('D', 'E')
        data = pd.read_csv(
            StringIO(data), header=None, delim_whitespace=True, 
            names=(PDG_CODE_NAME, POSITION_NAME_X, POSITION_NAME_Y, POSITION_NAME_Z,
                VELOCITY_NAME_X, VELOCITY_NAME_Y, VELOCITY_NAME_Z, ENERGY_NAME, WEIGHT_NAME,
                'time', 'c1', 'c2', 'c3', 'name', 'nocas', 'nobch', 'no')
            )
    elif inputFile.name.endswith(".mcpl") or inputFile.name.endswith(".mcpl.gz"):
        # Read the MCPL file (format used for Brightify input files)
        myfile=mcpl.MCPLFile(str(inputFile),blocklength=1e13) # the blocklength is so high, it will almost always load the entire file
        # Loop through each particle block in the MCPL file and store data
        # in a list of pandas DataFrames, and then combine all into one dataframe.
        data_frames = []
        for p in myfile.particle_blocks: # since we have a loop over the blocks, shouldn't each block's data be appended to self.data, instead of overwriting it each time?
            df = pd.DataFrame({
                PDG_CODE_NAME: p.pdgcode,  # Particle PDG code
                POSITION_NAME_X: p.x,  # x-position
                POSITION_NAME_Y: p.y,  # y-position
                POSITION_NAME_Z: p.z,  # z-position
                VELOCITY_NAME_X: p.ux,  # Velocity in x-direction
                VELOCITY_NAME_Y: p.uy,  # Velocity in y-direction
                VELOCITY_NAME_Z: p.uz,  # Velocity in z-direction
                ENERGY_NAME: p.ekin,  # Kinetic energy
                WEIGHT_NAME: p.weight  # Weight of the particle
            })
            data_frames.append(df)
        data = pd.concat(data_frames, ignore_index=True)
    else: raise ValueError("input could not be interpreted") 
    print("Data has been read from "+str(inputFile))
    return data  


#%%
class BrightifyModel:
    # This is the super-class.
    # The classes that inherit from it (via "class Subclass(Superclass)") I call sub-classes, e.g. here "Flat" and "Spherical".
    # The user is expected to start by calling "object = Subclass()" e.g. "flat_object = Flat()",
    # i.e. calling the constructor, and then work on that object by calling its variables/methods as e.g. "flat_object.variable1" or "flat_object.method1()".
    # 
    # The super-class is supposed to define variables and methods that any sub-class object would need in the same way.
    # They can be called just as "subclassobject.superclassvariable" or similarly for methods.
    # To enforce that a method MUST be implemented in all subclasses, the @abstractmethod decorator can be used.
    # For example, a variable that could be shared is "inputFile" or "data" since regardless of the Model type, that is needed.
    # A shared method would be for example "save_data" that just saves the data object as a dump file.
    # 
    # If a variable or method shall have special content/behavior for a sub-class object,
    # it should be defined in that sub class again, overwriting the functionality of that what is provided by the super-class.
    # It MUST be defined in the sub class if the superclass method has the @abstractmethod decorator.
    #
    # Regarding the constructors, we can defined a super-class constructor that does things that any sub-class constructor should do in the same way.
    # Defining sub-class constructors overwrites that super-class constructor, but we can call the super-class constructor as "super().__init__().
    # In general, any superclass thing can be called as "super().thing".
    # Note that the super class constructor should only contain parameters that all the sub-class constructors have in common.
    
    
    #  ----- 1. Object creation logic -----

    # Classical constructor
    def __init__(self, 
                 inputFile=None,
                 input_is_ascii=False,
                 outputFile=None,
                 primary_protons=None,
                 pCurrent=None,   
                 dir_size=None
                 ): 
        """
        Initializes the BrightifyModel class.
        
        Parameters:
        - inputFile: Path to the input file containing particle data (optional). Default accepted format is .mcpl(.gz).
        - input_is_ascii: Whether you are using ascii input separated by whitespace instead of .mcpl(.gz) files. 
        - outputFile: Optional path to file in which previous BrightifyModel was saved (pickle dump of attributes).
          Note: Please use brightify.load_model(outputFile) instead
        - primary_protons: Number of primary protons (used for brightness calculation). Can be found in standard output file of MC code
        - pCurrent: Proton current (used for brightness calculation).
        - dir_size: Size of the solid-angle window (in steradians) around the direction used for brightness calculation
        """

        # --- 0. catch (deprecated) case where user gives a dump file ---
        
        if outputFile is not None:
            raise ValueError("outputFile cannot be used here. Please use the class method load_model(outputFile) instead")
        self.outputFile = None

        # --- 1. Assign simple variables ---

        # inputFile
        if inputFile is not None:
            inputFile = Path(inputFile).resolve()
            if not inputFile.exists(): raise ValueError("inputFile path does not exist")
        self.inputFile = inputFile # Even though data is put into self.data, the user might want to check later what his inputFile has been.
        
        # input_is_ascii
        try: input_is_ascii = bool(input_is_ascii)
        except (ValueError, TypeError): raise ValueError("input_is_ascii must be bool-convertible")
        self.input_is_ascii = input_is_ascii
        
        # primary_protons
        try: primary_protons = int(primary_protons)
        except (ValueError, TypeError): raise ValueError("primary_protons must be int-convertible")
        self.primary_protons = primary_protons
        
        # pCurrent
        try: pCurrent = float(pCurrent)
        except (ValueError, TypeError): raise ValueError("pCurrent must be float-convertible")
        self.pCurrent = pCurrent
        
        # dir_size
        try: dir_size = float(dir_size)
        except (ValueError, TypeError): raise ValueError("dir_size must be float-convertible")
        self.dir_size = dir_size

        # --- 2. Data loading ---

        # this needs to be a pandas DataFrame
        self.data = read_data(inputFile=self.inputFile, input_is_ascii=self.input_is_ascii)
             
        # --- 3. Filter initialization ---
        
        # Reset filter to include all data initially    
        self.reset_filter()

        # --- 4. Initialize calculation method variable (holds most recent applied calculation method string e.g. "adaptive") ---
        self.calc_method = None

    # Alternative: construct object of model from saved dump
    @classmethod
    def load_model(cls, dumpFile, enforce_type=False):
        """
        Reads previously saved model and returns it.
        By default, the returned model is of the type that is saved in the file.
        - enforce_type: if True, then the type is the class from which you called this method.
        If you call brightify.load_model, then the enforced type would be the base BrightifyModel.
        """
        # This is a sort of alternative constructor.
        # It is a class method, meaning it can be called as Class.method (e.g. Flat.load_model())
        dumpFile = Path(dumpFile).resolve()
        if not dumpFile.suffix == FILE_ENDING:
            print("Provided file does not have the expected file extension '"+FILE_ENDING+"'. Looking for intended file...")
            dumpFile = dumpFile.with_suffix(FILE_ENDING)
        # this reads a previous dump of a model and returns the dictionary (does not change self object, and doesn't need it at all because it is static method)
        if not dumpFile.exists():
            raise FileNotFoundError("File not found at "+str(dumpFile))
        with open(dumpFile, "rb") as f:
            (modelType, attrDict) = pkl.load(f)
        modelType = cls if enforce_type else modelType
        model = modelType.__new__(modelType)
        model.__dict__.update(attrDict)
        print("Model of type "+modelType.__name__+" was loaded from "+str(dumpFile))
        print("The following attributes were filled with values: "+str(list(attrDict.keys())).strip("[]"))
        return model

    # ----- 2. Saving and loading logic -----

    SAVE_KEYS = ['inputFile',
                 'input_is_ascii',
                 'primary_protons',
                 'pCurrent',
                 'dir_size',
                 'calc_method']
        # This class variable is read upon calling the saving method on an object whose type is this class.
        # It determines which attributes of the object are saved when the saving method is called on that object.
        # DEFINE THE IDENTICALLY NAMED VARIABLE in the child classes; their lists will be added to this one,
        # automatically by the saving method logic.
    
    def get_save_keys(self):
        """
        Returns the names of attributes that are by default supposed to be included when saving.
        This relies on the lists of attribute names defined in the class variables named SAVE_KEYS.
        """
        lists_of_keys = []
        for cls in type(self).__mro__: # MRO = Method resolution order
            # this loop goes through all parent classes and looks for the SAVE_KEYS variable
            if hasattr(cls, "SAVE_KEYS"):
                # put keys at the beginning of the list of lists
                lists_of_keys.insert(0, cls.SAVE_KEYS)
        keys = [k for lok in lists_of_keys for k in lok]
        return keys

    def get_save_attributes(self):
        """
        Returns a dictionary of the attributes (names and data content) that 
        are included by default when saving.
        """
        _MISSING_ATTR = object() # a unique thing that no attribute will ever be
        keys = self.get_save_keys()
        attrDict = {key:val for key in keys if (val := getattr(self, key ,_MISSING_ATTR)) is not _MISSING_ATTR}
        # only existing attributes are included!
        return attrDict

    def save_model(self, dumpFile, overwrite=False, include_filter=False, include_data=False):
        """
        Saves model to dumpFile (pickle dump).
        File will get the default ending if it doesn't yet have it.
        The dumped object is a tuple: (type of model, attributes dictionary).
        - overwrite = True will overwrite existing file without asking, else saving is canceled.
        - include_filter = True will include the filter (length = number of particles in input data).
        - include_data = True will include the input data (self.data, which is a pandas dataFrame).
            WARNING: File might get large.
        """
        dumpFile = Path(dumpFile).resolve()
        # Prepare file path
        if not dumpFile.suffix == FILE_ENDING:
            dumpFile = dumpFile.with_suffix(dumpFile.suffix + FILE_ENDING)
        if dumpFile.exists() and not overwrite:
                raise FileExistsError("File was not written because it already exists at "+str(dumpFile))
        # get attributes and type and save tuple (type, dictionary) to the location
        modelType = type(self)
        attrDict = self.get_save_attributes()
        # include/exclude filter:
        if not include_filter and "filter" in attrDict.keys():
            attrDict.pop("filter") 
        elif include_filter and hasattr(self, "filter"):
            attrDict["filter"] = getattr(self, "filter")
        # include/exclude data:
        if not include_data and "data" in attrDict.keys():
            attrDict.pop("data") 
        elif include_data and hasattr(self, "data"):
            attrDict["data"] = getattr(self, "data")
        # save
        save_tuple = (modelType, attrDict)
        with open(dumpFile, "wb") as f:
                pkl.dump(save_tuple, f) 
        print("A dump file of the model of type "+self.__class__.__name__+" has been saved to "+str(dumpFile)
              +" and the following attributes were included: "+str(list(attrDict.keys())).strip("[]"))

    # ----- 3. Calculations logic -----

    # def coordinate_transform(self):
    #     - coordinate_system: This program expects the z coordinate direction to be the normal direction of the brightness maps. \
    #                 You may want to transform your data beforehand yourself, or input here three Euler angles (alpha, beta, gamma) (in degrees).
    #                 alpha determines how much the x axis needs to be 
    
    @property
    def data_filtered(self):
        """
        Property method to return filtered data based on the current filter
        Returns the data that matches the current filter conditions.
        """
        return self.data[self.filter]
    
    @property
    def brightness(self):
        """
        Property method to calculate brightness based on provided parameters.
        Calculates brightness based on particle current, primary protons,
        position size, and direction size.
        """
        return (self.window_weights * self.pCurrent / self.primary_protons / self.pos_size / self.dir_size)
    
    @property
    def r_err(self):
        """
        method to calculate relative error based on the standard deviation
        (sigma).
        """
        
        return self.sigma / np.sqrt(self.window_weights)

    def compute_sigma(self, slice_filt_i, Wt_sq_sum, Am, window_weight):
        """
        Compute the standard deviation (sigma) for a given filtered data slice.
        """
        if window_weight <= 1:
            return 0  # Prevent division by zero
    
        # Compute sigma using precomputed values
        return np.sqrt(
            np.abs(fc.sum_f(slice_filt_i[WEIGHT_NAME]**2) / Wt_sq_sum - Am * window_weight) /
            (window_weight - 1)
        )
    
    def reset_filter(self):
        """
        Method to reset the filter to include all data.
        Resets the filter so that no data is excluded by default.
        """
        self.filter = self.data.index >= 0 # this can just be a ones-array of same length as data        
    
    def apply_filter(self, particle=None, energy=None, x=None, y=None): # for large data would it make sense to filter while creating the object (i.e. in the constructor) so not all data is read into memory?
        """
        Applies filters to the particle data based on the given criteria:
        - particle: Type of particle (e.g., 'proton', 'neutron')
        - energy: Energy range to filter particles by
        - x: x-position range to filter particles by
        - y: y-position range to filter particles by
        
        If you want to reset the filter, use the method "reset_filter()".
        """

        # If filter hasn't been initialized yet:
        if not hasattr(self, "filter"):
            self.reset_filter()

        # If a particle type is specified, filter by the particle's PDG code
        if particle:
            self.particle = particle
            self.filter &= (self.data[PDG_CODE_NAME] == PDG_CODES[particle])
        
        # If an energy range is specified, filter particles within that range
        if energy:
            self.energy = (min(energy), max(energy))
            self.filter &= (self.data[ENERGY_NAME] >= self.energy[0])
            self.filter &= (self.data[ENERGY_NAME] <= self.energy[1])
            
        # If an x-position range is specified, filter particles within that range
        if x:
            self.xf = (min(x), max(x))
            self.filter &= (self.data[POSITION_NAME_X] >= self.xf[0])
            self.filter &= (self.data[POSITION_NAME_X] <= self.xf[1])
            
        # If a y-position range is specified, filter particles within that range
        if y:
            self.yf = (min(y), max(y))
            self.filter &= (self.data[POSITION_NAME_Y] >= self.yf[0])
            self.filter &= (self.data[POSITION_NAME_Y] <= self.yf[1])

        print("filter has been applied")


    def calculate_properties(self, data_filt, shape1D, shape2D, 
                             slice_filter_func, points, forced_direction=None):
        """
        Common calculations for total weights, mean weights, directions, and errors.
        
        Parameters
        ----------
        forced_direction : np.ndarray or None
            If provided, use this direction for all slices instead of computing mean.
        """
        Wt_sq_sum = fc.sum_f(self.data[WEIGHT_NAME]**2)
        Am = np.mean(self.data[WEIGHT_NAME]**2)
    
        self.total_neutrons = np.full(shape1D, np.nan)
        self.total_weights = np.full(shape1D, np.nan)
        self.window_weights = np.full(shape1D, np.nan)
        self.total_mean = np.full(shape1D, np.nan)
        self.relative_error = np.full(shape1D, np.nan)
        self.sigma = np.full(shape1D, np.nan)
        self.mean_directions = np.full(shape2D, np.nan)
        self.mean_directions_spher = np.full(shape2D, np.nan)
        self.mean_directions_window = np.full(shape2D, np.nan)
    
        for i, vec in tqdm(enumerate(points), desc='Calculating', total=shape1D, ascii=' #'):
            slice_i = slice_filter_func(data_filt, vec)
            if slice_i.empty:
                continue
    
            self.total_neutrons[i] = len(slice_i)
            self.total_weights[i] = fc.sum_f(slice_i.Wt)
            self.total_mean[i] = np.mean(slice_i.Wt)
    
            # Use forced direction if provided, else compute mean direction
            if forced_direction is not None:
                direction_vec = forced_direction
            else:
                direction_vec = fc.direction_f(fc.mean_f(slice_i[VELOCITY_NAMES_XYZ], slice_i.Wt))
    
            self.mean_directions[i] = direction_vec
            self.mean_directions_spher[i] = fc.cart2spher(*direction_vec)
    
            # Apply angular filter
            slice_filt_i = slice_i[fc.cos_f(slice_i[VELOCITY_NAMES_XYZ], direction_vec) >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Sigma and relative error
            self.sigma[i] = self.compute_sigma(slice_filt_i=slice_filt_i, Wt_sq_sum=Wt_sq_sum, Am=Am, window_weight=self.window_weights[i])
            self.relative_error[i] = self.r_err[i]
    
            # Weighted direction
            self.mean_directions_window[i] = direction_vec * self.window_weights[i]
          
    def calculate_properties_adaptive(self, data_filt, shape1D, shape2D, 
                                  slice_filter_func, points, bins_cart=35, threshold_frac=0.98):
        """
        Adaptive calculation of total weights, mean weights, directions, 
        and errors using Cartesian histogram to find maxima.
        
        Parameters
        ----------
        bins_cart : int
            Number of bins along each Cartesian axis (X, Y, Z) for histogram.
        threshold_frac : float
            Fraction of the global histogram maximum to define maxima voxels.
        """
        # Calculate mean squared weight (Am) and sum of squared weights (Wt_sq_sum)
        Wt_sq_sum = fc.sum_f(self.data[WEIGHT_NAME]**2)  # Global sum of squared weights
        Am = np.mean(self.data[WEIGHT_NAME]**2)  # Mean squared weight
    
        # Initialize arrays
        self.total_neutrons = np.full(shape1D, np.nan)
        self.total_weights = np.full(shape1D, np.nan)
        self.window_weights = np.full(shape1D, np.nan)
        self.total_mean = np.full(shape1D, np.nan)
        self.relative_error = np.full(shape1D, np.nan)
        self.sigma = np.full(shape1D, np.nan)
        self.mean_directions = np.full(shape2D, np.nan)
        self.mean_directions_spher = np.full(shape2D, np.nan)
        self.mean_directions_window = np.full(shape2D, np.nan)
        self.adaptive_dir = np.full(shape2D, np.nan)
        self.adaptive_dir_spher = np.full(shape2D, np.nan)
        self.hist = [None] * shape1D
        self.hist_edges = [None] * shape1D
    
        # Loop through the points
        for i, vec in tqdm(enumerate(points), desc='Calculating', total=shape1D, ascii=' #'):
            # Slice data
            slice_i = slice_filter_func(data_filt, vec)  
            if len(slice_i) == 0:
                continue
    
            # Store basic stats
            self.total_neutrons[i] = len(slice_i)
            self.total_weights[i] = fc.sum_f(slice_i.Wt)
            self.total_mean[i] = np.mean(slice_i.Wt)
    
            # Mean direction of slice
            self.mean_directions[i] = fc.direction_f(fc.mean_f(slice_i[VELOCITY_NAMES_XYZ], slice_i.Wt))         
            self.mean_directions_spher[i] = fc.cart2spher(*self.mean_directions[i])
    
            # ------------------ Cartesian histogram ------------------
            hist, edges = np.histogramdd(slice_i[VELOCITY_NAMES_XYZ].to_numpy(),
                                         bins=[bins_cart]*3,
                                         range=[[-1,1],[-1,1],[-1,1]],
                                         weights=slice_i[WEIGHT_NAME])
            self.hist[i] = hist
            self.hist_edges[i] = edges
    
            # Threshold maxima
            threshold = threshold_frac * hist.max()
            maxima_idx = np.argwhere(hist >= threshold)
    
            # Convert maxima voxel indices to Cartesian centers
            maxima_coords = []
            for idx in maxima_idx:
                x = 0.5 * (edges[0][idx[0]] + edges[0][idx[0]+1])
                y = 0.5 * (edges[1][idx[1]] + edges[1][idx[1]+1])
                z = 0.5 * (edges[2][idx[2]] + edges[2][idx[2]+1])
                maxima_coords.append([x, y, z])
            maxima_coords = np.array(maxima_coords)
    
            # Compute adaptive direction as mean of maxima vectors
            if len(maxima_coords) > 0:
                maxima_dirs_unit = maxima_coords / np.linalg.norm(maxima_coords, axis=1)[:, np.newaxis]
                adaptive_dir_vec = np.mean(maxima_dirs_unit, axis=0)
                adaptive_dir_vec /= np.linalg.norm(adaptive_dir_vec)
            else:
                adaptive_dir_vec = np.array([np.nan, np.nan, np.nan])
    
            self.adaptive_dir[i] = adaptive_dir_vec
            self.adaptive_dir_spher[i] = fc.cart2spher(*adaptive_dir_vec)
    
            # Filter directions within angular window
            cos_filter = fc.cos_f(slice_i[VELOCITY_NAMES_XYZ], self.adaptive_dir[i])
            slice_filt_i = slice_i[cos_filter >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Sigma and weighted mean direction
            self.sigma[i] = self.compute_sigma(slice_filt_i=slice_filt_i, Wt_sq_sum=Wt_sq_sum, Am=Am, window_weight=self.window_weights[i])
            self.relative_error[i] = self.r_err[i]
            self.mean_directions_window[i] = self.mean_directions[i] * self.window_weights[i]


class Spherical(BrightifyModel):

    #--------- Saving logic helper variable ----------
    SAVE_KEYS = ['pos_size',
                 'particle', 'energy', 'radius', 'spiral_points', 'spiral',
                 'pos_window', 'dir_window', 'total_neutrons', 'total_weights',
                 'window_weights', 'total_mean', 'relative_error', 'mean_directions',
                 'mean_directions_window'] 

    #------------------ Constructor -------------------------------#
    def __init__(self, 
                 inputFile=None, input_is_ascii=False,
                 outputFile=None,
                 primary_protons=None, pCurrent=None,
                 dir_size=None, 
                 pos_size=None,
                 radius=None,
                 spiral_points=None):
        """
        Initialize a Spherical model, child of BrightifyModel.
        Base parameters:
        - inputFile: Path to the input file containing particle data (optional). Default accepted format is .mcpl(.gz).
        - input_is_ascii: Whether you are using ascii input separated by whitespace instead of .mcpl(.gz) files. 
        - outputFile: Optional path to file in which previous BrightifyModel was saved (pickle dump of attributes).
            Note: Please use brightify.load_model(outputFile) instead
        - primary_protons: Number of primary protons (used for brightness calculation). Can be found in standard output file of MC code
        - pCurrent: Proton current (used for brightness calculation).
        - dir_size: Size of the solid-angle window (in steradians) around the direction used for brightness calculation            
        Additional Spherical class parameters:
        - pos_size: position window size.
        - radius: radius of the spherical shell of the dump source defined by user.
        - spiral_points: number of points for spiral scan.

        Returns
        -------
        None.

        """
        # --- 0. Initialize via the parent class (BrightifyModel) constructor with relevant arguments ---
        
        super().__init__(inputFile=inputFile, input_is_ascii=input_is_ascii,
                         outputFile=outputFile,
                         primary_protons=primary_protons, pCurrent=pCurrent,
                         dir_size=dir_size)

        # --- 1. Assigning simple variables --- 

        # pos_size
        try: pos_size = float(pos_size)
        except (ValueError, TypeError): raise ValueError("pos_size must be float-convertible")
        self.pos_size = pos_size

        # radius
        try: radius = float(radius)
        except (ValueError, TypeError): raise ValueError("radius must be float-convertible")
        self.radius = radius

        # spiral_points
        try: spiral_points = int(spiral_points)
        except (ValueError, TypeError): raise ValueError("spiral_points must be int-convertible")
        self.spiral_points = spiral_points

        # --- 2. Spiral generation ---
        self.spiral = Spherical.spiral_scan(spiral_points=spiral_points)
              
    def spherical_filter_func(self, data_filt, vec):
        """
        Filtering function specific to the Spherical subclass. 
        It filters `slice_i` based on the cosine of the angle between `vec` and `data_filt`.
        """
        return data_filt[fc.cos_f(data_filt[POSITION_NAMES_XYZ], vec) > self.pos_window]

    def calculate(self):
        data_filt = self.data_filtered
        self.pos_window = 1 - self.pos_size / (2 * np.pi * self.radius**2)  # Position window
        self.dir_window = 1 - self.dir_size / (2 * np.pi)  # Direction window (d(omega) = dS/r^2)
    
        return self.calculate_properties(data_filt=data_filt, shape1D=self.spiral.shape[0],
                                          shape2D=self.spiral.shape,
                                          slice_filter_func=self.spherical_filter_func,
                                          points=self.spiral)            

    @staticmethod
    def spiral_scan(spiral_points):
        """
        This function generates a spiral scan based on the golden spiral pattern.
        
        Parameters
        ----------
        spiral_points : int
            The number of points for the spiral scan.
        
        Returns
        -------
        np.ndarray
            An array of coordinates representing the points of the spiral in spherical coordinates (theta, phi).
        """
        assert type(int(spiral_points)) == int
        
        # Generate the spiral angles (theta, phi) based on the golden spiral pattern
        indices = np.arange(0, spiral_points, dtype=np.float64) + 0.5
        spiral_theta = np.arccos(1 - 2 * indices / spiral_points)
        spiral_phi = np.pi * (1 + 5**0.5) * indices
        
        # Convert the spherical coordinates (theta, phi) to Cartesian coordinates
        return fc.spher2cart(1, spiral_theta, spiral_phi)
    
    def plot(self):
        """
        Plots a 3D brightness map of the spherical shell.
        """
        
        ax = plt.figure().add_subplot(111, projection='3d')
        sc = ax.scatter(
            self.spiral[:,2], self.spiral[:,1], self.spiral[:,0],
            c=self.brightness, cmap=plt.cm.magma, vmin=0 #c=self.window_weights
            )
        ax.quiver(0, 0, 0, 1, 0, 0, color='r', length=1.0, normalize=True)
        ax.text(1.5, 0, 0, "Z", color='r')
        ax.quiver(0, 0, 0, 0, 1, 0, color='g', length=1.0, normalize=True)
        ax.text(0, 1.5, 0, "Y", color='g')
        ax.quiver(0, 0, 0, 0, 0, 1, color='b', length=1.0, normalize=True)
        ax.text(0, 0, 1.5, "X", color='b')
        ax.quiver(
            self.spiral[:,2], self.spiral[:,1], self.spiral[:,0], 
            self.mean_directions_window[:,2],
            self.mean_directions_window[:,1],
            self.mean_directions_window[:,0],
            length=0.0001
            )
        plt.colorbar(sc, ax=ax)

    def summary(self):
        """
        Returns a summary of the results of the spherical scan.
        
        Returns
        -------
        dict
            A dictionary containing spiral points, mean directions, and brightness.
        """
        summaryDict = {
            'spiral': [*self.spiral],
            'mean_directions': [*self.mean_directions],
            'brightness': [*self.brightness]
            }
        return summaryDict
        

class Flat(BrightifyModel):

    #--------- Saving logic helper variable ----------
    SAVE_KEYS = ['pos_size_x', 'pos_size_y', 'pos_size',
                 'particle', 'energy', 'dir_window', 'x_range', 'y_range',
                 'x_mesh', 'y_mesh', 'z_mesh', 'total_neutrons', 'total_weights',
                 'window_weights', 'total_mean', 'relative_error', 'sigma',
                 # for the case method == "mean" (and "normal"):
                 'mean_directions', 'mean_directions_window', 'mean_directions_spher'
                 # for the case method == "adaptive":
                 'mean_directions','mean_directions_spher','adaptive_dir', 'adaptive_dir_spher', 'hist','hist_edges',
                 ]
    
    #------------------ Constructor -------------------------------#
    def __init__(self, inputFile=None, input_is_ascii=False,
                 outputFile=None,
                 primary_protons=None, pCurrent=None, 
                 dir_size=None,
                 pos_size_x=None, pos_size_y=None):
        """
        Initialize a Flat model, child of BrightifyModel.
        Base parameters:
        - inputFile: Path to the input file containing particle data (optional). Default accepted format is .mcpl(.gz).
        - input_is_ascii: Whether you are using ascii input separated by whitespace instead of .mcpl(.gz) files. 
        - outputFile: Optional path to file in which previous BrightifyModel was saved (pickle dump of attributes).
            Note: Please use brightify.load_model(outputFile) instead
        - primary_protons: Number of primary protons (used for brightness calculation). Can be found in standard output file of MC code
        - pCurrent: Proton current (used for brightness calculation).
        - dir_size: Size of the solid-angle window (in steradians) around the direction used for brightness calculation            
        Additional Flat class parameters:
        - pos_size_x: Horizontal size of the position window (used for brightness calculation).
        - pos_size_y: Vertical size of the position window.
        """
        
        # --- 0. Initialize via the parent class (BrightifyModel) constructor with relevant arguments ---
        
        super().__init__(inputFile=inputFile, input_is_ascii=input_is_ascii,
                         outputFile=outputFile,
                         primary_protons=primary_protons, pCurrent=pCurrent, 
                         dir_size=dir_size)

        # --- 1. Initialize position window variables ---
        
        # pos_size_x
        try: pos_size_x = float(pos_size_x)
        except (ValueError, TypeError): raise ValueError("pos_size_x must be float-convertible")
        self.pos_size_x = pos_size_x

        # pos_size_y
        try: pos_size_y = float(pos_size_y)
        except (ValueError, TypeError): raise ValueError("pos_size_y must be float-convertible")
        self.pos_size_y = pos_size_y

        self.pos_size = pos_size_x * pos_size_y

    #------------------ Helper: prepare meshgrid ------------------#
    def prepare_vectors(self):
        data_filt = self.data_filtered
        self.x_min, self.x_max = data_filt.x.min(), data_filt.x.max()
        self.y_min, self.y_max = data_filt.y.min(), data_filt.y.max()
        self.x_range = np.arange(self.x_min, self.x_max + 1, self.pos_size_x/4)
        self.y_range = np.arange(self.y_min, self.y_max + 1, self.pos_size_y/4)
        self.x_mesh, self.y_mesh = np.meshgrid(self.x_range, self.y_range)
        self.z_mesh = np.full_like(self.x_mesh, 1.)
        self.vectors = np.stack((self.x_mesh.flatten(), self.y_mesh.flatten(),
                                 self.z_mesh.flatten()), axis=-1)

    #------------------ Filter function ---------------------------#
    def flat_filter_func(self, data_filt, vec):
        return data_filt[(data_filt[POSITION_NAME_X] >= vec[0] - self.pos_size_x/2) & 
                         (data_filt[POSITION_NAME_X] <  vec[0] + self.pos_size_x/2) & 
                         (data_filt[POSITION_NAME_Y] >= vec[1] - self.pos_size_y/2) & 
                         (data_filt[POSITION_NAME_Y] <  vec[1] + self.pos_size_y/2)]

    #------------------ Unified calculation -----------------------#
    def calculate(self, method="adaptive", bins_cart=35, threshold_frac=0.98):
        """
        Calculate brightness and relative error on the mesh.
        
        method: "mean", "adaptive", or "normal"
        """
        data_filt = self.data_filtered
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        self.prepare_vectors()
        result = None
        if method == "mean":
            result = self.calculate_properties(
                data_filt=data_filt,
                shape1D=self.vectors.shape[0],
                shape2D=self.vectors.shape,
                slice_filter_func=self.flat_filter_func,
                points=self.vectors
            )
        elif method == "adaptive":
            result = self.calculate_properties_adaptive(
                data_filt=data_filt,
                shape1D=self.vectors.shape[0],
                shape2D=self.vectors.shape,
                slice_filter_func=self.flat_filter_func,
                points=self.vectors,
                bins_cart=bins_cart,
                threshold_frac=threshold_frac
            )
        elif method == "normal":
            # Use forced surface-normal direction
            result = self.calculate_properties(
                data_filt=data_filt,
                shape1D=self.vectors.shape[0],
                shape2D=self.vectors.shape,
                slice_filter_func=self.flat_filter_func,
                points=self.vectors,
                forced_direction=np.array([0,0,1])
            )
        else:
            raise ValueError(f"Unknown calculation method: {method}")
        
        self.calc_method = method
        return result

    #------------------ Unified plotting ------------------#
    def plot_brightness_map(self, method=None, show_arrows=True):
        """
        Plot brightness map with optional auto-detection of method
        used in calculate(), and arrows for direction vectors.
    
        method: "mean", "adaptive", "normal", or None (auto-detect)
        show_arrows: whether to draw direction arrows
        """
    
        # Ensure bounds exist
        if not hasattr(self, "x_min") or not hasattr(self, "x_max"):
            self.x_min, self.x_max = self.x_range.min(), self.x_range.max()
            self.y_min, self.y_max = self.y_range.min(), self.y_range.max()
    
        # Auto-detect method if not provided
        if method is None: method = self.calc_method
    
        # Select data
        if method == "mean":
            brightness_map = self.brightness.reshape(len(self.y_range),
                                                     len(self.x_range))
            directions = self.mean_directions
            directions_spher = self.mean_directions_spher
    
        elif method == "adaptive":
            brightness_map = self.window_weights.reshape(len(self.y_range),
                                                         len(self.x_range))
            directions = self.adaptive_dir
            directions_spher = self.adaptive_dir_spher
    
        elif method == "normal":
            brightness_map = self.window_weights.reshape(len(self.y_range),
                                                         len(self.x_range))
            directions = np.tile(np.array([[0, 0, 1]]), (self.x_mesh.size, 1))
            directions_spher = np.stack((np.zeros(len(directions)),
                                         np.zeros(len(directions))), axis=1)
        else:
            raise ValueError(f"Unknown calculation method: {method}")
    
        circle_centers = np.stack((self.x_mesh.flatten(), self.y_mesh.flatten()), axis=-1)
        circle_radii = np.sin(directions_spher[:, 1])
    
        fig, ax = plt.subplots(figsize=(16, 12))
        im = ax.imshow(brightness_map, cmap=plt.cm.viridis, origin='lower',
                       extent=[self.x_min, self.x_max, self.y_min, self.y_max])
    
        if show_arrows:
            for center, radius, vector in zip(circle_centers, circle_radii, directions):
                ax.arrow(center[0], center[1], radius * vector[0], radius * vector[1],
                         head_width=0.1, head_length=0.2, fc='blue', ec='blue')
    
        ax.set_xlabel('x [cm]', fontsize=18)
        ax.set_ylabel('y [cm]', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        cbar = plt.colorbar(im, label='brightness')
        cbar.ax.tick_params(labelsize=18)
        cbar.ax.set_ylabel('brightness [n/s/cm$^2$/sr]', fontsize=18)
        plt.show()

    def plot_error_map(self, method=None, show_arrows=True):
        """
        Plot relative error map with optional auto-detection of method
        used in calculate(), and arrows for direction vectors.
        """
    
        # Ensure bounds exist
        if not hasattr(self, "x_min") or not hasattr(self, "x_max"):
            self.x_min, self.x_max = self.x_range.min(), self.x_range.max()
            self.y_min, self.y_max = self.y_range.min(), self.y_range.max()
    
        # Auto-detect method if not provided
        if method is None: method = self.calc_method
    
        # Select data
        if method == "mean":
            error_map = self.relative_error.reshape(len(self.y_range),
                                                    len(self.x_range))
            directions = self.mean_directions
            directions_spher = self.mean_directions_spher
    
        elif method == "adaptive":
            error_map = np.full_like(self.window_weights, np.nan).reshape(
                len(self.y_range), len(self.x_range)
            )
            directions = self.adaptive_dir
            directions_spher = self.adaptive_dir_spher
    
        elif method == "normal":
            error_map = self.relative_error.reshape(len(self.y_range),
                                                    len(self.x_range))
            directions = np.tile(np.array([[0, 0, 1]]), (self.x_mesh.size, 1))
            directions_spher = np.stack((np.zeros(len(directions)),
                                         np.zeros(len(directions))), axis=1)
        else:
            raise ValueError(f"Unknown calculation method: {method}")
    
        circle_centers = np.stack((self.x_mesh.flatten(), self.y_mesh.flatten()), axis=-1)
        circle_radii = np.sin(directions_spher[:, 1])
    
        fig, ax = plt.subplots(figsize=(16, 12))
        im = ax.imshow(error_map, cmap=plt.cm.viridis, origin='lower',
                       extent=[self.x_min, self.x_max, self.y_min, self.y_max])
    
        if show_arrows:
            for center, radius, vector in zip(circle_centers, circle_radii, directions):
                ax.arrow(center[0], center[1], radius * vector[0], radius * vector[1],
                         head_width=0.1, head_length=0.2, fc='blue', ec='blue')
    
        ax.set_xlabel('x [cm]', fontsize=18)
        ax.set_ylabel('y [cm]', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        cbar = plt.colorbar(im, label='relative error')
        cbar.ax.tick_params(labelsize=18)
        cbar.ax.set_ylabel('relative error', fontsize=18)
        plt.show()
     
    def surface_crossing(self, v_x, v_y, v_z, theta_D):
        """
        Method to calculate neutron surface crossing for PHITS comparison
        w.r.t a given direction
        
        """
        data_filt = self.data_filtered
        dir_D = np.array([v_x, v_y, v_z])
        slice_filt_i = data_filt[fc.cosRR_f(data_filt[VELOCITY_NAMES_XYZ], dir_D) >= np.cos(theta_D)]
        self.window_weights = fc.sum_f(slice_filt_i.Wt)
        
        # Compute statistical uncertainty
        Wt_sq_sum = fc.sum_f(self.data[WEIGHT_NAME]**2)  # Global sum of squared weights
        Am = np.mean(self.data[WEIGHT_NAME]**2)  # Mean squared weight
        self.sigma = self.compute_sigma(slice_filt_i=slice_filt_i, Wt_sq_sum=Wt_sq_sum, Am=Am, window_weight=self.window_weights)
    
            # Calculate relative error
        self.relative_error = self.r_err

        return self.brightness, self.relative_error

    def surface_crossing_all_surface(self):
        """   
        Method to calculate surface crossing for the entire surface
        w.r.t mean direction
        
        """
        data_filt = self.data_filtered
        
        # d(omega) = dS/r^2 = 3*6*1e-4/10^2 = 0.18
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        
        # Define the range of x and y values
        self.x_min, self.x_max = data_filt.x.min(), data_filt.x.max()
        self.y_min, self.y_max = data_filt.y.min(), data_filt.y.max()
        
        Wt_sq_sum = fc.sum_f(self.data[WEIGHT_NAME]**2)  # Global sum of squared weights
        Am = np.mean(self.data[WEIGHT_NAME]**2)  # Mean squared weight
        
        slice_i = data_filt
        self.total_neutrons = len(slice_i)
        self.total_weights = fc.sum_f(slice_i.Wt)
        self.total_mean = np.mean(slice_i.Wt)
        
        self.mean_directions = fc.direction_f(fc.mean_f(slice_i[VELOCITY_NAMES_XYZ], slice_i.Wt))
        
        slice_filt_i = slice_i[
            fc.cos_f(slice_i[VELOCITY_NAMES_XYZ], self.mean_directions) >= self.dir_window
            ]
        self.window_weights = fc.sum_f(slice_filt_i.Wt)
        
        self.sigma = self.compute_sigma(slice_filt_i=slice_i, Wt_sq_sum=Wt_sq_sum, Am=Am, window_weight=self.window_weights)
    
            # Calculate relative error
        self.relative_error = self.r_err

  
        return self.brightness         


    