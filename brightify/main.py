#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
import pandas as pd
import numpy as np
from os import path
import pickle as pkl
from tqdm import tqdm
import matplotlib.pyplot as plt
import mcpl
from brightify import fast_computation as fc

#%% Define particle types with their corresponding PDG codes
PARTICLE = {
    'proton': 2212,  # Proton PDG code
    'neutron': 2112,  # Neutron PDG code
    'photon': 22,  # Photon PDG code
    'electron': 11  # Electron PDG code
}

#%% List of column names for different types of data (position, velocity, etc.)
# These will be used to organize and access particle data
pos = ['x', 'y', 'z']  # Position components (x, y, z)
vel = ['Vx', 'Vy', 'Vz']  # Velocity components (Vx, Vy, Vz)
direction = ['dirX', 'dirY', 'dirZ']  # Direction components (cosines in 3D space)
image_dir = ['imX', 'imY', 'imZ']  # Image plane direction components
pos_w = ['Xw', 'Yw', 'Zw']  # Weighted position components
vel_w = ['Vxw', 'Vyw', 'Vzw']  # Weighted velocity components
pos_m = ['Xm', 'Ym', 'Zm']  # Mean position components
vel_m = ['Vxm', 'Vym', 'Vzm']  # Mean velocity components

#%%
class BrightifyModel:
    def __init__(self, inputFile=None, outputFile=None, primary_protons=None,
                 pCurrent=None, pos_size=None, dir_size=None):
        """
        Initializes the BrightifyModel class.
        
        Parameters
        ----------
        inputFile: str
            Path to the input file containing particle data.
        outputFile: str
            Path to the output file where results are saved. This may be the path to a previously saved outputFile.
            In that case, it will be attempted to load values from it.
        primary_protons: int
            Number of primary protons; Can be found in standard output file of MC code.
        pCurrent: float
            Proton current.
        pos_size: float
            Spatial size (area) of the position window.
        dir_size: float
            Size of the direction window in steradians.

        Returns
        -------
        None.
        """
                
        # --- 1. Assign attributes, converting to correct types (unless None, which is allowed until after the reading of outputFile) ---
        self.inputFile = str(inputFile) if inputFile is not None else None
        self.outputFile = str(outputFile) if outputFile is not None else None
        self.primary_protons = int(primary_protons) if primary_protons is not None else None
        self.pCurrent = float(pCurrent) if pCurrent is not None else None
        self.pos_size = float(pos_size) if pos_size is not None else None
        self.dir_size = float(dir_size) if dir_size is not None else None

        # --- 2. Check correctness (unless None) and calculate additional attributes ---
        if self.inputFile is None and self.outputFile is None: # Ensure that either inputFile or outputFile exists
            raise ValueError("both inputFile and outputFile are None")
        if self.primary_protons is not None and self.primary_protons < 1:
            raise ValueError("primary_protons incorrect")
        if self.pCurrent is not None and self.pCurrent <= 0:
            raise ValueError("pCurrent incorrect")
        if self.pos_size is not None and self.pos_size <= 0:
            raise ValueError("pos_size incorrect")
        if self.dir_size is not None and self.dir_size <= 0:
            raise ValueError("dir_size incorrect")
        
        # If an input file is provided, load the particle data from the file
        if self.inputFile is not None:
            self.data = self.read_mcpl_data(self.inputFile)
        
        # If an output file is provided, load previously saved attributes unless already initialized with a different value (except None)
        if self.outputFile is not None:
            try: # Try to open and read the outputFile
                with open(self.outputFile, 'rb') as f:
                    dataDict = pkl.load(f)
                old_inputFile = dataDict.pop("inputFile", None) # this removes the key (if it exists). In the below code, that is benefitial.
                if old_inputFile != self.inputFile:
                    print("Attention: The outputFile says that the previously used inputFile path is different from the entered inputFile parameter.")
                    print("The previously used inputFile was located at "+str(old_inputFile))
                if self.inputFile is None and old_inputFile is not None:
                    print("The previous inputFile will be used now, because the entered inputFile parameter is None.")
                    self.data = self.read_mcpl_data(old_inputFile)
                    self.inputFile = old_inputFile
                # find out if currently new parameters/attributes are different than saved parameters/attributes:
                parameters_are_different = False
                for k,v in self.__dict__.items(): # loop over all currently new parameters
                    if v is None: continue # We do not mind parameters that are None being overwritten by previously used values, right?
                    if k in dataDict and dataDict[k] != v: # if outputFile has the parameter but with a different value:
                        parameters_are_different = True
                        break
                if not parameters_are_different:
                    for k,v in dataDict.items():
                        self.__dict__[k] = v
                    print("Parameters/attributes from outputFile were loaded.")
                else:
                    print("No parameters/attributes from outputFile were loaded because the entered parameters are different.")
            except: # In this case, self.outputFile is just the location where future results are supposed to be stored
                print("Results will be saved by default at "+self.outputFile)
        
        # Reset filter to include all data initially    
        self.reset_filter()

        # --- 3. Final check that we now have all that we need (should do that in the child classes too for their special attributes) ---
        for attr in ["primary_protons", "pCurrent", "pos_size", "dir_size", "data", "filter"]:
            if getattr(self, attr, None) is None:
                raise ValueError("The model is missing a critical attribute: "+attr)

    @staticmethod
    def read_mcpl_data(inputFile):
        # Read the MCPL file (format used for Brightify input files)
        myfile=mcpl.MCPLFile(inputFile,blocklength=1e13)
        # Loop through each particle block in the MCPL file and store data
        # in a pandas DataFrame
        datas = []
        for p in myfile.particle_blocks:
            datas.append(pd.DataFrame({
                'Kf': p.pdgcode,  # Particle PDG code
                'x': p.x,  # x-position
                'y': p.y,  # y-position
                'z': p.z,  # z-position
                'Vx': p.ux,  # Velocity in x-direction
                'Vy': p.uy,  # Velocity in y-direction
                'Vz': p.uz,  # Velocity in z-direction
                'e': p.ekin,  # Kinetic energy
                'Wt': p.weight  # Weight of the particle
            }))
            data = pd.concat(datas, ignore_index=True)
        return data

    # @staticmethod
    # def read_ascii_data(inputFile):
    #     # if the ascii version is used as input
    #     with open(inputFile, 'r') as f:
    #         data = f.read().replace('D', 'E')
    #     data = pd.read_csv(
    #         StringIO(data), header=None, delim_whitespace=True, 
    #         names=('Kf', 'x', 'y', 'z', 'Vx', 'Vy', 'Vz', 'e', 'Wt',
    #                'time', 'c1', 'c2', 'c3', 'name', 'nocas', 'nobch', 'no'
    #                )
    #         )
    #     return data  

    def apply_rotation(self, z_vector=[0,0,1], x_vector=[1,0,0], matrix=None, replace=False): # includes filter reset
        """
        Rotates the coordinate system of the model's data such that:
        - z_vector ([x,y,z]): Its direction is the Z-axis in the new data. This is required for 2D models like Flat which need Z to be the normal direction.
        - x_vector ([x,y,z]): Its projection onto the new X-Y-plane is the new X-axis.
        - matrix: Pass here any 3x3 matrix to be applied to the data instead of the two vector parameters. You should yourself ensure it is a rotation matrix.
        - replace: If True, the rotation that's saved in the model will be the method input's rotation, instead of the total rotation that results in combination with previously applied rotations.

        ATTENTION: This resets any filter that had been applied. After rotation, spatial data filters you apply need to be written with the new coordinates in mind! 
        """
        # this relies on self.rotation
        # reminder: usually, we calculate v'^T = R * v^T where v^T is a column vector.
        # here, the vectors are all rows of data. So instead to above, we do: v' = v * R^T,
        # because v' = (v'^T)^T = (R * v^T)^T = (v^T)^T * R^T = v * R^T.
        # that means instead of R * v^T, we write v * R^T to apply R.
        # In order to reverse it, we write v' * R since that is (v * R^T)*R = v * 1 = v.
        if matrix is not None:
            R = matrix
        else:
            # 1. Define the new Z-axis from the normal vector
            z_new = np.array(z_vector, dtype=float)
            z_new /= np.linalg.norm(z_new)
            # 2. Define the new X-axis using the reference vector
            ref = np.array(x_vector, dtype=float)
            ref /= np.linalg.norm(ref)
            # If the magnitude of the cross product is nearly zero, they are parallel
            if np.linalg.norm(np.cross(z_new, ref)) < 1e-6:
                raise ValueError("The reference vector cannot be parallel to the normal vector.")
            # Remove any component of ref that points along z_new (Gram-Schmidt orthogonalization)
            x_new = ref - np.dot(ref, z_new) * z_new
            x_new /= np.linalg.norm(x_new)
            # 3. Complete the right-handed coordinate system by finding the new Y-axis
            y_new = np.cross(z_new, x_new)
            # 4. Construct the rotation matrix (where each row is one of the new axes)
            R = np.vstack([x_new, y_new, z_new])
            # 5. Apply the TRANSPOSE of the rotation to position and velocity columns
        self.data.loc[:, pos] = self.data[pos].values @ R.T
        self.data.loc[:, vel] = self.data[vel].values @ R.T
        self.reset_filter()
        print("A coordinate rotation has been applied, and therefore the filter was reset.")
        # store/update rotation
        if replace or not hasattr(self, "rotation"):
            self.rotation = R
        else:
            self.rotation = R @ self.rotation

    def reset_rotation(self, undo=True): # includes filter reset
        """
            Reset rotation of coordinates. Resets filter too.
            - undo: Set True if you want to apply the inverse of previous rotation(s) done to the data. 
            Set to False if you just want to delete the attribute that stores the applied rotations, and do not want to change the data.

            ATTENTION: After rotation-undoing, spatial data filters you apply need to be written with the new coordinates in mind!
        """
        if not hasattr(self, "rotation"):
            print("No rotation found to reset.")
            return
        if undo:
            print("Undoing rotation...")
            self.data.loc[:, pos] = self.data[pos].values @ self.rotation
            self.data.loc[:, vel] = self.data[vel].values @ self.rotation
        self.reset_filter()
        del self.rotation
        print("Rotation was reset.")


    @property
    def data_filter(self):
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
            np.abs(fc.sum_f(slice_filt_i['Wt']**2) / Wt_sq_sum - Am * window_weight) /
            (window_weight - 1)
        )
    
    def reset_filter(self):
        """
        Method to reset the filter to include all data.
        Resets the filter so that no data is excluded by default.
        """
        self.filter = self.data.index >= 0
    
    
    def apply_filter(self, particle=None, energy=None, x=None, y=None):
        """
        Applies filters to the particle data based on the given criteria:
        - particle: Type of particle (e.g., 'proton', 'neutron')
        - energy: Energy range to filter particles by
        - x: x-position range to filter particles by
        - y: y-position range to filter particles by
        """
        # If a particle type is specified, filter by the particle's PDG code
        if particle:
            self.particle = particle
            self.filter &= (self.data['Kf'] == PARTICLE[particle])
        
        # If an energy range is specified, filter particles within that range
        if energy:
            self.energy = (min(energy), max(energy))
            self.filter &= (self.data['e'] >= self.energy[0])
            self.filter &= (self.data['e'] <= self.energy[1])
            
        # If an x-position range is specified, filter particles within that range
        if x:
            self.xf = (min(x), max(x))
            self.filter &= (self.data['x'] >= self.xf[0])
            self.filter &= (self.data['x'] <= self.xf[1])
            
        # If a y-position range is specified, filter particles within that range
        if y:
            self.yf = (min(y), max(y))
            self.filter &= (self.data['y'] >= self.yf[0])
            self.filter &= (self.data['y'] <= self.yf[1])
            
    
    def calculate_properties(self, data_filt, shape1D, shape2D, 
                             slice_filter_func, points, forced_direction=None):
        """
        Common calculations for total weights, mean weights, directions, and errors.
        
        Parameters
        ----------
        forced_direction : np.ndarray or None
            If provided, use this direction for all slices instead of computing mean.
        """
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)
        Am = np.mean(self.data['Wt']**2)
    
        self.total_neutrons = np.full(shape1D, np.nan)
        self.total_weights = np.full(shape1D, np.nan)
        self.window_weights = np.full(shape1D, np.nan)
        self.total_mean = np.full(shape1D, np.nan)
        self.relative_error = np.full(shape1D, np.nan)
        self.sigma = np.full(shape1D, np.nan)
        self.mean_directions = np.full(shape2D, np.nan)
        self.mean_directions_spher = np.full(shape2D, np.nan)
        self.mean_directions_window = np.full(shape2D, np.nan)
    
        for i, vec in tqdm(enumerate(points), desc='common', total=shape1D, ascii=' #'):
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
                direction_vec = fc.direction_f(fc.mean_f(slice_i[vel], slice_i.Wt))
    
            self.mean_directions[i] = direction_vec
            self.mean_directions_spher[i] = fc.cart2spher(*direction_vec)
    
            # Apply angular filter
            slice_filt_i = slice_i[fc.cos_f(slice_i[vel], direction_vec) >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Sigma and relative error
            self.sigma[i] = self.compute_sigma(slice_filt_i, Wt_sq_sum, Am, self.window_weights[i])
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
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)  # Global sum of squared weights
        Am = np.mean(self.data['Wt']**2)  # Mean squared weight
    
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
        for i, vec in tqdm(enumerate(points), desc='common', total=shape1D, ascii=' #'):
            # Slice data
            slice_i = slice_filter_func(data_filt, vec)  
            if len(slice_i) == 0:
                continue
    
            # Store basic stats
            self.total_neutrons[i] = len(slice_i)
            self.total_weights[i] = fc.sum_f(slice_i.Wt)
            self.total_mean[i] = np.mean(slice_i.Wt)
    
            # Mean direction of slice
            self.mean_directions[i] = fc.direction_f(fc.mean_f(slice_i[vel], slice_i.Wt))         
            self.mean_directions_spher[i] = fc.cart2spher(*self.mean_directions[i])
    
            # ------------------ Cartesian histogram ------------------
            hist, edges = np.histogramdd(slice_i[vel].to_numpy(),
                                         bins=[bins_cart]*3,
                                         range=[[-1,1],[-1,1],[-1,1]],
                                         weights=slice_i['Wt'])
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
            cos_filter = fc.cos_f(slice_i[vel], self.adaptive_dir[i])
            slice_filt_i = slice_i[cos_filter >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Sigma and weighted mean direction
            self.sigma[i] = self.compute_sigma(slice_filt_i, Wt_sq_sum, Am, self.window_weights[i])
            self.relative_error[i] = self.r_err[i]
            self.mean_directions_window[i] = self.mean_directions[i] * self.window_weights[i]



class Spherical(BrightifyModel):
    def __init__(self, inputFile=None, outputFile=None, primary_protons=None, 
                 pCurrent=None, pos_size=None, dir_size=None, radius=None,
                 spiral_points=None):
        """
        Initialize the Spherical model, which is a subclass of BrightifyModel.
        
        Parameters
        ----------
        inputFile: str
            Path to the input file containing particle data.
        outputFile: str
            Path to the output file where results are saved.
        primary_protons: int
            Number of primary protons; Can be found in standard output file of MC code.
        pCurrent: float
            Proton current.
        pos_size: float
            Spatial size (area) of the position window.
        dir_size: float
            Size of the direction window in steradians.
        radius : float
            radius of the spherical shell of the dump source defined by user.
        spiral_points : int
            number of points for spiral scan.

        Returns
        -------
        None.
        """

        # --- 1. Assign attributes (only those that super() doesn't already assign), converting to correct types (unless None) --- 
        self.radius = float(radius) if radius is not None else None
        self.spiral_points = int(spiral_points) if spiral_points is not None else None

        # --- 2. Check correctness (unless None) and calculate additional attributes ---
        if self.radius is not None and self.radius <= 0:
            raise ValueError("radius incorrect")
        if self.spiral_points is not None:
            if self.spiral_points >= 1:
                self.spiral = Spherical.spiral_scan(self.spiral_points)
            else:
                raise ValueError("spiral_points incorrect")
               
        # Initialize the parent class (BrightifyModel) with relevant arguments
        super().__init__(inputFile=inputFile, outputFile=outputFile,
                         primary_protons=primary_protons, pCurrent=pCurrent, 
                         pos_size=pos_size,
                         dir_size=dir_size)
        
        # --- 3. Final check that we now have all that we need (no need to check parent attributes) ---
        for attr in ["radius", "spiral_points", "spiral"]:
            if getattr(self, attr, None) is None:
                raise ValueError("The model is missing a critical attribute: "+attr)

        
    def spherical_filter_func(self, data_filt, vec):
        """
        Filtering function specific to the Spherical subclass. 
        It filters `slice_i` based on the cosine of the angle between `vec` and `data_filt`.
        """
        return data_filt[fc.cos_f(data_filt[pos], vec) > self.pos_window]
    

    def calculate(self):
        data_filt = self.data_filter
        self.pos_window = 1 - self.pos_size / (2 * np.pi * self.radius**2)  # Position window
        self.dir_window = 1 - self.dir_size / (2 * np.pi)  # Direction window (d(omega) = dS/r^2)
    
        return self.calculate_properties(data_filt, self.spiral.shape[0],
                                          self.spiral.shape,
                                          self.spherical_filter_func,
                                          self.spiral)            
       
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

        Returns
        -------
        fig, ax
            The Matplotlib figure and 3D axes for further customization.
        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
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
        fig.colorbar(sc, ax=ax)
        if hasattr(self, "rotation"):
            print("Attention: This plot is displaying data in new coordinates compared to the input data. \
                    The input data was transformed by applying the following rotation matrix R in the sense `R * (x,y,z)^T` (column vectors):\n"
                    +str(self.rotation))
        return fig, ax


    def save(self, outputFile=None):
        """
        Saves the model's data to a pickle file for future use.
        
        Parameters
        ----------
        outputFile : str
            The output file path where the data will be saved. If None is provided, the outputFile that was specified at model creation is used.
        
        Returns
        -------
        None.
        """
        # Specify the keys of the attributes to save
        keys = ('inputFile', 'primary_protons', 'pCurrent', 'pos_size', 'dir_size',
                'particle', 'energy', 'radius', 'spiral_points', 'spiral',
                'pos_window', 'dir_window', 'total_neutrons', 'total_weights',
                'window_weights', 'total_mean', 'relative_error', 'mean_directions',
                'mean_directions_window'
                )
        
        # Create a dictionary of data to save
        dataDict = {k: self.__dict__[k] for k in keys if k in self.__dict__}

        if outputFile is None:
            outputFile = self.outputFile
        
        # Save the data dictionary to a pickle file
        with open(outputFile, 'xb') as f:
            pkl.dump(dataDict, f)        
    
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
    def __init__(self, inputFile=None, outputFile=None, primary_protons=None, 
                 pCurrent=None, pos_size_x=None, pos_size_y=None, dir_size=None):
        """
        Initialize the Flat model, subclass of BrightifyModel.
        
        Parameters
        ----------
        inputFile: str
            Path to the input file containing particle data.
        outputFile: str
            Path to the output file where results are saved. This may be the path to a previously saved outputFile.
            In that case, it will be attempted to load values from it.
        primary_protons: int
            Number of primary protons; Can be found in standard output file of MC code.
        pCurrent: float
            Proton current.
        pos_size_x: float
            Horizontal size of the position window.
        pos_size_y: float
            Vertical size of the position window.
        dir_size: float
            Size of the direction window in steradians.

        Returns
        -------
        None.
        """
        # --- 1. Assign attributes (only those that super() doesn't already assign), converting to correct types (unless None) ---
        self.pos_size_x = float(pos_size_x) if pos_size_x is not None else None
        self.pos_size_y = float(pos_size_y) if pos_size_y is not None else None
        
        # --- 2. Check correctness (unless None) and calculate additional attributes ---
        if self.pos_size_x is not None and self.pos_size_x <= 0:
            raise ValueError("pos_size_x incorrect")
        if self.pos_size_y is not None and self.pos_size_y <= 0:
            raise ValueError("pos_size_y incorrect")
        
        # Initialize the parent class (BrightifyModel) with relevant arguments
        pos_size = pos_size_x * pos_size_y if (pos_size_x is not None) and (pos_size_y is not None) else None
        super().__init__(inputFile=inputFile, outputFile=outputFile,
                         primary_protons=primary_protons, pCurrent=pCurrent, 
                         pos_size=pos_size,
                         dir_size=dir_size)

        # --- 3. Final check that we now have all that we need (no need to check parent attributes) ---
        for attr in ["pos_size_x", "pos_size_y", "pos_size"]:
            if getattr(self, attr, None) is None:
                raise ValueError("The model is missing a critical attribute: "+attr)


    #------------------ Helper: prepare meshgrid ------------------#
    def prepare_vectors(self):
        data_filt = self.data_filter
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
        return data_filt[(data_filt['x'] >= vec[0] - self.pos_size_x/2) & 
                         (data_filt['x'] <  vec[0] + self.pos_size_x/2) & 
                         (data_filt['y'] >= vec[1] - self.pos_size_y/2) &
                         (data_filt['y'] <  vec[1] + self.pos_size_y/2)]

    #------------------ Unified calculation -----------------------#
    def calculate(self, method="adaptive", bins_cart=35, threshold_frac=0.98):
        """
        Calculate brightness and relative error on the mesh.

        method: "mean", "adaptive", or "normal"
        """
        data_filt = self.data_filter
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        self.prepare_vectors()
        self.calculation_method = method

        if method == "mean":
            return self.calculate_properties(
                data_filt,
                self.vectors.shape[0],
                self.vectors.shape,
                self.flat_filter_func,
                self.vectors
            )
        elif method == "adaptive":
            return self.calculate_properties_adaptive(
                data_filt,
                self.vectors.shape[0],
                self.vectors.shape,
                self.flat_filter_func,
                self.vectors,
                bins_cart,
                threshold_frac
            )
        elif method == "normal":
            # Use forced surface-normal direction
            return self.calculate_properties(
                data_filt,
                self.vectors.shape[0],
                self.vectors.shape,
                self.flat_filter_func,
                self.vectors,
                forced_direction=np.array([0,0,1])
            )
        else:
            raise ValueError(f"Unknown calculation method: {method}")

    #------------------   save data   --------------------#
    def save(self, outputFile=None, method="adaptive"):
        """
        Saves the model's data to a pickle file for future use.
        
        Parameters
        ----------
        outputFile : str
            The output file path where the data will be saved. If None is provided, the outputFile that was specified at model creation is used.
        method : str
            The calculation method that was used ("adaptive", "mean" or "normal"). It determines which of the models attributes are included in the save.
        
        Returns
        -------
        None.
        """
        common_keys = (
            'inputFile', 'primary_protons', 'pCurrent', 'pos_size', 'dir_size',
            'particle', 'energy', 'dir_window', 'x_min', 'x_max', 'y_min',
            'y_max', 'x_range', 'y_range',
            'x_mesh', 'y_mesh', 'z_mesh', 'total_neutrons', 'total_weights',
            'window_weights', 'total_mean', 'relative_error', 'sigma'
            )
        if method == "mean":
            common_keys += (
                'mean_directions', 'mean_directions_window',
                'mean_directions_spher'
                )
        elif method == "adaptive":
            common_keys += (
                'mean_directions','mean_directions_spher','adaptive_dir',
                'adaptive_dir_spher', 'hist','hist_edges'
            )
        elif method != "normal":
            raise ValueError(f"Unknown save method: {method}")
        
        dataDict = {k: self.__dict__[k] for k in common_keys if k in self.__dict__}
        
        if outputFile is None:
            outputFile = self.outputFile

        with open(outputFile, 'xb') as f:
            pkl.dump(dataDict, f)

    #------------------   load data   --------------------#
    def load(self, inputFile):
        """
        Load saved brightness data into the Flat object.
        """
        with open(inputFile, 'rb') as f:
            dataDict = pkl.load(f)
    
        # Restore attributes to self
        for key, value in dataDict.items():
            setattr(self, key, value)
    
        return self

    #------------------ Unified plotting ------------------#
    @staticmethod
    def _centers_to_edges(centers, fallback_width):
        """Return cell edges whose midpoints are the supplied grid centers."""
        centers = np.asarray(centers, dtype=float)
        if centers.ndim != 1 or centers.size == 0:
            raise ValueError("Plot coordinates must be a non-empty 1D array")
        if centers.size == 1:
            half_width = fallback_width / 2
            return np.array([centers[0] - half_width,
                             centers[0] + half_width])

        differences = np.diff(centers)
        if np.any(differences <= 0):
            raise ValueError("Plot coordinates must be strictly increasing")
        midpoints = centers[:-1] + differences / 2
        return np.concatenate(([centers[0] - differences[0] / 2],
                               midpoints,
                               [centers[-1] + differences[-1] / 2]))

    def _plot_method(self, method):
        """Resolve the calculation method and its brightness direction."""
        if method is None:
            method = getattr(self, "calculation_method", None)
        if method is None:
            # Compatibility with result files written before calculation_method
            # was stored.
            if hasattr(self, "adaptive_dir"):
                method = "adaptive"
            elif hasattr(self, "mean_directions"):
                method = "mean"
            else:
                method = "normal"

        if method == "adaptive":
            if not hasattr(self, "adaptive_dir"):
                raise ValueError("Adaptive directions are not available; run "
                                 "calculate(method='adaptive') first")
            directions = self.adaptive_dir
            directions_spher = self.adaptive_dir_spher
        elif method == "mean":
            if not hasattr(self, "mean_directions"):
                raise ValueError("Mean directions are not available; run "
                                 "calculate(method='mean') first")
            directions = self.mean_directions
            directions_spher = self.mean_directions_spher
        elif method == "normal":
            directions = np.tile((0.0, 0.0, 1.0), (self.x_mesh.size, 1))
            directions_spher = np.zeros((self.x_mesh.size, 2))
        else:
            raise ValueError(f"Unknown plotting method: {method}")
        return (
            method,
            np.asarray(directions, dtype=float), 
            np.asarray(directions_spher, dtype=float),
            )

    def _plot_flat_map(self, values, colorbar_label, method, show_arrows,
                       figsize, square_axes, arrow_scale):
        """Plot flat-mesh values and projected directions on the same grid."""
        if not np.isscalar(arrow_scale) or not np.isfinite(arrow_scale):
            raise ValueError("arrow_scale must be a finite positive number")
        if arrow_scale <= 0:
            raise ValueError("arrow_scale must be a finite positive number")

        _, directions, directions_spher = self._plot_method(method)
        shape = (len(self.y_range), len(self.x_range))
        values = np.asarray(values).reshape(shape)
        if directions.shape != (values.size, 3):
            raise ValueError("Direction and map grids have different sizes")

        # x_range/y_range are sample centers.  Supplying their extrema to
        # imshow as an extent treats those centers as outer edges and shifts
        # every heat-map cell away from its direction arrow.  Explicit edges
        # also keep non-uniform and one-cell grids aligned correctly.
        x_edges = self._centers_to_edges(self.x_range, self.pos_size_x / 4)
        y_edges = self._centers_to_edges(self.y_range, self.pos_size_y / 4)

        fig, ax = plt.subplots(figsize=figsize)
        mesh = ax.pcolormesh(x_edges, y_edges, values, cmap=plt.cm.viridis,
                             shading="flat")
        if square_axes:
            # Fill a square plotting box even when the x and y ranges differ
            # substantially.  This makes long, narrow surfaces easier to read.
            ax.set_aspect("auto")
            ax.set_box_aspect(1)
        else:
            # Preserve equal physical scale in x and y.  The axes box will then
            # follow the aspect ratio of the filtered data.
            ax.set_aspect("equal")

        # prepare_vectors stores the bounds of the filtered particle data.
        # The regular calculation grid can extend beyond those bounds, so crop
        # the axes to the actual data instead of displaying empty edge cells.
        x_limits = (getattr(self, "x_min", self.x_range[0]),
                    getattr(self, "x_max", self.x_range[-1]))
        y_limits = (getattr(self, "y_min", self.y_range[0]),
                    getattr(self, "y_max", self.y_range[-1]))
        plot_x_limits = (x_limits if x_limits[0] < x_limits[1]
                         else (x_edges[0], x_edges[-1]))
        plot_y_limits = (y_limits if y_limits[0] < y_limits[1]
                         else (y_edges[0], y_edges[-1]))
        ax.set_xlim(plot_x_limits)
        ax.set_ylim(plot_y_limits)

        if show_arrows:
            centers = np.column_stack((self.x_mesh.ravel(),
                                       self.y_mesh.ravel()))
            projected = directions[:, :2]
            # The arrow length is proportional to sin of the polar angle.
            # The arrow orientation reflects the azimuthal angle.
            projected_magnitude = np.sin(directions_spher[:, 1])
            
            visible = (
                np.isfinite(projected).all(axis=1)
                & np.isfinite(projected_magnitude)
                & (np.linalg.norm(projected, axis=1) > 0)
            )
            
            if np.any(visible):
                # Preserve the in-plane direction and the arrow magnitude
                # from the spherical polar angle.
                arrow_vectors = (
                    projected[visible]
                    * projected_magnitude[visible, np.newaxis]
                )
                if square_axes:
                    # In a stretched data coordinate system, screen-space
                    # angles keep the displayed arrows faithful to the x/y
                    # direction components.  Size them relative to the
                    # narrowest displayed grid-cell dimension.
                    x_fraction = (np.min(np.diff(x_edges))
                                  / np.diff(plot_x_limits)[0])
                    y_fraction = (np.min(np.diff(y_edges))
                                  / np.diff(plot_y_limits)[0])
                    arrow_fraction = (0.9 * arrow_scale
                                      * min(x_fraction, y_fraction))
                    quiver_options = dict(angles="uv", scale_units="width",
                                          scale=1 / arrow_fraction)
                else:
                    cell_size = min(np.min(np.diff(x_edges)),
                                    np.min(np.diff(y_edges)))
                    arrow_vectors *= 0.9 * arrow_scale * cell_size
                    quiver_options = dict(angles="xy", scale_units="xy",
                                          scale=1)
                ax.quiver(centers[visible, 0], centers[visible, 1],
                          arrow_vectors[:, 0], arrow_vectors[:, 1],
                          pivot="middle", color="blue", width=0.005,
                          headwidth=4, headlength=5, headaxislength=4.5,
                          **quiver_options)

        ax.set_xlabel('x [cm]', fontsize=18)
        ax.set_ylabel('y [cm]', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.ax.tick_params(labelsize=18)
        cbar.ax.set_ylabel(colorbar_label, fontsize=18)
        return fig, ax

    def plot_brightness_map(self, method=None, show_arrows=True,
                            figsize=(10, 8), square_axes=True,
                            arrow_scale=1.0):
        """Plot brightness and its maximizing direction at each mesh point.

        ``method`` can be ``"mean"``, ``"adaptive"``, ``"normal"``, or
        ``None`` to use the method from the latest calculation.
        ``figsize`` is the Matplotlib figure size in inches.
        Set ``square_axes=False`` to preserve equal physical x/y scaling.
        ``arrow_scale`` multiplies arrow length; larger values make the arrows
        longer. Its default is ``1.0``.

        Returns
        -------
        fig, ax
            The Matplotlib figure and axes for further customization.
        """
        if hasattr(self, "rotation"):
            print("Attention: This plot is displaying data in new coordinates compared to the input data. \
                The input data was transformed by applying the following rotation matrix R in the sense `R * (x,y,z)^T` (column vectors):\n"
                +str(self.rotation))
        return self._plot_flat_map(self.brightness,
                                   'brightness [n/s/cm$^2$/sr]',
                                   method, show_arrows, figsize, square_axes,
                                   arrow_scale)

    def plot_error_map(self, method=None, show_arrows=True, figsize=(10, 8),
                       square_axes=True, arrow_scale=1.0):
        """Plot relative error and the selected direction at each mesh point.

        ``figsize`` is the Matplotlib figure size in inches.
        Set ``square_axes=False`` to preserve equal physical x/y scaling.
        ``arrow_scale`` multiplies arrow length; larger values make the arrows
        longer. Its default is ``1.0``.

        Returns
        -------
        fig, ax
            The Matplotlib figure and axes for further customization.
        """
        if hasattr(self, "rotation"):
            print("Attention: This plot is displaying data in new coordinates compared to the input data. \
                The input data was transformed by applying the following rotation matrix R in the sense `R * (x,y,z)^T` (column vectors):\n"
                +str(self.rotation))
        return self._plot_flat_map(self.relative_error, 'relative error',
                                   method, show_arrows, figsize, square_axes,
                                   arrow_scale)
     
    def surface_crossing(self, v_x, v_y, v_z, theta_D):
        """
        Method to calculate neutron surface crossing for PHITS comparison
        w.r.t a given direction
        
        """
        data_filt = self.data_filter
        dir_D = np.array([v_x, v_y, v_z])
        slice_filt_i = data_filt[fc.cosRR_f(data_filt[vel], dir_D) >= np.cos(theta_D)]
        self.window_weights = fc.sum_f(slice_filt_i.Wt)
        
        # Compute statistical uncertainty
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)  # Global sum of squared weights
        Am = np.mean(self.data['Wt']**2)  # Mean squared weight
        self.sigma = self.compute_sigma(slice_filt_i, Wt_sq_sum, Am, self.window_weights)
    
            # Calculate relative error
        self.relative_error = self.r_err

#        dataDict = {'r_err': r_err, 'brightness': brightness, 'value_D': weight}
        return self.brightness, self.relative_error


 
    def surface_crossing_all_surface(self):
        """   
        Method to calculate surface crossing for the entire surface
        w.r.t mean direction
        
        """
        data_filt = self.data_filter
        
        # d(omega) = dS/r^2 = 3*6*1e-4/10^2 = 0.18
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        
        # Define the range of x and y values
        self.x_min, self.x_max = data_filt.x.min(), data_filt.x.max()
        self.y_min, self.y_max = data_filt.y.min(), data_filt.y.max()
        
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)  # Global sum of squared weights
        Am = np.mean(self.data['Wt']**2)  # Mean squared weight
        
        slice_i = data_filt
        self.total_neutrons = len(slice_i)
        self.total_weights = fc.sum_f(slice_i.Wt)
        self.total_mean = np.mean(slice_i.Wt)
        
        self.mean_directions = fc.direction_f(fc.mean_f(slice_i[vel], slice_i.Wt))
        
        slice_filt_i = slice_i[
            fc.cos_f(slice_i[vel], self.mean_directions) >= self.dir_window
            ]
        self.window_weights = fc.sum_f(slice_filt_i.Wt)
        
        self.sigma = self.compute_sigma(slice_i, Wt_sq_sum, Am, self.window_weights)
    
            # Calculate relative error
        self.relative_error = self.r_err
                        
        return self.brightness, self.relative_error     
