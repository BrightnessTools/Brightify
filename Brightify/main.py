#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from os import path
import pickle as pkl
from tqdm import tqdm
import matplotlib.pyplot as plt
import mcpl
from . import FastComputation as fc


#%% Define particle types with their corresponding PDG codes
Particle = {
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
                 pCurrent=None, pos_size_x=None, pos_size_y=None, dir_size=None):
        """
        Initializes the BrightifyModel class.
        
        Parameters:
        - inputFile: Path to the input file containing particle data (optional).
        - outputFile: Path to the output file where results are saved (optional).
        - primary_protons: Number of primary protons (used for brightness calculation).
        - pCurrent: Particle current (used for brightness calculation).
        - pos_size_x: Horizontal size of the position window (used for brightness calculation).
        - pose_size_y: Vertical size of the position window.
        - dir_size: Size of the direction window (used for brightness calculation).
        """
        
        # Ensure that either inputFile or outputFile exists
#        assert path.exists(inputFile) or path.exists(outputFile)
        
        self.inputFile = inputFile
        self.primary_protons = primary_protons
        self.pCurrent = pCurrent
        self.pos_size_x = pos_size_x
        self.pos_size_y = pos_size_y
        self.pos_size = pos_size_x * pos_size_y
        self.dir_size = dir_size
        
        # If an input file is provided, load the particle data from the file
        if inputFile:
            # Check that provided values are of the correct type
            assert type(int(self.primary_protons)) == int
            assert type(self.pCurrent) == float
            assert type(pos_size_x) == float
            assert type(pos_size_y) == float
            assert type(dir_size) == float
            
# if the ascii version is used as input
#            with open(self.inputFile, 'r') as f:
#                 self.data = f.read().replace('D', 'E')
#                
#            self.data = pd.read_csv(
#                 StringIO(self.data), header=None, delim_whitespace=True, 
#                 names=('Kf', 'x', 'y', 'z', 'Vx', 'Vy', 'Vz', 'e', 'Wt',
#                        'time', 'c1', 'c2', 'c3', 'name', 'nocas', 'nobch', 'no'
#                        )
#                 )
            
            # Read the MCPL file (format used for Brightify input files)
            myfile=mcpl.MCPLFile(self.inputFile,blocklength=1e13)
            # Loop through each particle block in the MCPL file and store data
            # in a pandas DataFrame
            for p in myfile.particle_blocks:
                self.data = pd.DataFrame({
                    'Kf': p.pdgcode,  # Particle PDG code
                    'x': p.x,  # x-position
                    'y': p.y,  # y-position
                    'z': p.z,  # z-position
                    'Vx': p.ux,  # Velocity in x-direction
                    'Vy': p.uy,  # Velocity in y-direction
                    'Vz': p.uz,  # Velocity in z-direction
                    'e': p.ekin,  # Kinetic energy
                    'Wt': p.weight  # Weight of the particle
                })
                
            # Reset filter to include all data initially    
            self.reset_filter()
            
        # If an output file is provided, load previously saved data
        if outputFile:
            with open(outputFile, 'rb') as f:
                dataDict = pkl.load(f)
            
            for k,v in dataDict.items():
                self.__dict__[k] = v
        
        
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
            self.filter &= (self.data['Kf'] == Particle[particle])
        
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
        primary_protons : int
            number of primary protons. Can be found in standard output file of 
            MC code
        pCurrent : float
            proton current defined by user.
        pos_size : float
            position window size defined by user.
        dir_size : float
            direction window size defined by user.
        radius : float
            radius of the spherical shell of the dump source defined by user.
        spiral_points : int
            number of points for spiral scan.

        Returns
        -------
        None.

        """
        # Initialize the parent class (BrightifyModel) with relevant arguments
        super().__init__(inputFile, outputFile, primary_protons, pCurrent, 
                         pos_size, dir_size)
        
        # If input file is provided, validate and initialize additional parameters
        if inputFile:
            assert type(radius) == float
            
            self.radius = radius
            self.spiral_points = spiral_points
            self.spiral = Spherical.spiral_scan(spiral_points)
            
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


    def save(self, outputFile):
        """
        Saves the model's data to a pickle file for future use.
        
        Parameters
        ----------
        outputFile : str
            The output file path where the data will be saved.
        
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
        dataDict = {k: self.__dict__[k] for k in keys}
        
        # Save the data dictionary to a pickle file
        with open(outputFile, 'wb') as f:
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
        """
        super().__init__(inputFile, outputFile, primary_protons, pCurrent, 
                         pos_size_x, pos_size_y, dir_size)

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
    def save(self, outputFile, method="adaptive"):
        common_keys = (
            'inputFile', 'primary_protons', 'pCurrent', 'pos_size', 'dir_size',
            'particle', 'energy', 'dir_window', 'x_range', 'y_range',
            'x_mesh', 'y_mesh', 'z_mesh', 'total_neutrons', 'total_weights',
            'window_weights', 'total_mean', 'relative_error'
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
        dataDict = {k: self.__dict__[k] for k in common_keys}
        with open(outputFile, 'wb') as f:
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
    def plot_brightness_map(self, method=None, show_arrows=True):
        """
        Plot brightness map with optional auto-detection of method
        and arrows for direction vectors.
    
        method: "mean", "adaptive", "normal", or None (auto-detect)
        show_arrows: whether to draw direction arrows
        """
    
        # Ensure bounds exist
        if not hasattr(self, "x_min") or not hasattr(self, "x_max"):
            self.x_min, self.x_max = self.x_range.min(), self.x_range.max()
            self.y_min, self.y_max = self.y_range.min(), self.y_range.max()
    
        # Auto-detect method if not provided
        if method is None:
            if hasattr(self, "adaptive_dir"):
                method = "adaptive"
            elif hasattr(self, "mean_directions"):
                method = "mean"
            else:
                method = "normal"
    
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
            raise ValueError(f"Unknown plotting method: {method}")
    
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
        and arrows for direction vectors.
        """
    
        # Ensure bounds exist
        if not hasattr(self, "x_min") or not hasattr(self, "x_max"):
            self.x_min, self.x_max = self.x_range.min(), self.x_range.max()
            self.y_min, self.y_max = self.y_range.min(), self.y_range.max()
    
        # Auto-detect method if not provided
        if method is None:
            if hasattr(self, "adaptive_dir"):
                method = "adaptive"
            elif hasattr(self, "mean_directions"):
                method = "mean"
            else:
                method = "normal"
    
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
            raise ValueError(f"Unknown plotting method: {method}")
    
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


 
    def surface_crossing_all_surface(self, outputFile):
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
        
        
        keys = ('inputFile', 'primary_protons', 'pCurrent', 'dir_size',
        'particle', 'energy', 'dir_window', 'total_neutrons', 'total_weights',
        'window_weights', 'total_mean', 'relative_error', 'sigma', 'mean_directions'
        )
        
        dataDict = {k: self.__dict__[k] for k in keys}
        
        with open(outputFile, 'wb') as f:
            pkl.dump(dataDict, f)
        
        return dataDict, self.brightness         
