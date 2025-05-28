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
                 pCurrent=None, pos_size=None, dir_size=None):
        """
        Initializes the BrightifyModel class.
        
        Parameters:
        - inputFile: Path to the input file containing particle data (optional).
        - outputFile: Path to the output file where results are saved (optional).
        - primary_protons: Number of primary protons (used for brightness calculation).
        - pCurrent: Particle current (used for brightness calculation).
        - pos_size: Size of the position window (used for brightness calculation).
        - dir_size: Size of the direction window (used for brightness calculation).
        """
        
        # Ensure that either inputFile or outputFile exists
        assert path.exists(inputFile) or path.exists(outputFile)
        
        self.inputFile = inputFile
        self.primary_protons = primary_protons
        self.pCurrent = pCurrent
        self.pos_size = pos_size
        self.dir_size = dir_size
        
        # If an input file is provided, load the particle data from the file
        if inputFile:
            # Check that provided values are of the correct type
            assert type(int(self.primary_protons)) == int
            assert type(self.pCurrent) == float
            assert type(pos_size) == float
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
                             slice_filter_func, points):
        """
        Common calculations for total weights, mean weights, directions, and errors.
        This function can be reused in both `Spherical` and `Flat` models.
        """
        # Calculate mean squared weight (Am) and sum of squared weights (Wt_sq_sum)
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)  # Global sum of squared weights
        Am = np.mean(self.data['Wt']**2)  # Mean squared weight
    
        self.total_neutrons = np.full(shape1D, np.nan)
        self.total_weights = np.full(shape1D, np.nan)
        self.window_weights = np.full(shape1D, np.nan)
        self.total_mean = np.full(shape1D, np.nan)
        self.relative_error = np.full(shape1D, np.nan)
        self.sigma = np.full(shape1D, np.nan)
        self.mean_directions = np.full(shape2D, np.nan)
        self.mean_directions_spher = np.full(shape2D, np.nan)
        self.mean_directions_window = np.full(shape2D, np.nan)
        
    
        # Loop through the points and perform calculations
        for i, vec in tqdm(enumerate(points), desc='common', total=shape1D, ascii=' #'):
            # Apply the filtering function for slice_i
            slice_i = slice_filter_func(data_filt, vec)  # Callback for filtering slice_i
            
            # Store basic calculated values
            self.total_neutrons[i] = len(slice_i)
            self.total_weights[i] = fc.sum_f(slice_i.Wt)
            self.total_mean[i] = np.mean(slice_i.Wt)
    
            # Calculate mean direction based on velocity
            self.mean_directions[i] = fc.direction_f(fc.mean_f(slice_i[vel], slice_i.Wt))
            
            self.mean_directions_spher[i] = fc.cart2spher(*self.mean_directions[i])
    
            # Apply direction filter
            slice_filt_i = slice_i[fc.cos_f(slice_i[vel], self.mean_directions[i]) >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Calculate sigma (relative error)
            self.sigma[i] = self.compute_sigma(slice_filt_i, Wt_sq_sum, Am, self.window_weights[i])
    
            # Calculate relative error
            self.relative_error[i] = self.r_err[i]
    
            # Calculate the weighted direction
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
                 pCurrent=None, pos_size=None, dir_size=None):
        """
        Initialize the Flat model, which is a subclass of BrightifyModel.
        
        Parameters
        ----------
        primary_protons : int
            number of primary protons.
        pCurrent : float
            proton current defined by user.
        pos_size : float
            the number of points along each axis of the mesh grid.
        dir_size : float
            direction window size defined by user.

        Returns
        -------
        None.

        """
        super().__init__(inputFile, outputFile, primary_protons, pCurrent, 
                         pos_size, dir_size)         
    
  
    def flat_filter_func(self, data_filt, vec):
        """
        Filtering function specific to the Flat subclass. 
        It filters `slice_i` based on the x and y window ranges defined by
        `pos_size`.
        """
        return data_filt[(data_filt['x'] >= vec[0] - np.sqrt(self.pos_size)/2) & 
                         (data_filt['x'] <  vec[0] + np.sqrt(self.pos_size)/2) & 
                         (data_filt['y'] >= vec[1] - np.sqrt(self.pos_size)/2) & 
                         (data_filt['y'] <  vec[1] + np.sqrt(self.pos_size)/2)]

    def calculate(self):
        data_filt = self.data_filter
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        
        # Define x and y ranges for meshgrid
        self.x_min, self.x_max = data_filt.x.min(), data_filt.x.max()
        self.y_min, self.y_max = data_filt.y.min(), data_filt.y.max()
        self.x_range = np.arange(self.x_min, self.x_max + 1, np.sqrt(self.pos_size)/4)  
        self.y_range = np.arange(self.y_min, self.y_max + 1, np.sqrt(self.pos_size)/4)
        self.x_mesh, self.y_mesh = np.meshgrid(self.x_range, self.y_range)
        self.z_mesh = np.full_like(self.x_mesh, 1.)
        
        self.vectors = np.stack((self.x_mesh.flatten(), self.y_mesh.flatten(),
                                 self.z_mesh.flatten()), axis=-1)

        # Call the common calculation method from BrightifyModel
        return self.calculate_properties(data_filt, self.vectors.shape[0],
                                         self.vectors.shape,
                                         self.flat_filter_func,
                                         self.vectors)
 
    def plot_brightness_map(self):
        """
         Method to plot the brightness map
         
        """
        # put brightness values in the XY mesh
        brightness_map = self.brightness.reshape(len(self.y_range),
                                                 len(self.x_range))
#        brightness_map = self.relative_error.reshape(len(self.y_range),
#                                         len(self.x_range))
        # Calculate circle centers and radii
        circle_centers = np.stack((self.x_mesh.flatten(), 
                                   self.y_mesh.flatten()), axis=-1)
        circle_radii = np.sin(self.mean_directions_spher[:, 1])
        
        # Create a 2D plot
        fig, ax = plt.subplots(figsize=(16, 12))
    
        x_min = np.min(circle_centers[:, 0]) - np.sqrt(self.pos_size) / 2
        x_max = np.max(circle_centers[:, 0]) + np.sqrt(self.pos_size) / 2
        y_min = np.min(circle_centers[:, 1]) - np.sqrt(self.pos_size) / 2
        y_max = np.max(circle_centers[:, 1]) + np.sqrt(self.pos_size) / 2
        
        # Plot the brightness map
        im = ax.imshow(brightness_map, cmap=plt.cm.viridis, origin='lower',
                       extent=[x_min, x_max, y_min, y_max]) # vmin=4e13
        
        # Plot circles with vectors inside each mesh
        for center, radius, vector in zip(circle_centers, circle_radii, self.mean_directions):
            ax.arrow(center[0], center[1], radius * vector[0], radius * vector[1],
                     head_width=0.1, head_length=0.2, fc='blue', ec='blue')
        
        # Customize plot labels and ticks
        ax.set_xlabel('x [cm]', fontsize=18)
        ax.set_ylabel('y [cm]', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        
        # Colorbar
        cbar = plt.colorbar(im, label='brightness')
        cbar.ax.tick_params(labelsize=18)
        cbar.ax.set_ylabel('brightness  [n/s/cm$^2$/sr]', fontsize=18)
        
        # Show the plot
        plt.show()
        
    def plot_error_map(self):
        """
         Method to plot the brightness map
         
        """
        # put error values in the XY mesh
        error_map = self.relative_error.reshape(len(self.y_range),
                                         len(self.x_range))
        # Calculate circle centers and radii
        circle_centers = np.stack((self.x_mesh.flatten(), 
                                   self.y_mesh.flatten()), axis=-1)
        circle_radii = np.sin(self.mean_directions_spher[:, 1])
        
        # Create a 2D plot
        fig, ax = plt.subplots(figsize=(16, 12))
    
        x_min = np.min(circle_centers[:, 0]) - self.pos_size / 2
        x_max = np.max(circle_centers[:, 0]) + self.pos_size / 2
        y_min = np.min(circle_centers[:, 1]) - self.pos_size / 2
        y_max = np.max(circle_centers[:, 1]) + self.pos_size / 2
        
        # Plot the brightness map
        im = ax.imshow(error_map, cmap=plt.cm.viridis, origin='lower',
                       extent=[x_min, x_max, y_min, y_max]) # vmin=4e13
        
        # Plot circles with vectors inside each mesh
        for center, radius, vector in zip(circle_centers, circle_radii, self.mean_directions):
            ax.arrow(center[0], center[1], radius * vector[0], radius * vector[1],
                     head_width=0.1, head_length=0.2, fc='blue', ec='blue')
        
        # Customize plot labels and ticks
        ax.set_xlabel('x [cm]', fontsize=18)
        ax.set_ylabel('y [cm]', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        
        # Colorbar
        cbar = plt.colorbar(im, label='relative error')
        cbar.ax.tick_params(labelsize=18)
        cbar.ax.set_ylabel('relative error', fontsize=18)
        
        # Show the plot
        plt.show()
        
    def save(self, outputFile):
        """
        Method to save data into a binary pickle file 
            
        """
        keys = ('inputFile', 'primary_protons', 'pCurrent', 'pos_size', 'dir_size',
                'particle', 'energy', 'dir_window', 'x_range', 'y_range', 
                'x_mesh', 'y_mesh', 'z_mesh', 'total_neutrons', 'total_weights',
                'window_weights', 'total_mean', 'relative_error', 'mean_directions',
                'mean_directions_window'
                )
        
        # Create a dictionary of relevant attributes
        dataDict = {k: self.__dict__[k] for k in keys}
        
        # Save the dictionary to a file using pickle
        with open(outputFile, 'wb') as f:
            pkl.dump(dataDict, f)

 
    def surface_crossing(self, v_x, v_y, v_z, theta_D):
        """
        Method to calculate neutron surface crossing for PHITS comparison
        
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


    def surface_crossing_mesh(self, outputFile):
        """
        Method to calculate neutron surface crossing for PHITS comparison
        w.r.t surface normal (z direction)
        including position mesh
        
        """
        
        data_filt = self.data_filter
        self.dir_window = 1 - self.dir_size / (2 * np.pi)
        
        # Define x and y ranges for meshgrid
        self.x_min, self.x_max = data_filt.x.min(), data_filt.x.max()
        self.y_min, self.y_max = data_filt.y.min(), data_filt.y.max()
        self.x_range = np.arange(self.x_min, self.x_max + 1, np.sqrt(self.pos_size)/4)  
        self.y_range = np.arange(self.y_min, self.y_max + 1, np.sqrt(self.pos_size)/4)
        self.x_mesh, self.y_mesh = np.meshgrid(self.x_range, self.y_range)
        self.z_mesh = np.full_like(self.x_mesh, 1.)
        
        points = np.stack((self.x_mesh.flatten(), self.y_mesh.flatten(),
                                 self.z_mesh.flatten()), axis=-1)

        # Call the common calculation method from BrightifyModel
        # Calculate mean squared weight (Am) and sum of squared weights (Wt_sq_sum)
        Wt_sq_sum = fc.sum_f(self.data['Wt']**2)  # Global sum of squared weights
        Am = np.mean(self.data['Wt']**2)  # Mean squared weight
        shape1D = points.shape[0]
        dir_D = np.array([0, 0, 1])

        self.window_weights = np.full(shape1D, np.nan)
        self.relative_error = np.full(shape1D, np.nan)
        self.sigma = np.full(shape1D, np.nan)
        
        # Loop through the points and perform calculations
        for i, vec in tqdm(enumerate(points), desc='common', total=shape1D, ascii=' #'):
            # Apply the filtering function for slice_i
            slice_i = self.flat_filter_func(data_filt, vec)  # Callback for filtering slice_i

    
            # Apply direction filter
            slice_filt_i = slice_i[fc.cosRR_f(slice_i[vel], dir_D) >= self.dir_window]
            self.window_weights[i] = fc.sum_f(slice_filt_i.Wt)
    
            # Calculate sigma (relative error)
            self.sigma[i] = self.compute_sigma(slice_filt_i, Wt_sq_sum, Am, self.window_weights[i])
    
            # Calculate relative error
            self.relative_error[i] = self.r_err[i]
    

        keys = ('inputFile', 'primary_protons', 'pCurrent', 'pos_size', 'dir_size',
                'particle', 'energy', 'dir_window', 'x_range', 'y_range', 
                'x_mesh', 'y_mesh', 'z_mesh',
                'window_weights', 'relative_error'
                )
        
        # Create a dictionary of relevant attributes
        dataDict = {k: self.__dict__[k] for k in keys}
        
        # Save the dictionary to a file using pickle
        with open(outputFile, 'wb') as f:
            pkl.dump(dataDict, f)

 
    def surface_crossing_all_surface(self, outputFile):
        """   
        Method to calculate surface crossing for the entire surface
        
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

