

# Brightify

**Brightify** is a Python tool designed for neutron source designers who use Monte Carlo simulation tools to evaluate and optimize neutron **brightness** — a key measure of phase-space density that varies with both **position** and **direction**. Brightify identifies optimal locations and directions for neutron guide placement by scanning a user-defined surface and locating where the **maximum brightness** occurs.



## 🧩 Is Brightify for You?

Use the checklist below to determine if Brightify fits your needs:

1.  Is the distribution of neutrons from your source non-uniform or asymmetric?
2.  Do you need a high-resolution brightness map, i.e. your instrument needs a smaller region than the entire moderator surface?
3.  Are you designing a low-energy neutron source where guides can be placed close to the moderator?
4.  Do you want a more flexible and intuitive alternative to point tally or surface current tally methods?

If you answered **yes** to **any** of these, Brightify is for you!


##  Background

In Monte Carlo codes like **MCNP** and **PHITS**, brightness is not calculated directly. Instead, users rely on:

- **Surface current tallies**, which compute angular current in a single direction (e.g., x, y, or z). These allow position-resolved brightness mapping but only in one fixed direction.
- **Point estimator tallies**, which require precise placement and collimator definition. These must often be placed far from the source, making them unsuitable for compact systems.

Both methods demand **a priori** knowledge of where brightness is highest — a guess that works only for uniform or symmetric sources. When the neutron distribution is complex or when a **smaller phase-space volume** is needed (e.g., in compact neutron sources), these methods fall short.

**Brightify** addresses these limitations by enabling full scanning of both **spatial** and **angular** domains. It evaluates brightness directly from **MCPL**-formatted particle lists, with no dependency on the original simulation software.



## 🚀 Key Features

- **MCPL-Compatible**  
  Works with any Monte Carlo simulation that can generate MCPL files (e.g., PHITS, MCNP, Geant4, OpenMC).
  
-  **2D Brightness Scanning**  
  Evaluates brightness over both position and direction.

- **Visualization-Ready**  
  Outputs data suitable for generating high-resolution brightness maps.


## 📥 Input & 📤 Output

- **Input**: MCPL file with particle data (position, direction, energy, etc.), along with 6 other parameters defined by the user: type of the surface (flat, curved), type of the particle (neutron, proton, etc), number of scan points, position window size, direction window size and energy range. For creating MCPL files, you are invited to follow the instructions on the [`official MCPL webpage`](https://mctools.github.io/mcpl/hooks/).
- **Output**: Brightness map with optimal position and direction for neutron guide placement.
A flowchart is available in the paper published on Brightify. You can find the details in the citation section.


## 📦 Installation

Make sure you have installed and upgraded pip and git.
It is recommended to install Brightify in a virtual environment to avoid conflict with other existing packages, projects, etc.
Clone the repository and install Brightify:

    git clone https://github.com/BrightnessTools/Brightify.git
    cd Brightify
    pip install .
Or manually install dependencies:

    pip install -r requirements.txt

##  Dependencies

Brightify relies on the following Python packages:

-   `numpy`
    
-   `pandas`
    
-   `matplotlib`
    
-   `tqdm`
    
-   [`mcpl`](https://github.com/mctools/mcpl)
    

These are listed in `requirements.txt` for easy setup.

##  How to Use

First step is to import Brightify in your python environment:

    Import brightify 

An example of usage is detailed in

> /examples/verification.py

For running this example, you will need an MCPL file which you can download from below:

> https://zenodo.org/records/15537265?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImZiMDZkM2M5LTdhNjEtNDc3Mi1iNGZhLWU3MTAzNzMxNTBkMSIsImRhdGEiOnt9LCJyYW5kb20iOiIyZGZjNzk1ZjczYTZkODkwNWVkZjc2MWJlMDdjNzdiMiJ9.1wo5IMfgg7P7pZ0UNy04Z-d89Mu02v4Ma-qtA_R4K_8M2b9ZiyuUwFBaFsgy03gH-Xy1qHkP155rpdpR3wbqkQ

If you do simulations with PHITS, the PHITS input file for creating of this MCPL file is included in 

> /examples/gaussian_xy_centric/PHITS_simulation_files

This corresponds to the case (b) in the brightify paper.
The other two forlders include other sample cases corresponding to the case (a) and case (c) in the same paper and if you want the MCPL files, you can drop me an email.

##  Citation

If you use **Brightify** in your research, please cite the paper:

> M. Akhyani, L. Zanini, H. M. Ronnow  
> **Brightify: A tool for calculating brightness in neutron sources**.  
> Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment,
2025, 171037, ISSN 0168-9002.
> https://doi.org/10.1016/j.nima.2025.17103


##  License

This project is licensed under the **MIT License**. See LICENSE for details.



## 🤝 Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue or submit a pull request.



##  Acknowledgments

Developed as part of research at Laboratory of Qauntum Magnetism (LQM), EPFL, Lausanne, Switzerland.  
Supported by the **European Union’s Horizon 2020** research and innovation programme under the **Marie Skłodowska-Curie grant** agreement No. 754354.
