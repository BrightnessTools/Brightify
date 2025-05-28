

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

- **Input**: MCPL file with particle data (position, direction, energy, etc.), along with 6 other parameters defined by the user: type of the surface (flat, curved), type of the particle (neutron, proton, etc), number of scan points, position window size, direction window size and energy range.
- **Output**: Brightness map with optimal position and direction for neutron guide placement.
A flowchart is available in the paper published on Brightify. You can find the details in the citation section.


## 📦 Installation

Clone the repository and install Brightify:

    git clone https://github.com/yourusername/brightify.git
    cd brightify
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

First step:

    Import brightify 

An example of usage is detailed in /examples/verification.py. You can also find different MCPL files there in the `/examples` folder. These are examples which I have used for the arxiv paper which you can find below.

##  Citation

The paper for brightify is currently submitted and is waiting for the decision. Meanwhile, if you use **Brightify** in your research, please cite the arxiv version (DOI will be announced in near future):

> M. Akhyani, L. Zanini, H. M. Ronnow  
> **Brightify: A Tool for calculating brightness in neutron sources**.  
> arxiv, 2025.  
> DOI: 


##  License

This project is licensed under the **MIT License**. See LICENSE for details.



## 🤝 Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue or submit a pull request.



##  Acknowledgments

Developed as part of research at Laboratory of Qauntum Magnetism (LQM), EPFL, Lausanne, Switzerland.  
Supported by the **European Union’s Horizon 2020** research and innovation programme under the **Marie Skłodowska-Curie grant** agreement No. 754354.
