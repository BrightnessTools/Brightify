#!/usr/bin/env python3
# -*- coding: utf-8 -*-

path = './examples/uniform/'

from Brightify import Flat

flatModel = Flat(
    inputFile= path +'output.mcpl.gz',
    primary_protons =   47000000,
    pCurrent = 1.0,
    pos_size = 1.0,
    dir_size = 6e-3
    )

flatModel.apply_filter(
    particle = 'neutron',
    energy = (0.0, 2.0),
    )

flatModel.calculate()

flatModel.plot_brightness_map()

flatModel.plot_error_map()
