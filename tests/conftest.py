### This file holds "pytest fixtures" for any other test file

import pytest
from brightify.main import Flat

@pytest.fixture
def flat_model():
    return Flat(
        inputFile="./tests/data.mcpl.gz", # this holds 100 neutrons with pdg code 2112
        outputFile=None,
        primary_protons=47_000_000,
        pCurrent=1.0,
        pos_size_x=1.0,
        pos_size_y=1.0,
        dir_size=6e-3,
    )