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

def test_model_creation(flat_model):
    assert isinstance(flat_model, Flat)

def test_model_data(flat_model):
    assert len(getattr(flat_model, "data")) == 100

def test_model_particle_type(flat_model):
    assert flat_model.data["Kf"][0] == 2112

def test_model_filtering(flat_model):
    flat_model.apply_filter(energy=(0,1e-10))
    assert len(flat_model.data_filter) == 0