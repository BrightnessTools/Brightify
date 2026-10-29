import pytest
from brightify.main import Flat, BrightifyModel
import pandas as pd

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


def test_model_missing_inputs():
    # Test that providing neither inputFile nor outputFile raises an error
    with pytest.raises(ValueError, match="both inputFile and outputFile are None"):
        Flat(inputFile=None, outputFile=None)
        
def test_flat_invalid_pos_size():
    # Test that a negative pos_size_x raises a ValueError
    with pytest.raises(ValueError, match="pos_size_x incorrect"):
        Flat(
            inputFile="./tests/data.mcpl.gz",
            primary_protons=1000,
            pCurrent=1.0,
            pos_size_x=-1.0,  # Invalid
            pos_size_y=1.0,
            dir_size=6e-3,
        )

def test_read_mcpl_data():
    df = BrightifyModel.read_mcpl_data("./tests/data.mcpl.gz")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Kf" in df.columns
    assert "e" in df.columns

def test_save_default_output(tmp_path, flat_model):
    # Set an output file path using pytest's temporary directory fixture
    out_file = tmp_path / "output.bfy"
    flat_model.outputFile = str(out_file)
    
    # Call save without passing an argument (should use self.outputFile)
    flat_model.save(outputFile=None)
    
    assert out_file.exists()
