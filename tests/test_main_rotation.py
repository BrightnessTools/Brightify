import pytest
import numpy as np

def test_apply_rotation_with_matrix(flat_model):
    # Create a 90-degree rotation matrix around the Z-axis
    # [ cos(theta), -sin(theta), 0 ]
    # [ sin(theta),  cos(theta), 0 ]
    # [     0,           0,      1 ]
    theta = np.pi / 2
    c, s = np.cos(theta), np.sin(theta)
    rot_matrix = np.array([
        [c, -s, 0.0],
        [s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ])
    
    # Grab initial coordinates of the first particle for comparison
    initial_x = flat_model.data.loc[0, 'x']
    initial_y = flat_model.data.loc[0, 'y']
    
    # Apply the rotation matrix
    flat_model.apply_rotation(matrix=rot_matrix)
    
    # Verify the rotation attribute was stored
    assert hasattr(flat_model, "rotation")
    np.testing.assert_array_almost_equal(flat_model.rotation, rot_matrix)
    
    # Check that coordinates transformed as expected: (x', y') = (-y, x) for 90 deg
    new_x = flat_model.data.loc[0, 'x']
    new_y = flat_model.data.loc[0, 'y']
    assert np.isclose(new_x, -initial_y, atol=1e-7)
    assert np.isclose(new_y, initial_x, atol=1e-7)

def test_apply_rotation_vectors(flat_model):
    # Test valid vector rotation (e.g., pointing Z along [0, 0, 1] with custom X reference)
    flat_model.apply_rotation(z_vector=[0, 0, 1], x_vector=[1, 0, 0])
    assert hasattr(flat_model, "rotation")
    assert flat_model.rotation.shape == (3, 3)

def test_apply_rotation_parallel_vectors_raises_error(flat_model):
    # Test that parallel z_vector and x_vector raise a ValueError
    with pytest.raises(ValueError, match="The reference vector cannot be parallel to the normal vector"):
        flat_model.apply_rotation(z_vector=[0, 0, 1], x_vector=[0, 0, 2])

def test_reset_rotation_undo(flat_model):
    initial_x_coords = flat_model.data['x'].copy()
    
    # Apply rotation
    flat_model.apply_rotation(z_vector=[0, 1, 1], x_vector=[1, 0, 0])
    assert hasattr(flat_model, "rotation")
    
    # Reset rotation with undo=True
    flat_model.reset_rotation(undo=True)
    
    # Verify rotation attribute is deleted and coordinates are restored
    assert not hasattr(flat_model, "rotation")
    np.testing.assert_array_almost_equal(flat_model.data['x'], initial_x_coords)

def test_rotation_resets_filter(flat_model):
    # Apply a filter first
    flat_model.apply_filter(energy=(0, 1.0))
    # Modify data_filter length or state expectation
    assert len(flat_model.data_filter) <= len(flat_model.data)
    
    # Applying rotation should reset the filter
    flat_model.apply_rotation(z_vector=[0, 0, 1], x_vector=[1, 0, 0])
    
    # After reset, data_filter should cover all rows again
    assert len(flat_model.data_filter) == len(flat_model.data)