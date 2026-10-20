from brightify.main import BrightifyModel

def test_model_creation():
    bfy_model = BrightifyModel(pos_size_x=1.0, pos_size_y=1.0)
    assert bfy_model is not None