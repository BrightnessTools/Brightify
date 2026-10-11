

__version__ = '2.0.0-alpha.1'

# All that we define here (i.e. give a name),
# will be allowed to be accessed via brightify.name,
# as long as any import statement that contains brightify
# is done by the user in his code.

from brightify.main import BrightifyModel, Flat, Spherical # this allows to do e.g. brightify.Flat(...)
from brightify.main import read_data # this allows to do brightify.read_data(...)
load_model = BrightifyModel.load_model # This allows to do brightify.load_model(...)
