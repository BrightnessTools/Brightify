import unittest
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Brightify.main import Flat


class FlatPlottingTest(unittest.TestCase):
    def setUp(self):
        self.model = Flat(pos_size_x=4.0, pos_size_y=8.0)
        self.model.x_range = np.array([0.0, 1.0, 2.0])
        self.model.y_range = np.array([10.0, 12.0])
        self.model.x_mesh, self.model.y_mesh = np.meshgrid(
            self.model.x_range, self.model.y_range
        )
        self.model.x_min, self.model.x_max = 0.25, 1.75
        self.model.y_min, self.model.y_max = 10.25, 11.75
        self.model.dir_size = 1.0
        self.model.pCurrent = 2.0
        self.model.primary_protons = 4.0
        self.model.window_weights = np.arange(1.0, 7.0)
        self.model.relative_error = np.arange(11.0, 17.0)
        self.model.mean_directions = np.tile([1.0, 1.0, 1.0], (6, 1))
        self.model.adaptive_dir = np.tile([-1.0, 1.0, 1.0], (6, 1))

    def tearDown(self):
        plt.close("all")

    @patch("matplotlib.pyplot.show")
    def test_arrows_and_cells_share_grid_centers_for_direction_methods(self, _):
        for method, expected_x_sign in (("mean", 1), ("adaptive", -1)):
            with self.subTest(method=method):
                self.model.plot_brightness_map(method=method)
                axes = plt.gcf().axes[0]
                cell_coordinates = axes.collections[0].get_coordinates()
                arrows = axes.collections[1]

                np.testing.assert_allclose(cell_coordinates[0, :, 0],
                                           [-0.5, 0.5, 1.5, 2.5])
                np.testing.assert_allclose(cell_coordinates[:, 0, 1],
                                           [9.0, 11.0, 13.0])
                np.testing.assert_allclose(arrows.X, self.model.x_mesh.ravel())
                np.testing.assert_allclose(arrows.Y, self.model.y_mesh.ravel())
                self.assertEqual(arrows.pivot, "middle")
                self.assertTrue(np.all(np.sign(arrows.U) == expected_x_sign))
                self.assertTrue(np.all(arrows.V > 0))
                np.testing.assert_allclose(axes.get_xlim(), [0.25, 1.75])
                np.testing.assert_allclose(axes.get_ylim(), [10.25, 11.75])
                plt.close("all")

    @patch("matplotlib.pyplot.show")
    def test_custom_figure_size(self, _):
        self.model.plot_brightness_map(method="mean", show_arrows=False,
                                       figsize=(7, 5))
        np.testing.assert_allclose(plt.gcf().get_size_inches(), [7, 5])

    @patch("matplotlib.pyplot.show")
    def test_square_axes_are_default_and_can_be_disabled(self, _):
        self.model.x_max = 10.0
        self.model.plot_brightness_map(method="mean", show_arrows=False)
        axes = plt.gcf().axes[0]
        self.assertEqual(axes.get_box_aspect(), 1.0)

        plt.close("all")
        self.model.plot_brightness_map(method="mean", show_arrows=False,
                                       square_axes=False)
        axes = plt.gcf().axes[0]
        self.assertEqual(axes.get_aspect(), 1.0)

    @patch("matplotlib.pyplot.show")
    def test_adaptive_error_map_displays_computed_errors(self, _):
        self.model.plot_error_map(method="adaptive", show_arrows=False)
        plotted = np.asarray(plt.gcf().axes[0].collections[0].get_array())
        np.testing.assert_allclose(plotted.ravel(), self.model.relative_error)

    @patch("matplotlib.pyplot.show")
    def test_all_methods_display_physical_brightness(self, _):
        for method in ("mean", "adaptive", "normal"):
            with self.subTest(method=method):
                self.model.plot_brightness_map(method=method,
                                               show_arrows=False)
                plotted = np.asarray(
                    plt.gcf().axes[0].collections[0].get_array()
                )
                np.testing.assert_allclose(plotted.ravel(),
                                           self.model.brightness)
                plt.close("all")


if __name__ == "__main__":
    unittest.main()
