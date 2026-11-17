import unittest
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from brightify.main import Flat


class FlatPlottingTest(unittest.TestCase):
    def setUp(self):
        self.model = Flat(
            inputFile="./tests/data.mcpl.gz", # this holds 100 neutrons with pdg code 2112 (Dummy data for these tests here)
            outputFile=None,
            primary_protons=4.0,
            pCurrent=2.0,
            pos_size_x=4.0,
            pos_size_y=8.0,
            dir_size=1.0,
            )
        self.model.calculate(method="adaptive", bins_cart=10)

    def tearDown(self):
        plt.close("all")

    @patch("matplotlib.pyplot.show")
    def test_arrows_and_cells_share_grid_centers_for_direction_methods(self, _):
        for method in ("mean", "adaptive"):
            with self.subTest(method=method):
                self.model.plot_brightness_map(method=method)
                axes = plt.gcf().axes[0]
                cell_coordinates = axes.collections[0].get_coordinates()
                arrows = axes.collections[1]

                x_cell_centers = (
                    cell_coordinates[0, :-1, 0] + cell_coordinates[0, 1:, 0]
                ) / 2
                y_cell_centers = (
                    cell_coordinates[:-1, 0, 1] + cell_coordinates[1:, 0, 1]
                ) / 2
                np.testing.assert_allclose(x_cell_centers, self.model.x_range)
                np.testing.assert_allclose(y_cell_centers, self.model.y_range)
                directions = (
                    self.model.mean_directions if method == "mean"
                    else self.model.adaptive_dir
                )
                visible = (
                    np.isfinite(directions).all(axis=1)
                    & (np.linalg.norm(directions[:, :2], axis=1) > 0)
                )
                np.testing.assert_allclose(arrows.X,
                                           self.model.x_mesh.ravel()[visible])
                np.testing.assert_allclose(arrows.Y,
                                           self.model.y_mesh.ravel()[visible])
                self.assertEqual(arrows.pivot, "middle")
                self.assertTrue(np.all(np.isfinite(arrows.U)))
                self.assertTrue(np.all(np.isfinite(arrows.V)))
                np.testing.assert_array_equal(np.sign(arrows.U),
                                              np.sign(directions[visible, 0]))
                np.testing.assert_array_equal(np.sign(arrows.V),
                                              np.sign(directions[visible, 1]))
                expected_xlim = (
                    [self.model.x_min, self.model.x_max]
                    if self.model.x_min < self.model.x_max
                    else [cell_coordinates[0, 0, 0],
                          cell_coordinates[0, -1, 0]]
                )
                expected_ylim = (
                    [self.model.y_min, self.model.y_max]
                    if self.model.y_min < self.model.y_max
                    else [cell_coordinates[0, 0, 1],
                          cell_coordinates[-1, 0, 1]]
                )
                np.testing.assert_allclose(axes.get_xlim(), expected_xlim)
                np.testing.assert_allclose(axes.get_ylim(), expected_ylim)
                plt.close("all")

    @patch("matplotlib.pyplot.show")
    def test_custom_figure_size(self, _):
        self.model.plot_brightness_map(method="mean", show_arrows=False,
                                       figsize=(7, 5))
        np.testing.assert_allclose(plt.gcf().get_size_inches(), [7, 5])

    @patch("matplotlib.pyplot.show")
    def test_square_axes_are_default_and_can_be_disabled(self, _):
        self.model.plot_brightness_map(method="mean", show_arrows=False)
        axes = plt.gcf().axes[0]
        self.assertEqual(axes.get_box_aspect(), 1.0)

        plt.close("all")
        self.model.plot_brightness_map(method="mean", show_arrows=False,
                                       square_axes=False)
        axes = plt.gcf().axes[0]
        self.assertEqual(axes.get_aspect(), 1.0)

    @patch("matplotlib.pyplot.show")
    def test_arrow_scale_controls_length(self, _):
        self.model.plot_brightness_map(method="mean", arrow_scale=1.0)
        default_arrows = plt.gcf().axes[0].collections[1]
        default_scale = default_arrows.scale

        plt.close("all")
        self.model.plot_brightness_map(method="mean", arrow_scale=2.0)
        larger_arrows = plt.gcf().axes[0].collections[1]
        self.assertAlmostEqual(larger_arrows.scale, default_scale / 2)
        self.assertEqual(larger_arrows.width, 0.005)

    @patch("matplotlib.pyplot.show")
    def test_arrow_scale_must_be_positive(self, _):
        with self.assertRaisesRegex(ValueError, "arrow_scale"):
            self.model.plot_brightness_map(method="mean", arrow_scale=0)

    @patch("matplotlib.pyplot.show")
    def test_plot_methods_return_editable_figure_before_display(self, show):
        for plot_method in (self.model.plot_brightness_map,
                            self.model.plot_error_map):
            with self.subTest(plot_method=plot_method.__name__):
                fig, axes = plot_method(method="mean", show_arrows=False)
                self.assertIs(fig, plt.gcf())
                self.assertIs(axes, fig.axes[0])

                axes.set_xlim(0.5, 1.5)
                axes.set_ylim(10.5, 11.5)
                np.testing.assert_allclose(axes.get_xlim(), [0.5, 1.5])
                np.testing.assert_allclose(axes.get_ylim(), [10.5, 11.5])
                plt.close(fig)
        show.assert_not_called()

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
