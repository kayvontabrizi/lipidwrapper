## imports

# standard
import os

# custom
import numpy
import pytest

# local
from lipidwrapper import molecule
from lipidwrapper import output


## classes


class TestPrintOutMeshPoints:
    @pytest.fixture
    def sample_triangles(self):
        tri1_pts = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])
        tri2_pts = numpy.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 1.0, 0.0]])
        return [molecule.Triangle(tri1_pts), molecule.Triangle(tri2_pts)]

    def test_creates_grid_points_file(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_mesh_points(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "grid_points.pdb")
        assert os.path.exists(expected_file)

    def test_grid_points_file_content(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_mesh_points(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "grid_points.pdb")
        with open(expected_file, "r") as f:
            content = f.read()
        assert "ATOM" in content
        assert "X" in content

    def test_no_duplicate_points(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_mesh_points(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "grid_points.pdb")
        with open(expected_file, "r") as f:
            lines = f.readlines()
        assert len(lines) == 5

    def test_prints_to_stdout_when_no_directory(self, sample_triangles, capsys):
        params = {"output_directory": ""}
        output.print_out_mesh_points(sample_triangles, params)
        captured = capsys.readouterr()
        assert "ATOM" in captured.out


class TestPrintOutTriangleTclFile:
    @pytest.fixture
    def sample_triangles(self):
        tri1_pts = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])
        tri2_pts = numpy.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 1.0, 0.0]])
        return [molecule.Triangle(tri1_pts), molecule.Triangle(tri2_pts)]

    def test_creates_tcl_file(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_triangle_tcl_file(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "triangles.tcl")
        assert os.path.exists(expected_file)

    def test_tcl_file_has_draw_commands(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_triangle_tcl_file(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "triangles.tcl")
        with open(expected_file, "r") as f:
            content = f.read()
        assert "draw delete all" in content
        assert "draw color red" in content
        assert "draw triangle" in content

    def test_tcl_file_has_correct_number_of_triangles(
        self, sample_triangles, temp_directory
    ):
        params = {"output_directory": temp_directory}
        output.print_out_triangle_tcl_file(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "triangles.tcl")
        with open(expected_file, "r") as f:
            content = f.read()
        triangle_count = content.count("draw triangle")
        assert triangle_count == len(sample_triangles)

    def test_tcl_file_contains_coordinates(self, sample_triangles, temp_directory):
        params = {"output_directory": temp_directory}
        output.print_out_triangle_tcl_file(sample_triangles, params)
        expected_file = os.path.join(temp_directory, "triangles.tcl")
        with open(expected_file, "r") as f:
            content = f.read()
        assert "0.0" in content
        assert "1.0" in content
