## imports

# standard
import tempfile
import os

# custom
import pytest
import numpy

# local
from lipidwrapper import molecule


## fixtures


@pytest.fixture
def sample_quaternion_values():
    return {"s": 1.0, "x": 0.0, "y": 0.0, "z": 0.0}


@pytest.fixture
def sample_rotation_quaternion_values():
    angle = numpy.pi / 4
    return {
        "s": numpy.cos(angle / 2),
        "x": 0.0,
        "y": 0.0,
        "z": numpy.sin(angle / 2),
    }


@pytest.fixture
def simple_triangle_points():
    return numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])


@pytest.fixture
def simple_triangle(simple_triangle_points):
    return molecule.Triangle(simple_triangle_points)


@pytest.fixture
def equilateral_triangle_points():
    return numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, numpy.sqrt(3) / 2, 0.0]]
    )


@pytest.fixture
def equilateral_triangle(equilateral_triangle_points):
    return molecule.Triangle(equilateral_triangle_points)


@pytest.fixture
def sample_pdb_lines():
    return [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C",
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C",
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O",
        "ATOM      5  CB  ALA A   1       1.986  -0.760   1.217  1.00  0.00           C",
    ]


@pytest.fixture
def sample_molecule(sample_pdb_lines):
    mol = molecule.Molecule()
    mol.load_pdb_from_lines(sample_pdb_lines)
    return mol


@pytest.fixture
def temp_directory():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir + "/"


@pytest.fixture
def default_params(temp_directory):
    return {
        "output_directory": temp_directory,
        "memory_store_dir": temp_directory,
        "compress_output": "FALSE",
        "lipid_headgroup_marker": [(None, "", None, "P")],
        "clashing_potential_margin": 25.0,
        "triangle_center_proximity_cutoff_distance": 50.0,
    }
