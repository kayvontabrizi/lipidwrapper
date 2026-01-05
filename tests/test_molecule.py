## imports

# standard
import math

# custom
import numpy
import numpy.testing
import pytest

# local
from lipidwrapper import molecule


class TestUnitVector:
    def test_unit_vector_x_axis(self):
        vector = numpy.array([5.0, 0.0, 0.0])
        result = molecule.unit_vector(vector)
        expected = numpy.array([1.0, 0.0, 0.0])
        numpy.testing.assert_array_almost_equal(result, expected)

    def test_unit_vector_arbitrary(self):
        vector = numpy.array([3.0, 4.0, 0.0])
        result = molecule.unit_vector(vector)
        expected = numpy.array([0.6, 0.8, 0.0])
        numpy.testing.assert_array_almost_equal(result, expected)

    def test_unit_vector_magnitude_is_one(self):
        vector = numpy.array([1.0, 2.0, 3.0])
        result = molecule.unit_vector(vector)
        magnitude = numpy.linalg.norm(result)
        assert abs(magnitude - 1.0) < 1e-10

    def test_unit_vector_3d(self):
        vector = numpy.array([1.0, 1.0, 1.0])
        result = molecule.unit_vector(vector)
        expected_magnitude = 1.0 / math.sqrt(3)
        expected = numpy.array([expected_magnitude] * 3)
        numpy.testing.assert_array_almost_equal(result, expected)


class TestAngleBetween:
    def test_perpendicular_vectors(self):
        v1 = numpy.array([1.0, 0.0, 0.0])
        v2 = numpy.array([0.0, 1.0, 0.0])
        result = molecule.angle_between(v1, v2)
        expected = numpy.pi / 2
        assert abs(result - expected) < 1e-10

    def test_parallel_vectors(self):
        v1 = numpy.array([1.0, 0.0, 0.0])
        v2 = numpy.array([2.0, 0.0, 0.0])
        result = molecule.angle_between(v1, v2)
        assert result == "NORMALIZED VECTORS EQUAL!"

    def test_opposite_vectors(self):
        v1 = numpy.array([1.0, 0.0, 0.0])
        v2 = numpy.array([-1.0, 0.0, 0.0])
        result = molecule.angle_between(v1, v2)
        assert abs(result - numpy.pi) < 1e-10

    def test_45_degree_angle(self):
        v1 = numpy.array([1.0, 0.0, 0.0])
        v2 = numpy.array([1.0, 1.0, 0.0])
        result = molecule.angle_between(v1, v2)
        expected = numpy.pi / 4
        assert abs(result - expected) < 1e-10

    def test_3d_vectors(self):
        v1 = numpy.array([1.0, 0.0, 0.0])
        v2 = numpy.array([0.0, 0.0, 1.0])
        result = molecule.angle_between(v1, v2)
        expected = numpy.pi / 2
        assert abs(result - expected) < 1e-10


class TestMolecule:
    def test_init(self):
        mol = molecule.Molecule()
        assert mol.in_triangle_margin is True
        assert mol.in_triangle_submargin is False
        assert mol.headgroup_index is None

    def test_load_pdb_from_lines(self, sample_pdb_lines):
        mol = molecule.Molecule()
        mol.load_pdb_from_lines(sample_pdb_lines)
        assert len(mol.all_atoms_numpy) == 5
        assert mol.all_atoms_numpy.shape == (5, 3)

    def test_load_pdb_coordinates(self, sample_molecule):
        numpy.testing.assert_array_almost_equal(
            sample_molecule.all_atoms_numpy[0], [0.0, 0.0, 0.0]
        )
        numpy.testing.assert_array_almost_equal(
            sample_molecule.all_atoms_numpy[1], [1.458, 0.0, 0.0]
        )

    def test_load_pdb_atom_names(self, sample_molecule):
        assert sample_molecule.atom_inf_string_vals[0][2] == "N"
        assert sample_molecule.atom_inf_string_vals[1][2] == "CA"
        assert sample_molecule.atom_inf_string_vals[2][2] == "C"

    def test_load_pdb_resnames(self, sample_molecule):
        for i in range(5):
            assert sample_molecule.atom_inf_string_vals[i][1] == "ALA"

    def test_load_pdb_resids(self, sample_molecule):
        for i in range(5):
            assert sample_molecule.atom_inf_resids[i] == "1"

    def test_copy_of(self, sample_molecule):
        mol_copy = sample_molecule.copy_of()
        numpy.testing.assert_array_equal(
            mol_copy.all_atoms_numpy, sample_molecule.all_atoms_numpy
        )
        mol_copy.all_atoms_numpy[0][0] = 999.0
        assert sample_molecule.all_atoms_numpy[0][0] != 999.0

    def test_portion_of(self, sample_molecule):
        indices = numpy.array([0, 2, 4])
        portion = sample_molecule.portion_of(indices)
        assert len(portion.all_atoms_numpy) == 3
        numpy.testing.assert_array_almost_equal(
            portion.all_atoms_numpy[0], sample_molecule.all_atoms_numpy[0]
        )
        numpy.testing.assert_array_almost_equal(
            portion.all_atoms_numpy[1], sample_molecule.all_atoms_numpy[2]
        )

    def test_set_undo_and_undo(self, sample_molecule):
        original_coords = sample_molecule.all_atoms_numpy.copy()
        sample_molecule.set_undo_point()
        sample_molecule.all_atoms_numpy[0] = [100.0, 100.0, 100.0]
        sample_molecule.undo()
        numpy.testing.assert_array_almost_equal(
            sample_molecule.all_atoms_numpy, original_coords
        )

    def test_get_indices_of_mask_match_by_atomname(self, sample_molecule):
        masks = [(None, "", None, "CA")]
        indices = sample_molecule.get_indices_of_mask_match(masks)
        assert len(indices) == 1
        assert indices[0] == 1

    def test_get_indices_of_mask_match_by_resname(self, sample_molecule):
        masks = [(None, "ALA", None, "")]
        indices = sample_molecule.get_indices_of_mask_match(masks)
        assert len(indices) == 5

    def test_get_indices_of_mask_match_combined(self, sample_molecule):
        masks = [(None, "ALA", None, "N")]
        indices = sample_molecule.get_indices_of_mask_match(masks)
        assert len(indices) == 1
        assert indices[0] == 0

    def test_get_indices_of_mask_match_no_match(self, sample_molecule):
        masks = [(None, "GLY", None, "")]
        indices = sample_molecule.get_indices_of_mask_match(masks)
        assert len(indices) == 0

    def test_create_pdb_line(self, sample_molecule):
        line = sample_molecule.create_pdb_line(0)
        assert "ATOM" in line
        assert "N" in line
        assert "ALA" in line

    def test_create_pdb_line_with_custom_index(self, sample_molecule):
        line = sample_molecule.create_pdb_line(0, output_index=100)
        assert "100" in line

    def test_baseN_decimal(self, sample_molecule):
        result = sample_molecule.baseN(10, 10)
        assert result == "10"

    def test_baseN_binary(self, sample_molecule):
        result = sample_molecule.baseN(10, 2)
        assert result == "1010"

    def test_baseN_hex(self, sample_molecule):
        result = sample_molecule.baseN(255, 16)
        assert result.lower() == "ff"

    def test_rotate_mol_quat(self, sample_molecule):
        from lipidwrapper import numpy_extensions

        angle = numpy.pi / 2
        q = numpy_extensions.Quaternion(
            numpy.cos(angle / 2), 0.0, 0.0, numpy.sin(angle / 2)
        )
        original_coords = sample_molecule.all_atoms_numpy.copy()
        sample_molecule.rotate_mol_quat(q)
        assert not numpy.allclose(
            sample_molecule.all_atoms_numpy[1], original_coords[1]
        )


class TestTriangle:
    def test_init(self, simple_triangle_points):
        tri = molecule.Triangle(simple_triangle_points)
        numpy.testing.assert_array_equal(tri.points, simple_triangle_points)

    def test_getitem(self, simple_triangle):
        numpy.testing.assert_array_equal(
            simple_triangle[0], numpy.array([0.0, 0.0, 0.0])
        )
        numpy.testing.assert_array_equal(
            simple_triangle[1], numpy.array([1.0, 0.0, 0.0])
        )

    def test_center(self, simple_triangle):
        center = simple_triangle.center()
        expected = numpy.array([0.5, 1.0 / 3.0, 0.0])
        numpy.testing.assert_array_almost_equal(center, expected)

    def test_center_equilateral(self, equilateral_triangle):
        center = equilateral_triangle.center()
        expected = numpy.array([0.5, numpy.sqrt(3) / 6, 0.0])
        numpy.testing.assert_array_almost_equal(center, expected)

    def test_radii(self, simple_triangle):
        radii = simple_triangle.radii()
        assert radii.shape == (1, 3)
        assert all(r > 0 for r in radii[0])

    def test_max_radius(self, simple_triangle):
        max_r = simple_triangle.max_radius()
        assert max_r > 0
        radii = simple_triangle.radii()
        assert max_r == numpy.max(radii)

    def test_project_points_onto_triangle_coplanar(self, simple_triangle):
        points = numpy.array([[0.25, 0.25, 0.0], [0.5, 0.5, 0.0]])
        projected = simple_triangle.project_points_onto_triangle(points)
        numpy.testing.assert_array_almost_equal(projected, points)

    def test_project_points_onto_triangle_above(self, simple_triangle):
        points = numpy.array([[0.5, 0.5, 5.0]])
        projected = simple_triangle.project_points_onto_triangle(points)
        expected = numpy.array([[0.5, 0.5, 0.0]])
        numpy.testing.assert_array_almost_equal(projected, expected)

    def test_get_indices_of_points_within_triangle_boundaries_inside(
        self, simple_triangle
    ):
        points = numpy.array([[0.5, 0.3, 0.0]])
        indices = simple_triangle.get_indices_of_points_within_triangle_boundaries(
            points
        )
        assert len(indices) == 1
        assert indices[0] == 0

    def test_get_indices_of_points_within_triangle_boundaries_outside(
        self, simple_triangle
    ):
        points = numpy.array([[5.0, 5.0, 0.0]])
        indices = simple_triangle.get_indices_of_points_within_triangle_boundaries(
            points
        )
        assert len(indices) == 0

    def test_get_indices_of_points_within_triangle_boundaries_mixed(
        self, simple_triangle
    ):
        points = numpy.array([[0.5, 0.3, 0.0], [5.0, 5.0, 0.0], [0.3, 0.2, 0.0]])
        indices = simple_triangle.get_indices_of_points_within_triangle_boundaries(
            points
        )
        assert len(indices) == 2
        assert 0 in indices
        assert 2 in indices

    def test_new_triangle_expanded_by_margin_positive(self, simple_triangle):
        expanded = simple_triangle.new_triangle_expanded_by_margin(0.5)
        assert expanded.max_radius() > simple_triangle.max_radius()

    def test_new_triangle_expanded_by_margin_negative(self, simple_triangle):
        contracted = simple_triangle.new_triangle_expanded_by_margin(-0.1)
        assert contracted.max_radius() < simple_triangle.max_radius()

    def test_new_triangle_expanded_by_margin_zero(self, simple_triangle):
        same = simple_triangle.new_triangle_expanded_by_margin(0.0)
        numpy.testing.assert_array_almost_equal(same.center(), simple_triangle.center())

    def test_near_other_triangle_adjacent(self, default_params):
        tri1_pts = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])
        tri2_pts = numpy.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 1.0, 0.0]])
        tri1 = molecule.Triangle(tri1_pts)
        tri2 = molecule.Triangle(tri2_pts)
        assert tri1.near_other_triangle(tri2, default_params) is True

    def test_near_other_triangle_distant(self, default_params):
        tri1_pts = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])
        tri2_pts = numpy.array(
            [[100.0, 100.0, 0.0], [101.0, 100.0, 0.0], [100.5, 101.0, 0.0]]
        )
        tri1 = molecule.Triangle(tri1_pts)
        tri2 = molecule.Triangle(tri2_pts)
        assert tri1.near_other_triangle(tri2, default_params) is False

    def test_degenerate_triangle_line(self):
        pts = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        tri = molecule.Triangle(pts)
        points = numpy.array([[0.5, 0.0, 0.0]])
        indices = tri.get_indices_of_points_within_triangle_boundaries(points)
        assert len(indices) == 0

    def test_degenerate_triangle_point(self):
        pts = numpy.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        tri = molecule.Triangle(pts)
        points = numpy.array([[0.0, 0.0, 0.0]])
        indices = tri.get_indices_of_points_within_triangle_boundaries(points)
        assert len(indices) == 0
