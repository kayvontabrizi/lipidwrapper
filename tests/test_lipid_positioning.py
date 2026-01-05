## imports

# custom
import numpy
import numpy.testing
import pytest

# local
from lipidwrapper import lipid_positioning
from lipidwrapper import numpy_extensions


## classes


class TestGetTransformationData:
    def test_identity_transformation(self):
        points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        center_dynamic, rot_quat, center_static = (
            lipid_positioning.get_transformation_data(points, points)
        )
        numpy.testing.assert_array_almost_equal(center_dynamic, center_static)

    def test_translation_only(self):
        static_points = numpy.array(
            [[10.0, 10.0, 10.0], [11.0, 10.0, 10.0], [10.0, 11.0, 10.0]]
        )
        dynamic_points = numpy.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        )
        center_dynamic, rot_quat, center_static = (
            lipid_positioning.get_transformation_data(static_points, dynamic_points)
        )
        expected_dynamic_center = numpy.array([1.0 / 3, 1.0 / 3, 0.0])
        expected_static_center = numpy.array([31.0 / 3, 31.0 / 3, 10.0])
        numpy.testing.assert_array_almost_equal(center_dynamic, expected_dynamic_center)
        numpy.testing.assert_array_almost_equal(center_static, expected_static_center)

    def test_returns_quaternion(self):
        points1 = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        points2 = numpy.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
        _, rot_quat, _ = lipid_positioning.get_transformation_data(points1, points2)
        assert isinstance(rot_quat, numpy_extensions.Quaternion)

    def test_transformation_data_components(self):
        static_points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]])
        dynamic_points = numpy.array(
            [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.5, 6.0, 5.0]]
        )
        result = lipid_positioning.get_transformation_data(
            static_points, dynamic_points
        )
        assert len(result) == 3
        assert result[0].shape == (3,)
        assert isinstance(result[1], numpy_extensions.Quaternion)
        assert result[2].shape == (3,)


class TestApplyTransformation:
    def test_apply_transformation_translation(self):
        from lipidwrapper import molecule

        mol = molecule.Molecule()
        mol.all_atoms_numpy = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        center_dynamic = numpy.array([0.5, 0.0, 0.0])
        rot_quat = numpy_extensions.Quaternion(1.0, 0.0, 0.0, 0.0)
        center_static = numpy.array([10.5, 10.0, 10.0])
        transform_data = (center_dynamic, rot_quat, center_static)
        lipid_positioning.apply_transformation(mol, transform_data)
        expected = numpy.array([[10.0, 10.0, 10.0], [11.0, 10.0, 10.0]])
        numpy.testing.assert_array_almost_equal(mol.all_atoms_numpy, expected)

    def test_apply_transformation_preserves_relative_positions(self):
        from lipidwrapper import molecule

        mol = molecule.Molecule()
        mol.all_atoms_numpy = numpy.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        )
        original_distances = []
        for i in range(3):
            for j in range(i + 1, 3):
                dist = numpy.linalg.norm(
                    mol.all_atoms_numpy[i] - mol.all_atoms_numpy[j]
                )
                original_distances.append(dist)
        center_dynamic = numpy.mean(mol.all_atoms_numpy, axis=0)
        rot_quat = numpy_extensions.Quaternion(1.0, 0.0, 0.0, 0.0)
        center_static = numpy.array([50.0, 50.0, 50.0])
        transform_data = (center_dynamic, rot_quat, center_static)
        lipid_positioning.apply_transformation(mol, transform_data)
        new_distances = []
        for i in range(3):
            for j in range(i + 1, 3):
                dist = numpy.linalg.norm(
                    mol.all_atoms_numpy[i] - mol.all_atoms_numpy[j]
                )
                new_distances.append(dist)
        numpy.testing.assert_array_almost_equal(original_distances, new_distances)

    def test_apply_transformation_with_rotation(self):
        from lipidwrapper import molecule

        mol = molecule.Molecule()
        mol.all_atoms_numpy = numpy.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        angle = numpy.pi / 2
        rot_quat = numpy_extensions.Quaternion(
            numpy.cos(angle / 2), 0.0, 0.0, numpy.sin(angle / 2)
        )
        center_dynamic = numpy.array([1.5, 0.0, 0.0])
        center_static = numpy.array([1.5, 0.0, 0.0])
        transform_data = (center_dynamic, rot_quat, center_static)
        lipid_positioning.apply_transformation(mol, transform_data)
        expected_y_values = [0.5, -0.5]
        for i, expected_y in enumerate(expected_y_values):
            assert abs(mol.all_atoms_numpy[i][0] - 1.5) < 1e-5
            assert abs(abs(mol.all_atoms_numpy[i][1]) - 0.5) < 1e-5


class TestLoadMeshPointsAndTriangulations:
    @pytest.fixture
    def equation_params(self, temp_directory):
        return {
            "surface_filename": "",
            "surface_equation": "z = 0",
            "min_x": 0,
            "max_x": 10,
            "min_y": 0,
            "max_y": 10,
            "step_x": 5,
            "step_y": 5,
            "output_directory": temp_directory,
        }

    def test_load_from_equation_returns_triangles(self, equation_params):
        triangles = lipid_positioning.load_mesh_points_and_triangulations(
            equation_params
        )
        assert len(triangles) > 0

    def test_load_from_equation_flat_surface(self, equation_params):
        triangles = lipid_positioning.load_mesh_points_and_triangulations(
            equation_params
        )
        for tri in triangles:
            for point in tri.points:
                assert abs(point[2]) < 1e-10

    def test_load_from_equation_within_bounds(self, equation_params):
        triangles = lipid_positioning.load_mesh_points_and_triangulations(
            equation_params
        )
        for tri in triangles:
            for point in tri.points:
                assert point[0] >= equation_params["min_x"] - 1e-10
                assert point[0] <= equation_params["max_x"] + 1e-10
                assert point[1] >= equation_params["min_y"] - 1e-10
                assert point[1] <= equation_params["max_y"] + 1e-10

    def test_load_from_equation_varying_surface(self, equation_params):
        equation_params["surface_equation"] = "z = x * 0.1"
        equation_params["min_x"] = 0
        equation_params["max_x"] = 100
        equation_params["min_y"] = 0
        equation_params["max_y"] = 100
        equation_params["step_x"] = 25
        equation_params["step_y"] = 25
        triangles = lipid_positioning.load_mesh_points_and_triangulations(
            equation_params
        )
        assert len(triangles) > 0
        all_points = []
        for tri in triangles:
            for point in tri.points:
                all_points.append(point)
        x_values = [p[0] for p in all_points]
        assert max(x_values) > min(x_values)

    def test_triangles_have_correct_structure(self, equation_params):
        from lipidwrapper import molecule

        triangles = lipid_positioning.load_mesh_points_and_triangulations(
            equation_params
        )
        for tri in triangles:
            assert isinstance(tri, molecule.Triangle)
            assert tri.points.shape == (3, 3)
