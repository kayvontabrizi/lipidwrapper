## imports

# standard
import math

# custom
import pytest
import numpy
import numpy.testing

# local
from lipidwrapper import numpy_extensions


class TestQuaternion:
    def test_init(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        assert q.v[0] == 1.0
        assert q.v[1] == 0.0
        assert q.v[2] == 0.0
        assert q.v[3] == 0.0

    def test_str(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        result = str(q)
        assert "1.0" in result
        assert "0.0" in result

    def test_copy_of(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        q_copy = q.copy_of()
        numpy.testing.assert_array_equal(q.v, q_copy.v)
        q_copy.v[0] = 999.0
        assert q.v[0] != q_copy.v[0]

    def test_normalize_identity(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        q_norm = q.normalize()
        magnitude = math.sqrt(sum(v**2 for v in q_norm.v))
        assert abs(magnitude - 1.0) < 1e-10

    def test_normalize_arbitrary(self):
        q = numpy_extensions.Quaternion(2.0, 3.0, 4.0, 5.0)
        q_norm = q.normalize()
        magnitude = math.sqrt(sum(v**2 for v in q_norm.v))
        assert abs(magnitude - 1.0) < 1e-10

    def test_add(self):
        q1 = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        q2 = numpy_extensions.Quaternion(5.0, 6.0, 7.0, 8.0)
        result = q1.add(q2)
        assert result.v[0] == 6.0
        assert result.v[1] == 8.0
        assert result.v[2] == 10.0
        assert result.v[3] == 12.0

    def test_minus(self):
        q1 = numpy_extensions.Quaternion(5.0, 6.0, 7.0, 8.0)
        q2 = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        result = q1.minus(q2)
        assert result.v[0] == 4.0
        assert result.v[1] == 4.0
        assert result.v[2] == 4.0
        assert result.v[3] == 4.0

    def test_scale(self):
        q = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        result = q.scale(2.0)
        assert result.v[0] == 2.0
        assert result.v[1] == 4.0
        assert result.v[2] == 6.0
        assert result.v[3] == 8.0

    def test_invert(self):
        q = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        result = q.invert()
        assert result.v[0] == 1.0
        assert result.v[1] == -2.0
        assert result.v[2] == -3.0
        assert result.v[3] == -4.0

    def test_multiply_identity(self):
        q1 = numpy_extensions.Quaternion(1.0, 0.0, 0.0, 0.0)
        q2 = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        result = q1.multiply(q2)
        numpy.testing.assert_array_almost_equal(result.v, q2.v)

    def test_multiply_associative(self):
        q1 = numpy_extensions.Quaternion(1.0, 2.0, 3.0, 4.0)
        q2 = numpy_extensions.Quaternion(5.0, 6.0, 7.0, 8.0)
        q3 = numpy_extensions.Quaternion(9.0, 10.0, 11.0, 12.0)
        result1 = q1.multiply(q2).multiply(q3)
        result2 = q1.multiply(q2.multiply(q3))
        numpy.testing.assert_array_almost_equal(result1.v, result2.v)

    def test_to_matrix_identity(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        matrix = q.to_matrix()
        expected = numpy.eye(3)
        numpy.testing.assert_array_almost_equal(matrix, expected)

    def test_to_matrix_rotation(self, sample_rotation_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_rotation_quaternion_values["s"],
            sample_rotation_quaternion_values["x"],
            sample_rotation_quaternion_values["y"],
            sample_rotation_quaternion_values["z"],
        )
        matrix = q.to_matrix()
        assert matrix.shape == (3, 3)
        determinant = numpy.linalg.det(matrix)
        assert abs(determinant - 1.0) < 1e-10

    def test_rep_as_44_matrix(self, sample_quaternion_values):
        q = numpy_extensions.Quaternion(
            sample_quaternion_values["s"],
            sample_quaternion_values["x"],
            sample_quaternion_values["y"],
            sample_quaternion_values["z"],
        )
        matrix = q.rep_as_44_matrix()
        assert matrix.shape == (4, 4)

    def test_load_from_mat_identity(self):
        q = numpy_extensions.Quaternion(0.0, 0.0, 0.0, 0.0)
        identity_matrix = numpy.eye(3)
        q.load_from_mat(identity_matrix)
        q_norm = q.normalize()
        assert abs(q_norm.v[0]) > 0.99

    def test_load_from_mat_rotation_z(self):
        angle = numpy.pi / 2
        rotation_matrix = numpy.array(
            [
                [numpy.cos(angle), -numpy.sin(angle), 0],
                [numpy.sin(angle), numpy.cos(angle), 0],
                [0, 0, 1],
            ]
        )
        q = numpy_extensions.Quaternion(0.0, 0.0, 0.0, 0.0)
        q.load_from_mat(rotation_matrix)
        result_matrix = q.to_matrix()
        numpy.testing.assert_array_almost_equal(
            result_matrix, rotation_matrix, decimal=5
        )


class TestGetNumpySlice:
    def test_basic_slice(self):
        array = numpy.array([10, 20, 30, 40, 50])
        indices = numpy.array([0, 2, 4])
        result = numpy_extensions.get_numpy_slice(array, indices)
        expected = numpy.array([10, 30, 50])
        numpy.testing.assert_array_equal(result, expected)

    def test_empty_indices(self):
        array = numpy.array([10, 20, 30, 40, 50])
        indices = numpy.array([])
        result = numpy_extensions.get_numpy_slice(array, indices)
        assert len(result) == 0

    def test_2d_array(self):
        array = numpy.array([[1, 2], [3, 4], [5, 6]])
        indices = numpy.array([0, 2])
        result = numpy_extensions.get_numpy_slice(array, indices)
        expected = numpy.array([[1, 2], [5, 6]])
        numpy.testing.assert_array_equal(result, expected)

    def test_single_index(self):
        array = numpy.array([10, 20, 30, 40, 50])
        indices = numpy.array([2])
        result = numpy_extensions.get_numpy_slice(array, indices)
        expected = numpy.array([30])
        numpy.testing.assert_array_equal(result, expected)
