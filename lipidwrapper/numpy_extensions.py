## imports

# standard
import math

# custom
import numpy


## classes


class Quaternion:
    v: numpy.ndarray

    def __init__(self, s: float, x: float, y: float, z: float) -> None:
        self.v = numpy.empty(4)
        self.v[0] = s
        self.v[1] = x
        self.v[2] = y
        self.v[3] = z

    def __str__(self) -> str:
        return (
            ""
            + str(self.v[0])
            + "\t"
            + str(self.v[1])
            + "\t"
            + str(self.v[2])
            + "\t"
            + str(self.v[3])
        )

    def copy_of(self) -> "Quaternion":
        return Quaternion(self.v[0], self.v[1], self.v[2], self.v[3])

    def load_from_mat(self, m: numpy.ndarray) -> None:
        if m.shape[0] != 3 or m.shape[1] != 3:
            print("Could not load quaternion from matrix...size is not (3x3)")
            return

        # Check that matrix is orthogonal. m_T = m_inv
        if not numpy.array_equal(numpy.transpose(m), numpy.linalg.inv(m)):
            print("Load Quaternion error. Matrix is not orthogonal")
            return

        # Need to make sure that the matrix is special orthogonal
        if math.fabs(1 - numpy.linalg.det(m)) > 0.000001:  # Done for rounding errors
            print("Load Quaternion error.  Determinant is not 1")
            return

        # First calculate the sum of the diagonal elements
        t = m.trace()

        if t > 0:
            S = math.sqrt(t + 1.0) * 2
            self.v[0] = 0.25 * S
            self.v[1] = (m[2, 1] - m[1, 2]) / S
            self.v[2] = (m[0, 2] - m[2, 0]) / S
            self.v[3] = (m[1, 0] - m[0, 1]) / S
        elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
            S = math.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
            self.v[0] = (m[2, 1] - m[1, 2]) / S
            self.v[1] = 0.25 * S
            self.v[2] = (m[0, 1] + m[1, 0]) / S
            self.v[3] = (m[0, 2] + m[2, 0]) / S
        elif m[1, 1] > m[2, 2]:
            S = math.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
            self.v[0] = (m[0, 2] - m[2, 0]) / S
            self.v[1] = (m[0, 1] + m[1, 0]) / S
            self.v[2] = 0.25 * S
            self.v[3] = (m[2, 1] + m[1, 2]) / S
        else:
            S = math.sqrt(1.0) * 2
            self.v[0] = (m[1, 0] - m[0, 1]) / S
            self.v[1] = (m[0, 2] + m[2, 0]) / S
            self.v[2] = (m[2, 1] + m[1, 2]) / S
            self.v[3] = 0.25 * S

    def rep_as_44_matrix(self) -> numpy.ndarray:
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]

        return numpy.array(
            [
                [qw, qx, qy, qz],
                [-qx, qw, -qz, qy],
                [-qy, qz, qw, -qx],
                [-qz, -qy, qx, qw],
            ]
        )

    def to_matrix(self) -> numpy.ndarray:
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]
        return numpy.array(
            [
                [
                    1.0 - 2.0 * qy * qy - 2.0 * qz * qz,
                    2.0 * qx * qy - 2.0 * qz * qw,
                    2.0 * qx * qz + 2.0 * qy * qw,
                ],
                [
                    2.0 * qx * qy + 2.0 * qz * qw,
                    1.0 - 2.0 * qx * qx - 2.0 * qz * qz,
                    2.0 * qy * qz - 2.0 * qx * qw,
                ],
                [
                    2.0 * qx * qz - 2.0 * qy * qw,
                    2.0 * qy * qz + 2.0 * qx * qw,
                    1.0 - 2.0 * qy * qy - 2.0 * qx * qx,
                ],
            ]
        )

    def add(self, q2: "Quaternion") -> "Quaternion":
        return Quaternion(
            self.v[0] + q2.v[0],
            self.v[1] + q2.v[1],
            self.v[2] + q2.v[2],
            self.v[3] + q2.v[3],
        )

    def invert(self) -> "Quaternion":
        return Quaternion(self.v[0], -1 * self.v[1], -1 * self.v[2], -1 * self.v[3])

    def minus(self, q2: "Quaternion") -> "Quaternion":
        return Quaternion(
            self.v[0] - q2.v[0],
            self.v[1] - q2.v[1],
            self.v[2] - q2.v[2],
            self.v[3] - q2.v[3],
        )

    def multiply(self, q2: "Quaternion") -> "Quaternion":
        return Quaternion(
            self.v[0] * q2.v[0]
            - self.v[1] * q2.v[1]
            - self.v[2] * q2.v[2]
            - self.v[3] * q2.v[3],
            self.v[1] * q2.v[0]
            + self.v[0] * q2.v[1]
            + self.v[2] * q2.v[3]
            - self.v[3] * q2.v[2],
            self.v[0] * q2.v[2]
            - self.v[1] * q2.v[3]
            + self.v[2] * q2.v[0]
            + self.v[3] * q2.v[1],
            self.v[0] * q2.v[3]
            + self.v[1] * q2.v[2]
            - self.v[2] * q2.v[1]
            + self.v[3] * q2.v[0],
        )

    def normalize(self) -> "Quaternion":
        n = math.sqrt(
            math.pow(self.v[0], 2)
            + math.pow(self.v[1], 2)
            + math.pow(self.v[2], 2)
            + math.pow(self.v[3], 2)
        )

        return Quaternion(self.v[0] / n, self.v[1] / n, self.v[2] / n, self.v[3] / n)

    def scale(self, scalar: float) -> "Quaternion":
        return Quaternion(
            self.v[0] * scalar,
            self.v[1] * scalar,
            self.v[2] * scalar,
            self.v[3] * scalar,
        )


## methods


def get_numpy_slice(
    numpy_array: numpy.ndarray, indices: numpy.ndarray
) -> numpy.ndarray:
    try:
        return numpy_array[indices]
    except:
        if len(indices) == 0:
            return numpy.array([])
        else:
            print("Error!")
