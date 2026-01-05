## imports

# standard
import gc
import math
import typing

# custom
import numpy
import scipy.spatial.distance


## methods


def angle_between(v1: numpy.ndarray, v2: numpy.ndarray) -> typing.Union[float, str]:
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    if numpy.array_equal(v1_u, v2_u):
        return "NORMALIZED VECTORS EQUAL!"
    if numpy.linalg.norm(v1_u - v2_u) < 1e-7:
        return "NORMALIZED VECTORS EQUAL!"

    if numpy.array_equal(v1_u, -v2_u):
        return numpy.pi
    if numpy.linalg.norm(v1_u + v2_u) < 1e-7:
        return numpy.pi

    angle = numpy.arccos(numpy.dot(v1_u, v2_u))

    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return numpy.pi

    return angle


def unit_vector(vector: numpy.ndarray) -> numpy.ndarray:
    return vector / numpy.linalg.norm(vector)


## classes


class Molecule:
    in_triangle_margin: bool
    in_triangle_submargin: bool
    headgroup_index: typing.Optional[int]
    all_atoms_numpy: numpy.ndarray
    atom_inf_string_vals: numpy.ndarray
    atom_inf_resids: numpy.ndarray

    def __init__(self) -> None:
        self.in_triangle_margin = True
        self.in_triangle_submargin = False
        self.headgroup_index = None

    def get_headgroup_index(self, lipid_headgroup_marker: list[tuple]) -> int:
        if self.headgroup_index == None:
            self.headgroup_index = self.get_indices_of_mask_match(
                lipid_headgroup_marker
            )[0]
        return self.headgroup_index

    def load_pdb(self, filename: str) -> None:
        file = open(filename, "r")
        lines = file.readlines()
        file.close()
        self.load_pdb_from_lines(lines)

    def load_pdb_from_lines(self, lines: list[str]) -> None:
        self.__init__()

        gc.disable()

        self.atom_inf_string_vals = numpy.empty((len(lines), 4), dtype="U9")
        self.atom_inf_resids = numpy.empty(len(lines), dtype="U4")
        self.all_atoms_numpy = numpy.empty((len(lines), 3))

        count = 0
        for t in range(0, len(lines)):
            line = lines[t]
            if len(line) >= 7:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    count = count + 1

                    self.all_atoms_numpy[t][0] = float(line[30:38])
                    self.all_atoms_numpy[t][1] = float(line[38:46])
                    self.all_atoms_numpy[t][2] = float(line[46:54])

                    resname = line[16:21].strip()
                    atomname = line[11:16].strip()

                    try:
                        resid = line[22:26].strip()
                    except:
                        resid = "0"

                    self.atom_inf_string_vals[t][0] = line[21:22].strip()
                    self.atom_inf_string_vals[t][1] = resname
                    self.atom_inf_string_vals[t][2] = atomname
                    self.atom_inf_string_vals[t][3] = resname + "_" + atomname

                    self.atom_inf_resids[t] = resid

        gc.enable()

        self.atom_inf_string_vals = self.atom_inf_string_vals[:count]
        self.atom_inf_resids = self.atom_inf_resids[:count]
        self.all_atoms_numpy = self.all_atoms_numpy[:count]

    def save_pdb(self, filename: str) -> None:
        file = open(filename, "w")
        for index in range(len(self.all_atoms_numpy)):
            file.write(self.create_pdb_line(index) + "\n")
        file.close()

    def set_undo_point(self) -> None:
        self.all_atoms_numpy_undo = numpy.copy(self.all_atoms_numpy)

    def undo(self) -> None:
        self.all_atoms_numpy = numpy.copy(self.all_atoms_numpy_undo)

    def rotate_mol_quat(self, rot_quat: "numpy_extensions.Quaternion") -> None:
        rot_mat = rot_quat.to_matrix()
        self.all_atoms_numpy = numpy.dot(self.all_atoms_numpy, rot_mat)

    def baseN(
        self,
        num: int,
        b: int,
        numerals: str = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ",
    ) -> str:
        return ((num == 0) and numerals[0]) or (
            self.baseN(num // b, b, numerals).lstrip(numerals[0]) + numerals[num % b]
        )

    def create_pdb_line(
        self,
        index: int,
        output_index: typing.Optional[int] = None,
        output_resid: typing.Optional[int] = None,
    ) -> str:
        if output_index is None:
            output_index_str = str(index)
        else:
            output_index_str = str(output_index)

        if len(output_index_str) >= 7:
            output_index_str = "******"

        if output_resid is None:
            output_resid_str = self.atom_inf_resids[index]
        else:
            output_resid_str = str(output_resid)

        if len(output_resid_str) >= 5:
            output_resid_str = self.baseN(int(output_resid_str) + 2373280, 62)

        output = "ATOM "
        output = (
            output
            + str(output_index_str).rjust(6)
            + self.atom_inf_string_vals[index][2].rjust(5)
            + self.atom_inf_string_vals[index][1].rjust(5)
            + self.atom_inf_string_vals[index][0].rjust(1)
            + output_resid_str.rjust(4)
        )
        output = output + ("%.3f" % self.all_atoms_numpy[index][0]).rjust(12)
        output = output + ("%.3f" % self.all_atoms_numpy[index][1]).rjust(8)
        output = output + ("%.3f" % self.all_atoms_numpy[index][2]).rjust(8)

        return output

    def copy_of(self) -> "Molecule":
        new = Molecule()

        new.atom_inf_string_vals = self.atom_inf_string_vals.copy()
        new.atom_inf_resids = self.atom_inf_resids.copy()
        new.all_atoms_numpy = self.all_atoms_numpy.copy()
        new.headgroup_index = self.headgroup_index

        return new

    def portion_of(self, list_of_indices: numpy.ndarray) -> "Molecule":
        new = Molecule()

        new.atom_inf_string_vals = self.atom_inf_string_vals[list_of_indices]
        new.atom_inf_resids = self.atom_inf_resids[list_of_indices]
        new.all_atoms_numpy = self.all_atoms_numpy[list_of_indices]

        return new

    def get_indices_of_mask_match(self, masks: list[tuple]) -> numpy.ndarray:
        indices = numpy.array([], dtype=int)

        for mask in masks:
            chain = mask[0]
            resname = mask[1]
            resid = mask[2]
            atomname = mask[3]

            # get all the indices of the ones that have the same resname
            if chain == "" or chain is None or chain == -9999:
                indices_of_ones_with_this_chain = numpy.array(
                    list(range(len(self.atom_inf_string_vals)))
                )  # so it can be anything
            else:
                indices_of_ones_with_this_chain = numpy.nonzero(
                    self.atom_inf_string_vals[:, 0] == chain
                )[0]

            if resname == "" or resname is None or resname == -9999:
                indices_of_ones_with_this_resname = numpy.array(
                    list(range(len(self.atom_inf_string_vals)))
                )  # so it can be anything
            else:
                indices_of_ones_with_this_resname = numpy.nonzero(
                    self.atom_inf_string_vals[:, 1] == resname
                )[0]

            if resid == "" or resid is None or resid == -9999 or resid == "-9999":
                indices_of_ones_with_this_resid = numpy.array(
                    list(range(len(self.atom_inf_resids)))
                )  # so it can be anything
            else:
                indices_of_ones_with_this_resid = numpy.nonzero(
                    self.atom_inf_resids == resid
                )[0]

            if atomname == "" or atomname is None or atomname == -9999:
                indices_of_ones_with_this_atomname = numpy.array(
                    list(range(len(self.atom_inf_string_vals)))
                )  # so it can be anything
            else:
                indices_of_ones_with_this_atomname = numpy.nonzero(
                    self.atom_inf_string_vals[:, 2] == atomname
                )[0]

            # the intersection is the one that has both
            indices_in_all = numpy.intersect1d(
                indices_of_ones_with_this_chain,
                indices_of_ones_with_this_resname,
                assume_unique=True,
            )
            indices_in_all = numpy.intersect1d(
                indices_in_all, indices_of_ones_with_this_atomname, assume_unique=True
            )
            indices_in_all = numpy.intersect1d(
                indices_in_all, indices_of_ones_with_this_resid, assume_unique=True
            )
            indices = numpy.union1d(indices, indices_in_all)

        return indices


class Triangle:
    points: numpy.ndarray

    def __init__(self, pts: numpy.ndarray) -> None:
        self.points = pts

    def __getitem__(self, index: int) -> numpy.ndarray:
        return self.points[index]

    def center(self) -> numpy.ndarray:
        try:
            return self.center_pt
        except:
            self.center_pt = numpy.average(self.points, 0)
            return self.center_pt

    def radii(self) -> numpy.ndarray:
        try:
            return self.radii_lengths
        except:
            self.radii_lengths = scipy.spatial.distance.cdist(
                numpy.array([self.center()]), self.points, "euclidean"
            )
            return self.radii_lengths

    def max_radius(self) -> float:
        try:
            return self.radius_length
        except:
            self.radius_length = numpy.amax(self.radii())
            return self.radius_length

    def project_points_onto_triangle(self, pts: numpy.ndarray) -> numpy.ndarray:
        AB = self.points[1] - self.points[0]
        AC = self.points[2] - self.points[0]
        normal_to_plane = numpy.cross(AB, AC)
        normal_to_plane = normal_to_plane / numpy.linalg.norm(normal_to_plane)

        return pts - (
            (numpy.transpose([numpy.dot(pts - self.points[1], normal_to_plane)]))
            * normal_to_plane
        )

    def get_indices_of_points_within_triangle_boundaries(
        self, pts: numpy.ndarray
    ) -> numpy.ndarray:
        if (
            numpy.array_equal(self.points[0], self.points[1])
            or numpy.array_equal(self.points[0], self.points[2])
            or numpy.array_equal(self.points[1], self.points[2])
        ):
            return numpy.array([], dtype=numpy.intp)

        if numpy.isnan(self.points).any():
            return numpy.array([], dtype=numpy.intp)

        tri_min = numpy.min(self.points, 0) - numpy.array([1e-6, 1e-6, 1e-6])
        tri_max = numpy.max(self.points, 0) + numpy.array([1e-6, 1e-6, 1e-6])

        pts_not_in_triangle = (pts < tri_min).any(1)
        pts_not_in_triangle = numpy.logical_or(
            pts_not_in_triangle, (pts > tri_max).any(1)
        )
        pts_potentially_in_triangle = numpy.logical_not(pts_not_in_triangle)

        indices_of_pts_potentially_in_triangle = numpy.nonzero(
            pts_potentially_in_triangle
        )[0]

        indices_to_keep = []
        for t in indices_of_pts_potentially_in_triangle:
            t_v1 = self.points[0] - pts[t]
            t_v2 = self.points[1] - pts[t]
            t_v3 = self.points[2] - pts[t]

            angle1 = angle_between(t_v1, t_v2)
            angle2 = angle_between(t_v1, t_v3)
            angle3 = angle_between(t_v2, t_v3)

            if (
                angle1 == "NORMALIZED VECTORS EQUAL!"
                or angle2 == "NORMALIZED VECTORS EQUAL!"
                or angle3 == "NORMALIZED VECTORS EQUAL!"
            ):
                continue

            if math.fabs(angle1 + angle2 + angle3 - 2 * math.pi) < 0.01:
                indices_to_keep.append(t)

        return numpy.array(indices_to_keep, dtype=numpy.intp)

    def new_triangle_expanded_by_margin(self, margin: float) -> "Triangle":
        if margin < 0.0:
            if numpy.array_equal(self.points[0], self.points[1]) and numpy.array_equal(
                self.points[0], self.points[2]
            ):
                return Triangle(self.points)

        side_center_1 = numpy.average(self.points[[0, 1]], 0)
        side_center_2 = numpy.average(self.points[[1, 2]], 0)
        side_center_3 = numpy.average(self.points[[0, 2]], 0)
        side_centers = numpy.array([side_center_1, side_center_2, side_center_3])

        center_to_center_vec = side_centers - self.center()
        old_lengths = numpy.apply_along_axis(numpy.linalg.norm, 1, center_to_center_vec)
        new_lengths = old_lengths + margin

        if new_lengths[0] < 0.0 or new_lengths[1] < 0.0 or new_lengths[2] < 0.0:
            return Triangle(numpy.array([self.center(), self.center(), self.center()]))

        center_to_center_vec_normalized = (center_to_center_vec.T / old_lengths).T
        center_to_center_vec_new_length = (
            self.center() + (center_to_center_vec_normalized.T * new_lengths).T
        )

        second_points_on_new_side = (
            center_to_center_vec_new_length + self.points - side_centers
        )

        def seg_intersect(
            a1: numpy.ndarray,
            a2: numpy.ndarray,
            b1: numpy.ndarray,
            b2: numpy.ndarray,
        ) -> numpy.ndarray:
            pt1 = a1
            vec1 = a2 - a1

            pt2 = b1
            vec2 = b2 - b1

            coeffs = numpy.vstack((vec2, -vec1)).T
            best_sol_all = numpy.linalg.lstsq(coeffs, pt1 - pt2, rcond=None)
            best_sol = best_sol_all[0]

            if best_sol_all[1][0] == 0.0:
                return vec1 * best_sol[1] + pt1
            else:
                close_pt1 = vec1 * best_sol[1] + pt1
                close_pt2 = vec2 * best_sol[0] + pt2
                return (close_pt1 + close_pt2) * 0.5

        pt1 = seg_intersect(
            center_to_center_vec_new_length[0],
            second_points_on_new_side[0],
            center_to_center_vec_new_length[1],
            second_points_on_new_side[1],
        )
        pt2 = seg_intersect(
            center_to_center_vec_new_length[0],
            second_points_on_new_side[0],
            center_to_center_vec_new_length[2],
            second_points_on_new_side[2],
        )
        pt3 = seg_intersect(
            center_to_center_vec_new_length[2],
            second_points_on_new_side[2],
            center_to_center_vec_new_length[1],
            second_points_on_new_side[1],
        )

        new_pts = numpy.vstack((pt1, pt2, pt3))

        return Triangle(new_pts)

    def near_other_triangle(self, other_triangle: "Triangle", params: dict) -> bool:
        dists_between_triangle_pts = scipy.spatial.distance.cdist(
            self.points, other_triangle.points
        )
        if True in (dists_between_triangle_pts == 0.0):
            return True

        dist_between_center_points = numpy.linalg.norm(
            self.center() - other_triangle.center()
        )
        if (
            dist_between_center_points
            < params["triangle_center_proximity_cutoff_distance"]
        ):
            return True

        return False
