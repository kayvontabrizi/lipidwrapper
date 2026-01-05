## imports

# standard
import gc
import math

# custom
import numpy
import scipy.spatial.distance


## methods


def angle_between(v1, v2):
    """Calculates the angle between two vectors

    Arguments:
    v1 -- A 1x3 numpy array, the first vector
    v2 -- A 1x3 numpy array, the second vector

    Returns:
    A float, the angle between the two vectors in radians.
    If the two vectors are equal, the string "NORMALIZED VECTORS EQUAL!" is returned instead.

    """

    # get the unit vectors
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    # if the unit vectors are the same, you'll get an error.
    # return this string instead.
    if numpy.array_equal(v1_u, v2_u):
        return "NORMALIZED VECTORS EQUAL!"
    if numpy.linalg.norm(v1_u - v2_u) < 1e-7:
        return "NORMALIZED VECTORS EQUAL!"

    # if two vectors are pointing in the opposite directory, just return pi
    # This check is needed because sometimes numpy.dot(v1_u, v2_u) is actually slightly more than -1.0, giving an error
    if numpy.array_equal(v1_u, -v2_u):
        return numpy.pi
    if numpy.linalg.norm(v1_u + v2_u) < 1e-7:
        return numpy.pi

    # calculate the angle
    angle = numpy.arccos(numpy.dot(v1_u, v2_u))

    # if there's an error, modify the output
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return numpy.pi

    return angle


def unit_vector(vector):
    """Take a vector and scales it so its length is 1.0

    Arguments:
    vector -- A 1x3 numpy array, the vector to be scaled

    Returns:
    A 1x3 numpy array, the scaled vector

    """

    return vector / numpy.linalg.norm(vector)


## classes


class Molecule:
    """Loads, saves, and manupulates molecuar models."""

    def __init__(self):
        """Initializes the variables of the Molecule class."""

        self.in_triangle_margin = True
        self.in_triangle_submargin = False
        self.headgroup_index = None

    def get_headgroup_index(self, lipid_headgroup_marker):
        """Get's the indices of the current molecule's headgroup

        Arguments:
        lipid_headgroup_marker -- A tuple of the form (chain, resname, resid, atomname) specifying the headgroup

        Returns:
        An integer, the index of this molecule's headgroup

        """

        if self.headgroup_index == None:
            self.headgroup_index = self.get_indices_of_mask_match(
                lipid_headgroup_marker
            )[
                0
            ]  # so calculate it only if it's never been calculated before
        return self.headgroup_index

    def load_pdb(self, filename: str):
        """Loads a PDB file into the current Molecule object from a file

        Arguments:
        filename -- a string, the name of the file to load

        """

        # Now load the file into a list
        file = open(filename, "r")
        lines = file.readlines()
        file.close()

        # load the molecule from the list
        self.load_pdb_from_lines(lines)

    def load_pdb_from_lines(self, lines: list):
        """Loads a PDB file into the current Molecule object from a list of PDB lines

        Arguments:
        lines -- a list, containing the PDB lines to load into the current object

        """

        self.__init__()

        gc.disable()  # because appending objects slows down code if garbage collection turned on

        # set up the numpy arrays to store the data
        self.atom_inf_string_vals = numpy.empty(
            (len(lines), 4), dtype="U9"
        )  # chain, resname, atomname, id_keys
        self.atom_inf_resids = numpy.empty(len(lines), dtype="U4")
        self.all_atoms_numpy = numpy.empty((len(lines), 3))

        # read in the data from the lines
        count = 0
        for t in range(0, len(lines)):
            line = lines[t]
            if len(line) >= 7:
                if (
                    line[0:4] == "ATOM" or line[0:6] == "HETATM"
                ):  # Load atom data (coordinates, etc.)
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

                    self.atom_inf_string_vals[t][0] = line[21:22].strip()  # chain
                    self.atom_inf_string_vals[t][1] = resname  # resname
                    self.atom_inf_string_vals[t][2] = atomname  # atomname
                    self.atom_inf_string_vals[t][3] = (
                        resname + "_" + atomname
                    )  # id_keys

                    self.atom_inf_resids[t] = resid

        gc.enable()

        # now resize the array, cutting out bottom parts that were never populated
        self.atom_inf_string_vals = self.atom_inf_string_vals[:count]
        self.atom_inf_resids = self.atom_inf_resids[:count]
        self.all_atoms_numpy = self.all_atoms_numpy[:count]

    def save_pdb(self, filename: str):
        """Saves data to a PDB file.

        Arguments:
        filename -- A string, the filename to be written.

        """

        toprint = ""

        file = open(filename, "w")
        for index in range(len(self.all_atoms_numpy)):
            file.write(self.create_pdb_line(index) + "\n")
        file.close()

    def set_undo_point(self):
        """Sets ("saves") the undo point of all atoms. Any subsequent manipulations of atomic coordinates can be "undone" by reseting to this configuration."""

        self.all_atoms_numpy_undo = numpy.copy(self.all_atoms_numpy)

    def undo(self):
        """Resets the coordinates of all atoms to those saved using the set_undo_point function."""

        self.all_atoms_numpy = numpy.copy(self.all_atoms_numpy_undo)

    def rotate_mol_quat(self, rot_quat):
        """Support function that rotates a molecule according to a rotation quaternion

        Arguments:
        mol -- A molecule to be rotated in matrix format
        rot_quat -- Quaternion to rotate molecule

        """

        rot_mat = rot_quat.to_matrix()
        self.all_atoms_numpy = numpy.dot(self.all_atoms_numpy, rot_mat)

    def baseN(
        self,
        num: int,
        b: int,
        numerals: str = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ",
    ):
        """Return the value of a number in another base

        Arguments:
        num -- An integer, the number in base 10.
        b -- An integer, the number of the new base
        numerals -- An optional string, containing the numerals to use in the new base

        Returns:
        A string, the representation of the original integer, now in the specified base

        """

        return ((num == 0) and numerals[0]) or (
            self.baseN(num // b, b, numerals).lstrip(numerals[0]) + numerals[num % b]
        )

    def create_pdb_line(self, index: int, output_index=None, output_resid=None):
        """Create a string formatted according to the PDB standard from the atomic information contained in this atom class.

        Arguments:
        index -- An integer, the index of the atom in the Molecule object.
        output_index -- An optional integer, the index to use in the PDB-line output. If not specified, index is used.
        output_resid -- An optional integer, the resid to use in the PDB-line output. If not specified, the existing resid is used.

        Returns:
        A string, formatted according to the PDB standard.

        """

        # use user-specified index if provided
        if output_index is None:
            output_index = str(index)
        else:
            output_index = str(output_index)

        # PDB format is fixed column, so if the index is too big just turn it into stars
        if len(output_index) >= 7:
            output_index = "******"

        # use the user-specified resid if provided
        if output_resid is None:
            output_resid = self.atom_inf_resids[index]
        else:
            output_resid = str(output_resid)

        # PDB format is fixed column, so if the resid is too big, switch over to a string identifier that is unique to each residue
        if (
            len(output_resid) >= 5
        ):  # you need to start using letters in the resid after 9999
            # 2383280 is "a001" in base 62.
            # so map 10000 to 2383280 and convert to base 62.
            output_resid = self.baseN(int(output_resid) + 2373280, 62)
            # max using this method is 35999 residues

        # create the PDB line
        output = "ATOM "
        output = (
            output
            + str(output_index).rjust(6)
            + self.atom_inf_string_vals[index][2].rjust(5)
            + self.atom_inf_string_vals[index][1].rjust(5)
            + self.atom_inf_string_vals[index][0].rjust(1)
            + output_resid.rjust(4)
        )
        output = output + ("%.3f" % self.all_atoms_numpy[index][0]).rjust(12)
        output = output + ("%.3f" % self.all_atoms_numpy[index][1]).rjust(8)
        output = output + ("%.3f" % self.all_atoms_numpy[index][2]).rjust(8)

        return output

    def copy_of(self):
        """Create a copy of the current molecule

        Returns:
        A Molecule object, a copy of the current one

        """

        new = Molecule()

        new.atom_inf_string_vals = self.atom_inf_string_vals.copy()
        new.atom_inf_resids = self.atom_inf_resids.copy()
        new.all_atoms_numpy = self.all_atoms_numpy.copy()
        new.headgroup_index = self.headgroup_index

        return new

    def portion_of(self, list_of_indices):
        """Get a portion of the current molecule

        Arguments:
        list_of_indices -- A numpy array of integers, corresponding to the atoms in the current Molecule that should be in the return portion

        Returns:
        A Molecule object, containing only the atoms specified in list_of_indices

        """

        new = Molecule()

        new.atom_inf_string_vals = self.atom_inf_string_vals[list_of_indices]
        new.atom_inf_resids = self.atom_inf_resids[list_of_indices]
        new.all_atoms_numpy = self.all_atoms_numpy[list_of_indices]

        return new

    def get_indices_of_mask_match(self, masks):
        """Get the indices of atoms that match queries (masks)

        Arguments:
        masks -- A list of tuples. Each tuple is a mask/query of the format (chain, resname, resid, atomname). '', None, or -9999 all mean wildcard. For example: [('A', 'CHL1', 57, 'O3'), ('B', '', 783, 'P')]

        Returns:
        A numpy array of integers, containing the indices of the atoms that match the mask

        """

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
    """A class describing a triangle in three dimensions"""

    def __init__(self, pts):
        """Create a Triangle object.

        Arguments:
        pts -- A 3x3 numpy array containing the points of the triangle

        """

        self.points = pts

    def __getitem__(self, index: int):
        """Get one of the triangel points

        Arguments:
        index -- An integer, the index of the point

        Returns:
        A (1x3) numpy array, the coordinates of the requested point

        """

        return self.points[index]

    def center(self):
        """Get the triangle center

        Returns:
        A (1x3) numpy array, the coordinates of the triangle center

        """

        try:
            return self.center_pt
        except:
            self.center_pt = numpy.average(self.points, 0)
            return self.center_pt

    def radii(self):
        """Get the set of distances from the trangle center to each of its points (i.e., the "radii" of the triangle)

        Returns:
        A 1x3 numpy array, containing the distances/radii

        """

        try:
            return self.radii_lengths
        except:
            self.radii_lengths = scipy.spatial.distance.cdist(
                numpy.array([self.center()]), self.points, "euclidean"
            )
            return self.radii_lengths

    def max_radius(self):
        """Returns a float, the maximum distance from the triangle center to any of the contituent points"""

        try:
            return self.radius_length
        except:
            self.radius_length = numpy.amax(self.radii())
            return self.radius_length

    def project_points_onto_triangle(self, pts):
        """Projects a series of points onto the triangle plane

        Arguments:
        pts -- An nx3 numpy array, containing the points to be projected

        Returns:
        An nx3 numpy array, the coordinates of the requested point now projected onto the triangle plane

        """

        # define the triangle plane
        AB = self.points[1] - self.points[0]
        AC = self.points[2] - self.points[0]
        normal_to_plane = numpy.cross(AB, AC)
        normal_to_plane = normal_to_plane / numpy.linalg.norm(
            normal_to_plane
        )  # so it's a unit vector now

        # project the headgroups onto the triangle plane
        return pts - (
            (numpy.transpose([numpy.dot(pts - self.points[1], normal_to_plane)]))
            * normal_to_plane
        )

    def get_indices_of_points_within_triangle_boundaries(self, pts):
        """For a set of points that are coplanar with the current triangle, identify the indices of the points that are within the triangle boundaries

        Arguments:
        pts -- An nx3 numpy array, containing the points to be evaluated

        Returns:
        A numpy array, containing the indices of the points that fall within the triangle boundaries

        """

        # if the triangle is really just a line or a point, return an empty list
        if (
            numpy.array_equal(self.points[0], self.points[1])
            or numpy.array_equal(self.points[0], self.points[2])
            or numpy.array_equal(self.points[1], self.points[2])
        ):
            return numpy.array([])

        # some error has occurred previously, return an empty list
        if numpy.isnan(self.points).any():
            return numpy.array([])

        # get bounding box
        tri_min = numpy.min(self.points, 0) - numpy.array(
            [1e-6, 1e-6, 1e-6]
        )  # subtract a little to avoid rounding errors
        tri_max = numpy.max(self.points, 0) + numpy.array([1e-6, 1e-6, 1e-6])

        # identify points that couldn't possibly be in the triangle because they're outside the box
        pts_not_in_triangle = (pts < tri_min).any(1)
        pts_not_in_triangle = numpy.logical_or(
            pts_not_in_triangle, (pts > tri_max).any(1)
        )
        pts_potentially_in_triangle = numpy.logical_not(pts_not_in_triangle)

        # get the indices of the ones that could possibly be inside the triangle
        indices_of_pts_potentially_in_triangle = numpy.nonzero(
            pts_potentially_in_triangle
        )[0]

        # verify which ones really are in the triangle
        indices_to_keep = []
        for t in indices_of_pts_potentially_in_triangle:

            # calculate three vectors from the triangle verticies to the projection
            t_v1 = self.points[0] - pts[t]
            t_v2 = self.points[1] - pts[t]
            t_v3 = self.points[2] - pts[t]

            # get the appropriate angles
            angle1 = angle_between(t_v1, t_v2)
            angle2 = angle_between(t_v1, t_v3)
            angle3 = angle_between(t_v2, t_v3)

            # sometimes, if a triangle is small and the comparison point is very far away,
            # two of the vectors can end up being the same, especially after normalization.
            # Inevitably, in this case the point is not in the triangle.
            # we should account for that.
            if (
                angle1 == "NORMALIZED VECTORS EQUAL!"
                or angle2 == "NORMALIZED VECTORS EQUAL!"
                or angle3 == "NORMALIZED VECTORS EQUAL!"
            ):
                continue

            if (
                math.fabs(angle1 + angle2 + angle3 - 2 * math.pi) < 0.01
            ):  # it's inside the triangle
                indices_to_keep.append(t)

        return numpy.array(indices_to_keep)

    def new_triangle_expanded_by_margin(self, margin: float):
        """Return a triangle that is larger or smaller than the current one. The triangle is not simply scaled. A new triangle is carefully constructed such that each of its edges and the edges of the original triangle are a user-specified distance apart.

        Arguments:
        margin -- A float, the distance between the edges of the new triangle and the corresponding edges of the original triangle

        Returns:
        A Triangle object, the new triangle

        """

        # first, if this triangle is already a single point and the margin is negative, you can't collapse it further
        if margin < 0.0:
            if numpy.array_equal(self.points[0], self.points[1]) and numpy.array_equal(
                self.points[0], self.points[2]
            ):
                return Triangle(self.points)

        # get the centers of each side
        side_center_1 = numpy.average(self.points[[0, 1]], 0)
        side_center_2 = numpy.average(self.points[[1, 2]], 0)
        side_center_3 = numpy.average(self.points[[0, 2]], 0)
        side_centers = numpy.array([side_center_1, side_center_2, side_center_3])

        # get the vectors from the triangle center to each side center
        center_to_center_vec = side_centers - self.center()
        old_lengths = numpy.apply_along_axis(numpy.linalg.norm, 1, center_to_center_vec)
        new_lengths = old_lengths + margin

        # sanity check. None of the new lengths should be negative.
        # if any of them are, just set the whole triangle to a point
        # at the old triangle's center
        if new_lengths[0] < 0.0 or new_lengths[1] < 0.0 or new_lengths[2] < 0.0:
            return Triangle(numpy.array([self.center(), self.center(), self.center()]))

        # extend (or contract) these vectors farther out
        center_to_center_vec_normalized = (center_to_center_vec.T / old_lengths).T
        center_to_center_vec_new_length = (
            self.center() + (center_to_center_vec_normalized.T * new_lengths).T
        )

        # get an additional point on the extended side of the triangle
        second_points_on_new_side = (
            center_to_center_vec_new_length + self.points - side_centers
        )

        def seg_intersect(a1, a2, b1, b2):
            """Identify (or approximate) the point where two lines in 3D space intersect

            Arguments:
            a1 -- A 3x1 numpy array, a point on the first line
            a2 -- A 3x1 numpy array, a second point on the first line
            b1 -- A 3x1 numpy array, a point on the second line
            b2 -- A 3x1 numpy array, a second on the second line

            Returns:
            A 1x3 numpy array, the point (or an approximation thereof when an exact solution is not possible) where the defined lines intersect

            """

            # first, define the lines from the provided points
            pt1 = a1
            vec1 = a2 - a1

            pt2 = b1
            vec2 = b2 - b1

            # now get the points on the lines that are closest to each other
            coeffs = numpy.vstack((vec2, -vec1)).T
            best_sol_all = numpy.linalg.lstsq(coeffs, pt1 - pt2, rcond=None)
            best_sol = best_sol_all[0]

            if (
                best_sol_all[1][0] == 0.0
            ):  # an exact solution because the lines intersect
                return vec1 * best_sol[1] + pt1
            else:  # return the average pt of the two points that are closest to each other
                close_pt1 = vec1 * best_sol[1] + pt1
                close_pt2 = vec2 * best_sol[0] + pt2

                return (close_pt1 + close_pt2) * 0.5  # return the average pt

        # get the corners of the new triangle
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

    def near_other_triangle(self, other_triangle, params: dict):
        """Determines if another triangle is near this one

        Arguments:
        other_triangle -- A Triangle object, the other triangle to be evaluated
        params -- A dictionary, the user-specified command-line parameters

        Returns:
        A boolean, True of the two triangles are near each other, False otherwise

        """

        # check if the two triangles share a corner
        dists_between_triangle_pts = scipy.spatial.distance.cdist(
            self.points, other_triangle.points
        )
        if True in (dists_between_triangle_pts == 0.0):
            return True  # so they are adjacent

        # check if the distance to their centers is too close as well
        dist_between_center_points = numpy.linalg.norm(
            self.center() - other_triangle.center()
        )
        if (
            dist_between_center_points
            < params["triangle_center_proximity_cutoff_distance"]
        ):
            return True

        return False
