## imports

# custom
import numpy


## methods


def print_out_mesh_points(all_triangles: list, params: dict):
    """Save the mesh points to a PDB file

    Arguments:
    all_triangles -- A list of Triangle objects, the tesselation/triangulation
    params -- A dictionary, the user-specified command-line parameters

    """

    def create_pdb_line(numpy_array, letter: str):
        """Create a string formatted according to the PDB standard from the atomic information contained in this atom class.

        Arguments:
        numpy_array -- A numpy array, containing the atom coordinates
        letter -- A string, which will serve as the atom name, residue name, chain, etc.

        Returns:
        A string, formatted according to the PDB standard.

        """

        output = "ATOM "
        output = (
            output
            + "0".rjust(6)
            + letter.rjust(5)
            + "0".rjust(4)
            + letter.rjust(2)
            + "0".rjust(4)
        )
        output = output + ("%.3f" % numpy_array[0]).rjust(12)
        output = output + ("%.3f" % numpy_array[1]).rjust(8)
        output = output + ("%.3f" % numpy_array[2]).rjust(8)
        output = output + letter.rjust(24)
        return output

    toprint = []
    point_already_shown = []
    for tile_triangle in all_triangles:
        key = (
            str(tile_triangle[0][0])
            + "_"
            + str(tile_triangle[0][1])
            + "_"
            + str(tile_triangle[0][2])
        )
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[0], "X"))

        key = (
            str(tile_triangle[1][0])
            + "_"
            + str(tile_triangle[1][1])
            + "_"
            + str(tile_triangle[1][2])
        )
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[1], "X"))

        key = (
            str(tile_triangle[2][0])
            + "_"
            + str(tile_triangle[2][1])
            + "_"
            + str(tile_triangle[2][2])
        )
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[2], "X"))
    if params["output_directory"] == "":
        print "\n".join(toprint)
    else:
        f = open(params["output_directory"] + "grid_points.pdb", "w")
        f.write("\n".join(toprint))
        f.close()


def print_out_triangle_tcl_file(all_triangles: list, params: dict):
    """Save the tesselation/triangulation to a TCL file

    Arguments:
    all_triangles -- A list of Triangle objects, the tesselation/triangulation
    params -- A dictionary, the user-specified command-line parameters

    """

    # draw triangles
    f = open(params["output_directory"] + "triangles.tcl", "w")
    f.write("draw delete all\n")
    f.write("draw color red\n")

    for triangle in all_triangles:
        ia_pt = triangle[0]
        ib_pt = triangle[1]
        ic_pt = triangle[2]

        f.write(
            "draw triangle {"
            + str(ia_pt[0])
            + " "
            + str(ia_pt[1])
            + " "
            + str(ia_pt[2])
            + "} {"
            + str(ib_pt[0])
            + " "
            + str(ib_pt[1])
            + " "
            + str(ib_pt[2])
            + "} {"
            + str(ic_pt[0])
            + " "
            + str(ic_pt[1])
            + " "
            + str(ic_pt[2])
            + "}\n"
        )

    f.close()
