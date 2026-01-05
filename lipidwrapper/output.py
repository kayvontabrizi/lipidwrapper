## imports

# custom
import numpy


## methods


def print_out_mesh_points(all_triangles: list, params: dict) -> None:
    def create_pdb_line(array: numpy.ndarray, letter: str) -> str:
        output = "ATOM "
        output = (
            output
            + "0".rjust(6)
            + letter.rjust(5)
            + "0".rjust(4)
            + letter.rjust(2)
            + "0".rjust(4)
        )
        output = output + ("%.3f" % array[0]).rjust(12)
        output = output + ("%.3f" % array[1]).rjust(8)
        output = output + ("%.3f" % array[2]).rjust(8)
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
        print("\n".join(toprint))
    else:
        f = open(params["output_directory"] + "grid_points.pdb", "w")
        f.write("\n".join(toprint))
        f.close()


def print_out_triangle_tcl_file(all_triangles: list, params: dict) -> None:
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
