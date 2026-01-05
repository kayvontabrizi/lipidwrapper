"""LipidWrapper - A tool for wrapping lipid bilayers around 3D surfaces.

Copyright (c) 2014, Jacob D. Durrant
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
"""

## imports

__all__ = [
    "clash_removal",
    "file_io",
    "hole_filling",
    "lipid_positioning",
    "molecule",
    "multiprocessing_utils",
    "numpy_extensions",
    "output",
]

# standard
import gc
import os
import pathlib
import platform
import shutil
import time
import typing

# custom
import typer

# local
from . import clash_removal
from . import file_io
from . import hole_filling
from . import lipid_positioning
from . import molecule
from . import multiprocessing_utils
from . import numpy_extensions
from . import output


## constants

version = "2.0.0"

app = typer.Typer(
    name="lipidwrapper",
    help="Wrap lipid bilayers around 3D surfaces.",
    add_completion=False,
)


## methods


@app.command()
def main(
    lipid_pdb_filename: pathlib.Path = typer.Option(
        ...,
        help="PDB file containing an all-atom model of a planar lipid bilayer.",
    ),
    lipid_headgroup_marker: str = typer.Option(
        "_P,CHL1_O3",
        help="Comma-separated list of atom specifications (RESNAME_ATOMNAME) identifying lipid headgroups.",
    ),
    surface_filename: str = typer.Option(
        "",
        help="Surface mesh file (PDB, DAE, or image file like PNG).",
    ),
    surface_equation: str = typer.Option(
        "z = 100*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)",
        help="Python equation defining z given x and y. Supports math, numpy, and scipy functions.",
    ),
    min_x: float = typer.Option(500, help="Minimum x value for mesh generation."),
    max_x: float = typer.Option(750, help="Maximum x value for mesh generation."),
    min_y: float = typer.Option(500, help="Minimum y value for mesh generation."),
    max_y: float = typer.Option(750, help="Maximum y value for mesh generation."),
    step_x: float = typer.Option(30, help="X-distance between adjacent mesh points."),
    step_y: float = typer.Option(30, help="Y-distance between adjacent mesh points."),
    max_height: float = typer.Option(
        0, help="Height of bilayer model at white regions (for image-based surfaces)."
    ),
    delete_clashing_lipids: bool = typer.Option(
        False, help="Delete lipids that sterically clash at triangle interfaces."
    ),
    clash_cutoff: float = typer.Option(
        2.0, help="Distance in Angstroms that constitutes a steric clash."
    ),
    clashing_potential_margin: float = typer.Option(
        25.0,
        help="Distance from triangle edges in Angstroms to check for clashes and holes.",
    ),
    fill_holes: bool = typer.Option(
        True, help="Fill holes left by deleting clashing lipids."
    ),
    fill_hole_exhaustiveness: int = typer.Option(
        10, help="How long LipidWrapper should try to fill holes."
    ),
    very_distant_lipids_cutoff: float = typer.Option(
        50.0,
        help="Skip clash checks for lipids further apart than this distance in Angstroms.",
    ),
    triangle_center_proximity_cutoff_distance: float = typer.Option(
        50.0,
        help="Distance cutoff for checking clashes between non-adjacent triangles.",
    ),
    memory_optimization_factor: int = typer.Option(
        1,
        help="Divide atom lists into chunks for pairwise comparisons to reduce memory usage.",
    ),
    number_of_processors: int = typer.Option(
        1, help="Number of processors to use for parallel processing."
    ),
    show_grid_points: bool = typer.Option(
        False, help="Append mesh grid points as atoms named 'X' to output."
    ),
    create_triangle_tcl_file: bool = typer.Option(
        False, help="Generate triangles.tcl file for VMD visualization."
    ),
    output_directory: str = typer.Option(
        "", help="Directory to save all output and intermediate files."
    ),
    use_disk_instead_of_memory: bool = typer.Option(
        False, help="Store growing model on disk instead of memory for large systems."
    ),
    compress_output: bool = typer.Option(
        False, help="Compress output files using gzip."
    ),
) -> None:
    if platform.system().lower() == "windows" and number_of_processors > 1:
        typer.echo(
            "REMARK WARNING: Use of multiple processors is only supported on Linux and OS X."
        )
        number_of_processors = 1

    if not lipid_pdb_filename.exists():
        typer.echo(
            f"ERROR: The file specified by --lipid-pdb-filename ({lipid_pdb_filename}) does not exist.\n"
        )
        raise typer.Exit(1)

    if output_directory and not output_directory.endswith("/"):
        output_directory = output_directory + "/"

    headgroup_markers = [
        (None, marker.strip().split("_")[0], None, marker.strip().split("_")[1])
        for marker in lipid_headgroup_marker.split(",")
    ]

    params = {
        "surface_filename": surface_filename,
        "surface_equation": surface_equation,
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y,
        "step_x": step_x,
        "step_y": step_y,
        "max_height": max_height,
        "lipid_pdb_filename": str(lipid_pdb_filename),
        "lipid_headgroup_marker": headgroup_markers,
        "show_grid_points": "TRUE" if show_grid_points else "FALSE",
        "create_triangle_tcl_file": "TRUE" if create_triangle_tcl_file else "FALSE",
        "delete_clashing_lipids": "TRUE" if delete_clashing_lipids else "FALSE",
        "use_disk_instead_of_memory": "TRUE" if use_disk_instead_of_memory else "FALSE",
        "clash_cutoff": clash_cutoff,
        "fill_holes": "TRUE" if fill_holes else "FALSE",
        "output_directory": output_directory,
        "fill_hole_exhaustiveness": fill_hole_exhaustiveness,
        "number_of_processors": number_of_processors,
        "clashing_potential_margin": clashing_potential_margin,
        "triangle_center_proximity_cutoff_distance": triangle_center_proximity_cutoff_distance,
        "memory_optimization_factor": memory_optimization_factor,
        "very_distant_lipids_cutoff": very_distant_lipids_cutoff,
        "compress_output": "TRUE" if compress_output else "FALSE",
        "memory_store_dir": output_directory + "store_in_memory.tmp/",
    }

    parameter_remarks = [
        "REMARK Parameters: (use the --help command-line parameter for further explanation)"
    ]
    for param_name, param_value in params.items():
        if param_name != "lipid_headgroup_marker":
            parameter_remarks.append(f"REMARK      {param_name}: {param_value}")
        else:
            parameter_remarks.append(
                f"REMARK      {param_name}: {lipid_headgroup_marker}"
            )
    parameter_remarks.append("")

    if not output_directory:
        typer.echo("\n".join(parameter_remarks))
    else:
        try:
            os.mkdir(output_directory)
        except FileExistsError:
            pass
        with open(output_directory + "parameters.input", "w") as parameter_file:
            parameter_file.write("\n".join(parameter_remarks))

    run_with_params(params)


def run_with_params(params: dict[str, typing.Any]) -> None:
    starttime = time.time()
    current_step = 0

    print("\nREMARK      LipidWrapper " + version + "\n")

    if params["use_disk_instead_of_memory"] == "TRUE":
        if os.path.exists(params["memory_store_dir"]):
            shutil.rmtree(params["memory_store_dir"])
        os.mkdir(params["memory_store_dir"])

    print("REMARK      Loading/creating and triangulating the mesh...")
    all_triangles = lipid_positioning.load_mesh_points_and_triangulations(params)

    print(
        "REMARK      Loading the original lipid-bilayer model ("
        + params["lipid_pdb_filename"]
        + ")..."
    )
    lipid, min_headgroups, max_headgroups = lipid_positioning.load_lipid_model(params)

    print("REMARK      Position copies of the lipid bilayer on the trianguled mesh...")
    molecules_by_triangle = (
        lipid_positioning.position_lipid_model_on_triangulated_tiles(
            params, lipid, all_triangles, min_headgroups, max_headgroups
        )
    )

    if params["output_directory"] != "":
        print("REMARK      Saving positioned lipid bilayers...")
        current_step = current_step + 1

        dir_pathname = (
            params["output_directory"]
            + "step_"
            + str(current_step)
            + ".positioned_lipid_triangles"
            + "/"
        )
        if not os.path.exists(dir_pathname):
            os.mkdir(dir_pathname)

        class save_positioned_lipids_multiprocessing(
            multiprocessing_utils.general_task
        ):
            def value_func(
                self, item: tuple, results_queue: typing.Optional[typing.Any]
            ) -> None:
                lipids = item[0]
                i = item[1]
                current_step = item[2]
                params = item[3]
                dir_pathname = item[4]

                self.print_star_if_appropriate(i)

                f = file_io.openfile(
                    dir_pathname
                    + "step_"
                    + str(current_step)
                    + ".original_positioned_lipid_triangle."
                    + str(i)
                    + ".pdb",
                    "w",
                    params,
                )

                if params["use_disk_instead_of_memory"] == "TRUE":
                    lipids = file_io.load_pickle(lipids, params)

                for lipid in lipids:
                    for index in range(len(lipid.all_atoms_numpy)):
                        f.write(lipid.create_pdb_line(index) + "\n")
                f.close()

        some_input = []
        i = 1
        gc.disable()
        for discard, lipids in molecules_by_triangle:
            i = i + 1
            some_input.append((lipids, i, current_step, params, dir_pathname))
        gc.enable()

        multiprocessing_utils.multi_threading(
            some_input,
            params["number_of_processors"],
            save_positioned_lipids_multiprocessing,
            params,
            "REMARK ",
        )

    if params["delete_clashing_lipids"] == "TRUE":
        print("REMARK      Deleting lipids that sterically clash...")
        clash_removal.remove_steric_clashes(molecules_by_triangle, params)

        if params["output_directory"] != "":
            current_step = current_step + 1
            dir_pathname = (
                params["output_directory"]
                + "step_"
                + str(current_step)
                + ".remove_lipids_with_clashes"
                + "/"
            )
            if not os.path.exists(dir_pathname):
                os.mkdir(dir_pathname)

            print(
                "REMARK            Saving the lipids that were not deleted for reference..."
            )

            class save_nondeleted_lipids_multiprocessing(
                multiprocessing_utils.general_task
            ):
                def value_func(
                    self, item: tuple, results_queue: typing.Optional[typing.Any]
                ) -> None:
                    triangle_index = item[0]
                    dir_pathname = item[1]
                    current_step = item[2]
                    lipids = item[3]

                    self.print_star_if_appropriate(triangle_index)

                    f = file_io.openfile(
                        dir_pathname
                        + "step_"
                        + str(current_step)
                        + ".retained_lipids_no_clash."
                        + str(triangle_index + 1)
                        + ".pdb",
                        "w",
                        params,
                    )
                    if params["use_disk_instead_of_memory"] == "TRUE":
                        triangle_lipids = file_io.load_pickle(lipids, params)
                    else:
                        triangle_lipids = lipids

                    for lipid in triangle_lipids:
                        for index in range(len(lipid.all_atoms_numpy)):
                            f.write(lipid.create_pdb_line(index) + "\n")
                    f.close()

            some_input = []
            gc.disable()
            for triangle_index in range(len(molecules_by_triangle)):
                some_input.append(
                    (
                        triangle_index,
                        dir_pathname,
                        current_step,
                        molecules_by_triangle[triangle_index][1],
                    )
                )
            gc.enable()
            multiprocessing_utils.multi_threading(
                some_input,
                params["number_of_processors"],
                save_nondeleted_lipids_multiprocessing,
                params,
                "REMARK ",
            )

        if params["fill_holes"] == "TRUE":
            print(
                "REMARK      Filling holes in the bilayer with additional lipid molecules..."
            )
            positioned_lipids_by_triangle = hole_filling.fill_in_lipid_holes(
                molecules_by_triangle, params
            )

            print("REMARK            Removing added lipids that clash...")
            clash_removal.remove_steric_clashes(positioned_lipids_by_triangle, params)

            if params["output_directory"] != "":
                print(
                    "REMARK            Saving the lipids that were added for reference..."
                )
                current_step = current_step + 1
                dir_pathname = (
                    params["output_directory"]
                    + "step_"
                    + str(current_step)
                    + ".lipid_holes_plugged"
                    + "/"
                )
                if not os.path.exists(dir_pathname):
                    os.mkdir(dir_pathname)

                class save_plugging_lipids_multiprocessing(
                    multiprocessing_utils.general_task
                ):
                    def value_func(
                        self, item: tuple, results_queue: typing.Optional[typing.Any]
                    ) -> None:
                        triangle_lipids = item[0]
                        index = item[1]
                        dir_pathname = item[2]
                        current_step = item[3]
                        params = item[4]

                        self.print_star_if_appropriate(index)

                        if params["use_disk_instead_of_memory"] == "TRUE":
                            triangle_lipids = file_io.load_pickle(
                                triangle_lipids, params
                            )

                        f = file_io.openfile(
                            dir_pathname
                            + "step_"
                            + str(current_step)
                            + ".lipids_added_into_bilayer_holes."
                            + str(index + 1)
                            + ".pdb",
                            "w",
                            params,
                        )
                        for lipid in triangle_lipids:
                            for i in range(len(lipid.all_atoms_numpy)):
                                f.write(lipid.create_pdb_line(i) + "\n")
                        f.close()

                some_input = []
                gc.disable()
                for ignore_var, lipids, index in positioned_lipids_by_triangle:
                    some_input.append(
                        (lipids, index, dir_pathname, current_step, params)
                    )
                gc.enable()

                multiprocessing_utils.multi_threading(
                    some_input,
                    params["number_of_processors"],
                    save_plugging_lipids_multiprocessing,
                    params,
                    "REMARK ",
                )

            print(
                "REMARK            Adding the hole-filling lipids to the original models..."
            )
            if params["use_disk_instead_of_memory"] == "TRUE":

                class add_plugging_lipids_multiprocessing(
                    multiprocessing_utils.general_task
                ):
                    def value_func(
                        self, item: tuple, results_queue: typing.Optional[typing.Any]
                    ) -> None:
                        existing_lipids_pickle_id = item[0]
                        position_lipids_pickle_id = item[1]
                        params = item[3]

                        self.print_star_if_appropriate(item[2])

                        tmp = file_io.load_pickle(existing_lipids_pickle_id, params)
                        tmp2 = file_io.load_pickle(position_lipids_pickle_id, params)
                        tmp.extend(tmp2)
                        file_io.save_pickle(tmp, params, existing_lipids_pickle_id)

                some_input = []
                gc.disable()
                for var_not_needed, lipids, index in positioned_lipids_by_triangle:
                    some_input.append(
                        (molecules_by_triangle[index][1], lipids, index, params)
                    )
                gc.enable()
                multiprocessing_utils.multi_threading(
                    some_input,
                    params["number_of_processors"],
                    add_plugging_lipids_multiprocessing,
                    params,
                    "REMARK ",
                )

            else:
                for var_not_needed, lipids, index in positioned_lipids_by_triangle:
                    molecules_by_triangle[index][1].extend(lipids)

            if params["output_directory"] != "":
                print(
                    "REMARK            Saving the bilayers with holes plugged for reference..."
                )

                class save_plugged_lipids_multiprocessing(
                    multiprocessing_utils.general_task
                ):
                    def value_func(
                        self, item: tuple, results_queue: typing.Optional[typing.Any]
                    ) -> None:
                        index = item[0]
                        dir_pathname = item[1]
                        current_step = item[2]
                        params = item[3]
                        lipids = item[4]

                        self.print_star_if_appropriate(index)

                        f = file_io.openfile(
                            dir_pathname
                            + "step_"
                            + str(current_step)
                            + ".all_lipids_with_holes_plugged."
                            + str(index + 1)
                            + ".pdb",
                            "w",
                            params,
                        )

                        if params["use_disk_instead_of_memory"] == "TRUE":
                            triangle_lipids = file_io.load_pickle(lipids, params)
                        else:
                            triangle_lipids = lipids

                        for lipid in triangle_lipids:
                            for i in range(len(lipid.all_atoms_numpy)):
                                f.write(lipid.create_pdb_line(i) + "\n")
                        f.close()

                some_input = []
                gc.disable()
                for index in range(len(molecules_by_triangle)):
                    some_input.append(
                        (
                            index,
                            dir_pathname,
                            current_step,
                            params,
                            molecules_by_triangle[index][1],
                        )
                    )
                gc.enable()
                multiprocessing_utils.multi_threading(
                    some_input,
                    params["number_of_processors"],
                    save_plugged_lipids_multiprocessing,
                    params,
                    "REMARK ",
                )

    print("REMARK      Printing out or saving all lipids to a single file...")
    if params["output_directory"] != "":
        current_step = current_step + 1
        dir_pathname = (
            params["output_directory"]
            + "step_"
            + str(current_step)
            + ".final_lipid_triangles"
            + "/"
        )
        if not os.path.exists(dir_pathname):
            os.mkdir(dir_pathname)
        f = file_io.openfile(
            params["output_directory"]
            + "step_"
            + str(current_step + 1)
            + ".full_bilayer.pdb",
            "w",
            params,
        )

        class save_final_lipids_multiprocessing(multiprocessing_utils.general_task):
            def value_func(
                self, item: tuple, results_queue: typing.Optional[typing.Any]
            ) -> None:
                params = item[0]
                lipids_lists = item[1]
                dir_pathname = item[2]
                current_step = item[3]
                triangle_index = item[4]

                self.print_star_if_appropriate(triangle_index)

                if params["use_disk_instead_of_memory"] == "TRUE":
                    triangle_lipids = file_io.load_pickle(lipids_lists, params)
                else:
                    triangle_lipids = lipids_lists

                g = file_io.openfile(
                    dir_pathname
                    + "step_"
                    + str(current_step)
                    + ".final_lipid_triangle."
                    + str(triangle_index + 1)
                    + ".pdb",
                    "w",
                    params,
                )

                for lipid in triangle_lipids:
                    for index in range(len(lipid.all_atoms_numpy)):
                        g.write(lipid.create_pdb_line(index) + "\n")

                g.close()

        some_input = []
        gc.disable()
        for triangle_index in range(len(molecules_by_triangle)):
            some_input.append(
                (
                    params,
                    molecules_by_triangle[triangle_index][1],
                    dir_pathname,
                    current_step,
                    triangle_index,
                )
            )
        gc.enable()
        multiprocessing_utils.multi_threading(
            some_input,
            params["number_of_processors"],
            save_final_lipids_multiprocessing,
            params,
            "REMARK ",
        )

    atomindex = 0
    resindex = 0

    for triangle_index in range(len(molecules_by_triangle)):

        if params["use_disk_instead_of_memory"] == "TRUE":
            triangle_lipids = file_io.load_pickle(
                molecules_by_triangle[triangle_index][1], params
            )
        else:
            triangle_lipids = molecules_by_triangle[triangle_index][1]

        for lipid in triangle_lipids:
            resindex = resindex + 1
            for index in range(len(lipid.all_atoms_numpy)):
                atomindex = atomindex + 1
                if params["output_directory"] != "":
                    f.write(lipid.create_pdb_line(index, atomindex, resindex) + "\n")
                else:
                    print(lipid.create_pdb_line(index, atomindex, resindex))

    if params["output_directory"] != "":
        f.close()

    if params["show_grid_points"] == "TRUE":
        print("REMARK      Printing out or saving the grid points for reference...")
        output.print_out_mesh_points(all_triangles, params)
    if params["output_directory"] != "" or params["create_triangle_tcl_file"] == "TRUE":
        print("REMARK      Creating a VMD TCL file showing the triangulations...")
        output.print_out_triangle_tcl_file(all_triangles, params)

    if params["use_disk_instead_of_memory"] == "TRUE":
        shutil.rmtree(params["memory_store_dir"])

    print("REMARK      Execution time: " + str(time.time() - starttime) + " seconds")
