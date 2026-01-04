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

# standard
import sys
import time
import gc
import os
import shutil

# local
from . import multiprocessing_utils
from . import numpy_extensions
from . import molecule
from . import cli
from . import file_io
from . import lipid_positioning
from . import clash_removal
from . import hole_filling
from . import output


## constants

version = "1.15"


## methods


def run_program(argv: list):
    starttime = time.time()  # to keep track of execution time

    current_step = 0  # used for output filenames

    print "\nREMARK      LipidWrapper " + version + "\n"

    # check for Tkinter
    try:
        import Tkinter

        print "REMARK      The Tkinter python module is available. You may prefer to"
        print "REMARK      use the LipidWrapper graphical user interface"
        print "REMARK      (lipidwrapperGUI.py).\n"
        del Tkinter
    except:
        pass  # GUI not available

    params = cli.get_commandline_parameters(argv)  # get the commandline parameters

    # if you're going to be storing the growing model on the disk, make the temporary directory
    if params["use_disk_instead_of_memory"] == "TRUE":
        if os.path.exists(params["memory_store_dir"]):
            shutil.rmtree(params["memory_store_dir"])
        os.mkdir(params["memory_store_dir"])

    # load mesh points and generate triangulation
    print "REMARK      Loading/creating and triangulating the mesh..."
    all_triangles = lipid_positioning.load_mesh_points_and_triangulations(
        params
    )  # get the triangulations

    # load in the user-specified planar bilayer model
    print "REMARK      Loading the original lipid-bilayer model (" + params[
        "lipid_pdb_filename"
    ] + ")..."
    lipid, min_headgroups, max_headgroups = lipid_positioning.load_lipid_model(
        params
    )  # get the lipid molecule object, properly centered on x-y plane, as well as bounding-box coordinates

    # fill the tessellated triangles with bilayer sections
    print "REMARK      Position copies of the lipid bilayer on the trianguled mesh..."
    molecules_by_triangle = (
        lipid_positioning.position_lipid_model_on_triangulated_tiles(
            params, lipid, all_triangles, min_headgroups, max_headgroups
        )
    )  # position the lipids on the triangles

    # save the bilayer sections if user requested
    if params["output_directory"] != "":
        print "REMARK      Saving positioned lipid bilayers..."
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

        class save_positioned_lipids_multiprocessing(multiprocessing_utils.general_task):
            """A class for saving the lipid molecules associated with each triangle"""

            def value_func(self, item, results_queue):
                """Save lipid molecules associated with a triangle

                Arguments:
                item -- A list or tuple, the input data required for the calculation
                results_queue -- A multiprocessing.Queue() object for storing the calculation output

                """

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

    # delete clashing lipids if user requested
    if params["delete_clashing_lipids"] == "TRUE":
        print "REMARK      Deleting lipids that sterically clash..."
        clash_removal.remove_steric_clashes(
            molecules_by_triangle, params
        )  # remove steric clashes between lipids of adjacent tiles

        # save work from this step if user requested
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

            # print out remaining lipids
            print "REMARK            Saving the lipids that were not deleted for reference..."

            class save_nondeleted_lipids_multiprocessing(
                multiprocessing_utils.general_task
            ):
                """A class for saving the lipid molecules that were not deleted"""

                def value_func(self, item, results_queue):
                    """Save lipid molecules that were not deleted

                    Arguments:
                    item -- A list or tuple, the input data required for the calculation
                    results_queue -- A multiprocessing.Queue() object for storing the calculation output

                    """

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

        # fill holes in bilayer if user requested
        if params["fill_holes"] == "TRUE":

            # fill the holes
            print "REMARK      Filling holes in the bilayer with additional lipid molecules..."
            positioned_lipids_by_triangle = hole_filling.fill_in_lipid_holes(
                molecules_by_triangle, params
            )

            # remove filling lipids that clash
            print "REMARK            Removing added lipids that clash..."
            clash_removal.remove_steric_clashes(positioned_lipids_by_triangle, params)

            # save work from this step if user requested
            if params["output_directory"] != "":
                print "REMARK            Saving the lipids that were added for reference..."
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
                    """A class for saving the lipid molecules that were placed in lipid holes"""

                    def value_func(self, item, results_queue):
                        """Save lipid molecules that were placed in lipid holes

                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output

                        """

                        triangle_lipids = item[0]
                        index = item[1]
                        dir_pathname = item[2]
                        current_step = item[3]
                        params = item[4]

                        self.print_star_if_appropriate(index)

                        if params["use_disk_instead_of_memory"] == "TRUE":
                            triangle_lipids = file_io.load_pickle(triangle_lipids, params)

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

            # now add the positioned ligands into the lipid list associated with the original triangle
            print "REMARK            Adding the hole-filling lipids to the original models..."
            if params["use_disk_instead_of_memory"] == "TRUE":

                class add_plugging_lipids_multiprocessing(
                    multiprocessing_utils.general_task
                ):
                    """A class for adding the lipid molecules used to plug lipid holes to the list of lipids associated with the relevant triangle"""

                    def value_func(self, item, results_queue):
                        """Add lipid molecules used to plug holes to their associated triangles

                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output

                        """

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

            else:  # your going to have to merge all the lists into the main one regardless, so I see no utility in using multiple processors
                for var_not_needed, lipids, index in positioned_lipids_by_triangle:
                    molecules_by_triangle[index][1].extend(lipids)

            # save work from this step if user requested
            if params["output_directory"] != "":
                # print out all lipids
                print "REMARK            Saving the bilayers with holes plugged for reference..."

                class save_plugged_lipids_multiprocessing(
                    multiprocessing_utils.general_task
                ):
                    """A class for saving the lipid molecules used to plug lipid holes"""

                    def value_func(self, item, results_queue):
                        """Save lipid molecules used to plug lipid holes

                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output

                        """

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

    # now print out all the final molecules
    print "REMARK      Printing out or saving all lipids to a single file..."
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
            """A class for saving the final lipid models"""

            def value_func(self, item, results_queue):
                """Save the final lipid models

                Arguments:
                item -- A list or tuple, the input data required for the calculation
                results_queue -- A multiprocessing.Queue() object for storing the calculation output

                """

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

    # print out single files
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
                    print lipid.create_pdb_line(index, atomindex, resindex)

    if params["output_directory"] != "":
        f.close()

    # optional output files
    if params["show_grid_points"] == "TRUE":
        print "REMARK      Printing out or saving the grid points for reference..."
        output.print_out_mesh_points(all_triangles, params)
    if params["output_directory"] != "" or params["create_triangle_tcl_file"] == "TRUE":
        print "REMARK      Creating a VMD TCL file showing the triangulations..."
        output.print_out_triangle_tcl_file(all_triangles, params)

    # if the disk was used instead of memory, delete the temporary directory
    if params["use_disk_instead_of_memory"] == "TRUE":
        shutil.rmtree(params["memory_store_dir"])

    # tell the user how long it took for the program to execute
    print "REMARK      Execution time: " + str(time.time() - starttime) + " seconds"


if __name__ == "__main__":
    run_program(sys.argv)
