## imports

# standard
import sys
import getopt
import textwrap
import os
import platform


## methods


def get_commandline_parameters(argv: list):
    """Get the user-defined command-line parameters

    Returns:
    A dictionary, the user-specified command-line parameters

    """

    # first check if the user has requested the help file
    if "--help" in [t.lower() for t in argv]:
        print()
        print("The initial lipid model")
        print("=======================")
        print()
        print(
            textwrap.fill(
                "--lipid_pdb_filename: This parameter specifies a PDB file containing an all-atom model of a planar lipid bilayer. LipidWrapper will wrap this lipid around the user-generated mesh. Example: --lipid_pdb_filename lipid.pdb",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                '--lipid_headgroup_marker: A unique atom representing the headgroup of each lipid residue must be specified. The --lipid_headgroup_marker accepts a comma-separated lists of atom specifications (RESNAME_ATOMNAME). If either RESNAME or ATOMNAME is omitted, any value will be accepted. By default, LipidWrapper identifies lipid headgroups by looking for any atom named "P" (_P) or any atom named "O3" belonging to a cholesterol molecule (CHL1_O3). Example: --lipid_headgroup_marker "_P,CHL1_O3"',
                70,
                subsequent_indent="      ",
            )
        )
        print()
        print("Methods for creating a surface mesh")
        print("===================================")
        print()
        print(
            textwrap.fill(
                '--surface_equation: Generate a surface mesh from a python-formatted equation defining z, given x and y. The --min_x, --max_x, --min_y, and --max_y parameters are used to specify the region over which the function should be evaluated. The --step_x and --step_y parameters define the x-y distance between adjacent points. Python functions from the math, numpy, and scipy modules can be used. Example: --surface_equation "z = 250*numpy.sin(x*x/60000 +y*y/60000)"',
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--surface_filename: If this parameter specifies a file with the PDB extension, a surface mesh is generated from the coordinates of the PDB atoms. Example: --surface_filename mymesh.pdb",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--surface_filename: If this parameter specifies a file with the DAE extension, the mesh points and triangulations will be taken from the file. Example: --surface_filename mymodel.dae",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--surface_filename: If this parameter specifies a file that does not have the PDB extension, the file is assumed to be a gray-scale image, where black represents regions that are topologically low, and white represents regions that are topologically high. The --min_x, --max_x, --min_y, and --max_y parameters are used to specify the region where the mesh should be generated. The --step_x and --step_y parameters define the x-y distance between adjacent points. The --max_height parameter determines the height of the bilayer model at those locations where the image is white; black regions are assigned a height of 0. This feature is only available if the python PIL module has been installed on your system. Example: --surface_filename mymesh.png",
                70,
                subsequent_indent="      ",
            )
        )
        print()
        print("Methods for resolving lipid clashes")
        print("===================================")
        print()
        print(
            textwrap.fill(
                "--delete_clashing_lipids: It's common for lipids to sterically clash at the interface of two adjacent surface-mesh tessellated triangles. If this parameter is set to TRUE, any clashing lipids are deleted. Example: --delete_clashing_lipids TRUE",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--clash_cutoff: If you do choose to delete clashing lipids, this parameter determines how close two atoms must be (in Angstroms) to constitute a steric clash. Example: --clash_cutoff 2.0",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--clashing_potential_margin: Lipid clashes occur at the edges of adjacent tessellated triangles. If these triangles are very large, it's faster to only check for clashes and holes near the triangle edges. This variable specifies how far from the edges, in Angstroms, that LipidWrapper should look for clashes and holes. Example: --clashing_potential_margin 25.0",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--fill_holes: Deleting lipids often leaves holes in the membrane. If this parameter is set to TRUE, LipidWrapper tries to fill the hole. Example: --fill_holes TRUE",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--very_distant_lipids_cutoff: LipidWrapper determines if two lipids clash by comparing the distance between every atom in the first lipid with every atom in the second lipid. This can be computationally expensive. However, sometimes two lipids are so distant from each other, that it's obvious there are no clashes, making the pair-wise comparison unnecessary. Before performing this expensive pair-wise comparison, LipidWrapper calculates the distance between one atom of each lipid. If this distance is greater than this user-specified cutoff, the program will simply assume there are no clashes. WARNING: Remember to consider the width of your lipid bilayer when choosing this value. Adjacent lipids on opposite sides of the bilayer can seem distant when considering the distance between their headgroups, for example. Example: --very_distant_lipids_cutoff 50.0",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--triangle_center_proximity_cutoff_distance: Lipid steric clashes/holes typically occur between lipids that belong to adjacent tessellated triangles. However, if tessellated triangles are small enough, clashes are possible between lipids that belong to non-adjacent triangles as well. Consequently, in addition to checking for adjacency, LipidWrapper also checks the distance between the triangle centers, using this user-specified value as a cutoff. Example: --triangle_center_proximity_cutoff_distance 50.0",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--fill_hole_exhaustiveness: Essentially, how long LipidWrapper should try to fill the holes. Example: --fill_hole_exhaustiveness 10",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--memory_optimization_factor: When the tessellated triangles are very large and consequently contain many individual lipids, the extensive pairwise distance comparisons required can result in memory errors. This parameter tells lipid Wrapper to divide the list of atoms being compared into smaller chunks. The pairwise distance comparison is performed piecewise on each chunk-chunk pair and so uses less memory, albeit at the expensive of speed. Only increase the value of this parameter if you run into memory errors. Example: --memory_optimization_factor 1",
                70,
                subsequent_indent="      ",
            )
        )
        print()
        print("Additional options")
        print("==================")
        print()
        print(
            textwrap.fill(
                "--number_of_processors: Using multiple processors can significantly increase the speed of the LipidWrapper algorithm. Example: --number_of_processors 8",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                '--show_grid_points: Aside from producing PDB coordinates for lipid atoms, additional coordinates will be appended to the bottom of the output containing "atoms" named "X" that specify the location of the surface mesh points. Example: --show_grid_points TRUE',
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                '--create_triangle_tcl_file: A separate file named "triangles.tcl" will be generated containing a tcl script that can be run in VMD to visualize the mesh surface. Example: --create_triangle_tcl_file TRUE',
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--output_directory: If an output directory is specified, all LipidWrapper output files, as well as additional files representing the intermediate steps required to build the final bilayer, will be saved in that directory. Example: --output_directory ./my_output/",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--use_disk_instead_of_memory: For very large systems, storing the growing model in memory can be problematic. If this parameter is set to TRUE, the growing model will be stored on the hard disk instead. However, expect longer execution times if this parameter is set to TRUE. Example: --use_disk_instead_of_memory TRUE",
                70,
                subsequent_indent="      ",
            )
        )
        print(
            textwrap.fill(
                "--compress_output: Depending on the user options selected, LipidWrapper output can require a lot of disk space. If this parameter is set to TRUE, the output will be automatically compressed using the gzip algorithm (Lempel-Ziv coding LZ77). The files can be uncompressed with the UNIX gunzip utility, or similar Windows-based packages. Example: --compress_output TRUE",
                70,
                subsequent_indent="      ",
            )
        )
        print()
        print("Example")
        print("=======")
        print()
        print(
            textwrap.fill(
                'python lipidwrapper.py --surface_equation "z = 250*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)" --min_x 500 --max_x 1000 --min_y 500 --max_y 1000 --step_x 25 --step_y 25 --lipid_pdb_filename lipid.pdb --lipid_headgroup_marker "_P,CHL1_O3" --delete_clashing_lipids TRUE --clash_cutoff 1.0 --fill_holes TRUE --fill_hole_exhaustiveness 10 > lipid_model.pdb',
                70,
                subsequent_indent="      ",
            )
        )
        print()
        sys.exit(0)

    # defaults
    params = {}
    params["surface_filename"] = (
        ""  # could be a PDB or image file, depending on surface_source value
    )
    params["surface_equation"] = (
        "z = 100*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)"  # used if surface_source is set to "EQUATION"
    )
    params["min_x"] = 500  # used if surface_source is PNG or EQUATION
    params["max_x"] = 750  # used if surface_source is PNG or EQUATION
    params["min_y"] = 500  # used if surface_source is PNG or EQUATION
    params["max_y"] = 750  # used if surface_source is PNG or EQUATION
    params["step_x"] = 30  # used if surface_source is PNG or EQUATION
    params["step_y"] = 30  # used if surface_source is PNG or EQUATION
    params["max_height"] = 0  # used if surface_source is PNG
    params["lipid_pdb_filename"] = (
        ""  # the filename containing the small, planar lipid model
    )
    params["lipid_headgroup_marker"] = (
        "_P,CHL1_O3"  # by default, any phosphate atom is considered a marker for the lipid headgroup, and also any O3 atom belonging to a cholesterol
    )
    params["show_grid_points"] = "FALSE"
    params["create_triangle_tcl_file"] = "FALSE"
    params["delete_clashing_lipids"] = "FALSE"
    params["use_disk_instead_of_memory"] = "FALSE"
    params["clash_cutoff"] = 2.0
    params["fill_holes"] = "TRUE"
    params["output_directory"] = ""
    params["fill_hole_exhaustiveness"] = 10
    params["number_of_processors"] = 1
    params["clashing_potential_margin"] = 25.0
    params["triangle_center_proximity_cutoff_distance"] = 50.0
    params["memory_optimization_factor"] = 1
    params["very_distant_lipids_cutoff"] = 50.0
    params["compress_output"] = "FALSE"

    # get commandline parameters
    options, remainder = getopt.getopt(
        argv[1:],
        "",
        [
            "surface_filename=",
            "surface_equation=",
            "min_x=",
            "max_x=",
            "min_y=",
            "max_y=",
            "step_x=",
            "step_y=",
            "max_height=",
            "lipid_pdb_filename=",
            "lipid_headgroup_marker=",
            "show_grid_points=",
            "create_triangle_tcl_file=",
            "delete_clashing_lipids=",
            "clash_cutoff=",
            "fill_holes=",
            "fill_hole_exhaustiveness=",
            "output_directory=",
            "number_of_processors=",
            "use_disk_instead_of_memory=",
            "clashing_potential_margin=",
            "triangle_center_proximity_cutoff_distance=",
            "memory_optimization_factor=",
            "very_distant_lipids_cutoff=",
            "compress_output=",
        ],
    )

    # set parameters to variables
    params_string = [
        "compress_output",
        "surface_filename",
        "surface_equation",
        "lipid_pdb_filename",
        "lipid_headgroup_marker",
        "show_grid_points",
        "create_triangle_tcl_file",
        "delete_clashing_lipids",
        "fill_holes",
        "output_directory",
        "use_disk_instead_of_memory",
    ]
    params_floats = [
        "very_distant_lipids_cutoff",
        "memory_optimization_factor",
        "triangle_center_proximity_cutoff_distance",
        "clashing_potential_margin",
        "min_x",
        "max_x",
        "min_y",
        "max_y",
        "step_x",
        "step_y",
        "max_height",
        "clash_cutoff",
        "fill_hole_exhaustiveness",
        "number_of_processors",
    ]

    for opt, arg in options:
        opt = opt.replace("-", "")
        if opt in params_floats:
            arg = float(arg)
        params[opt] = arg

    # some parameters should be integers
    params["fill_hole_exhaustiveness"] = int(params["fill_hole_exhaustiveness"])
    params["number_of_processors"] = int(params["number_of_processors"])
    params["memory_optimization_factor"] = int(params["memory_optimization_factor"])

    # directories should end in / or \ (depending on os)
    if params["output_directory"] != "" and params["output_directory"][-1:] != "/":
        params["output_directory"] = params["output_directory"] + "/"

    # check if running windows. If so, you can only use one processor
    if platform.system().lower() == "windows" and params["number_of_processors"] > 1:
        print(
            "REMARK WARNING: Use of multiple processors is only supported on Linux and OS X."
        )
        params["number_of_processors"] = 1

    # Print out header
    toprint = []
    toprint.append(
        "REMARK Parameters: (use the --help command-line parameter for further explanation)"
    )
    for param in list(params.keys()):
        toprint.append("REMARK      " + param + ": " + str(params[param]))
    toprint.append("")

    # create the output directory if necessary, and write the parameters used to a file
    if params["output_directory"] == "":
        print("\n".join(toprint))
    else:
        try:
            os.mkdir(params["output_directory"])
        except:
            pass

        f = open(params["output_directory"] + "parameters.input", "w")
        f.write("\n".join(toprint))
        f.close()

    # in the case of the lipid_headgroup_marker, split it by the comma.
    params["lipid_headgroup_marker"] = [
        t.strip() for t in params["lipid_headgroup_marker"].split(",")
    ]
    params["lipid_headgroup_marker"] = [
        (None, t.split("_")[0], None, t.split("_")[1])
        for t in params["lipid_headgroup_marker"]
    ]  # not that chain and resid are set to nothing, so any chain or resid will do

    # TRUE/FALSE answers need to be in caps.
    params["show_grid_points"] = params["show_grid_points"].upper()
    params["create_triangle_tcl_file"] = params["create_triangle_tcl_file"].upper()
    params["delete_clashing_lipids"] = params["delete_clashing_lipids"].upper()
    params["fill_holes"] = params["fill_holes"].upper()
    params["use_disk_instead_of_memory"] = params["use_disk_instead_of_memory"].upper()
    params["compress_output"] = params["compress_output"].upper()

    # specify the location of the temporary directory programatically
    params["memory_store_dir"] = params["output_directory"] + "store_in_memory.tmp/"

    # now check a few of the parameters to make sure they're valid
    if not os.path.exists(params["lipid_pdb_filename"]):
        print(
            (
                "ERROR: The file specified by the --lipid_pdb_filename parameter ("
                + params["lipid_pdb_filename"]
                + ") does not exist.\n"
            )
        )
        sys.exit(0)

    return params
