## imports

# standard
import gc

# custom
import numpy
import scipy.spatial.distance

# local
from . import multiprocessing_utils
from . import file_io


## methods


def remove_steric_clashes(molecules_by_triangle: list, params: dict):
    """Remove lipids that have steric clashes

    Arguments:
    molecules_by_triangle -- A list of tuples, where each tuple contains a Triangle object and a list of lipid molecules (Molecule objects) that belong to that triangle
    params -- A dictionary, the user-specified command-line parameters

    """

    class remove_clashes_multiprocessing(multiprocessing_utils.general_task):
        """A class for identifying lipids that have steric clashes"""

        def value_func(self, item, results_queue):
            """Identify lipids that have steric clashes

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output

            """

            triangle_index1 = item[0]
            triangle_index2 = item[1]
            triangle1 = item[2]
            triangle1_lipids = item[3]
            triangle2 = item[4]
            triangle2_lipids = item[5]
            params = item[6]

            self.print_star_if_appropriate(item[7])

            if triangle1.near_other_triangle(
                triangle2, params
            ):  # so only the lipids of proximate triangles are considered

                if params["use_disk_instead_of_memory"] == "TRUE":
                    triangle1_lipids = file_io.load_pickle(triangle1_lipids, params)
                if params["use_disk_instead_of_memory"] == "TRUE":
                    triangle2_lipids = file_io.load_pickle(triangle2_lipids, params)

                clash_map = {}

                # now generate a numpy array containing all the headgroups of the lipids of each triangle
                triangle1_headgroups = numpy.empty((len(triangle1_lipids), 3))
                for lipid_index, lipid in enumerate(triangle1_lipids):
                    headgroup_loc = lipid.all_atoms_numpy[
                        lipid.get_headgroup_index(params["lipid_headgroup_marker"])
                    ]
                    triangle1_headgroups[lipid_index][0] = headgroup_loc[0]
                    triangle1_headgroups[lipid_index][1] = headgroup_loc[1]
                    triangle1_headgroups[lipid_index][2] = headgroup_loc[2]

                triangle2_headgroups = numpy.empty((len(triangle2_lipids), 3))
                for lipid_index, lipid in enumerate(triangle2_lipids):
                    headgroup_loc = lipid.all_atoms_numpy[
                        lipid.get_headgroup_index(params["lipid_headgroup_marker"])
                    ]
                    triangle2_headgroups[lipid_index][0] = headgroup_loc[0]
                    triangle2_headgroups[lipid_index][1] = headgroup_loc[1]
                    triangle2_headgroups[lipid_index][2] = headgroup_loc[2]

                # get the indices of all the lipids in the margin of the first lipid
                indices_in_margin1 = []
                for idx, lip in enumerate(triangle1_lipids):
                    if lip.in_triangle_margin == True:
                        indices_in_margin1.append(idx)
                indices_in_margin1 = numpy.array(indices_in_margin1)

                if (
                    len(indices_in_margin1) > 0
                ):  # so there are some lipids in the margin

                    # get the indices of all the lipids in the margin of the second lipid
                    indices_in_margin2 = []
                    for idx, lip in enumerate(triangle2_lipids):
                        if lip.in_triangle_margin == True:
                            indices_in_margin2.append(idx)
                    indices_in_margin2 = numpy.array(indices_in_margin2)

                    if (
                        len(indices_in_margin2) > 0
                    ):  # so there are some lipids in the margin

                        # now, look at distances between all headgroups in margin to identify ones that are close enough to potentially clash
                        dists = scipy.spatial.distance.cdist(
                            triangle1_headgroups[indices_in_margin1],
                            triangle2_headgroups[indices_in_margin2],
                        )
                        dists = dists < params["very_distant_lipids_cutoff"]
                        indices_to_look_at1 = indices_in_margin1[
                            numpy.nonzero(dists)[0]
                        ]
                        indices_to_look_at2 = indices_in_margin2[
                            numpy.nonzero(dists)[1]
                        ]

                        # now do a pairwise comparison of all lipids, looking for clashes
                        for t in range(len(indices_to_look_at1)):
                            lipid_index1 = indices_to_look_at1[t]
                            lipid1 = triangle1_lipids[lipid_index1]
                            if lipid1.in_triangle_margin == True:
                                lipid_index2 = indices_to_look_at2[t]
                                lipid2 = triangle2_lipids[lipid_index2]
                                if lipid2.in_triangle_margin == True:
                                    if two_lipids_clash(
                                        lipid1,
                                        lipid2,
                                        params["clash_cutoff"],
                                        1,
                                        params,
                                        False,
                                    ):

                                        # there's a clash. update the clash map

                                        try:
                                            clash_map[
                                                (lipid_index1, triangle_index1)
                                            ].append((lipid_index2, triangle_index2))
                                        except:
                                            clash_map[
                                                (lipid_index1, triangle_index1)
                                            ] = []
                                            clash_map[
                                                (lipid_index1, triangle_index1)
                                            ].append((lipid_index2, triangle_index2))

                                        try:
                                            clash_map[
                                                (lipid_index2, triangle_index2)
                                            ].append((lipid_index1, triangle_index1))
                                        except:
                                            clash_map[
                                                (lipid_index2, triangle_index2)
                                            ] = []
                                            clash_map[
                                                (lipid_index2, triangle_index2)
                                            ].append((lipid_index1, triangle_index1))
                        self.results.append(clash_map)

    # generate a clash map, whcih specifies which lipids clash with each other
    some_input = []
    gc.disable()  # because appending complex objects to a list
    t = 0
    for triangle_index1 in range(len(molecules_by_triangle) - 1):
        for triangle_index2 in range(triangle_index1 + 1, len(molecules_by_triangle)):

            t = t + 1

            triangle1 = molecules_by_triangle[triangle_index1][0]
            triangle1_lipids = molecules_by_triangle[triangle_index1][1]

            triangle2 = molecules_by_triangle[triangle_index2][0]
            triangle2_lipids = molecules_by_triangle[triangle_index2][1]

            some_input.append(
                (
                    triangle_index1,
                    triangle_index2,
                    triangle1,
                    triangle1_lipids,
                    triangle2,
                    triangle2_lipids,
                    params,
                    t,
                )
            )
    gc.enable()

    tmp = multiprocessing_utils.multi_threading(
        some_input,
        params["number_of_processors"],
        remove_clashes_multiprocessing,
        params,
        "REMARK            (step 1) ",
    )

    # now combine all the clash maps into one
    clash_map = {}
    for amap in tmp.results:
        for akey in list(amap.keys()):
            try:
                clash_map[akey].extend(amap[akey])
            except:
                clash_map[akey] = []
                clash_map[akey].extend(amap[akey])

    # how we actually delete molecules
    # keep deleting molecules that clash until everything's resolved.
    # this part runs on one processor, but clash detection ran on multiple ones

    # start eliminating molecules until all clashes in the clashmap are resolved, starting with the molecule that has the most clashes
    gc.disable()
    lipids_to_delete = {}
    while len(clash_map) > 0:  # keep going until all clashes are resolved

        # identify the molecule that makes the most clashes
        most_clashes = 0
        most_clashes_mol_index = -1
        for mol_index in list(clash_map.keys()):
            num_clashes = len(clash_map[mol_index])
            if num_clashes > most_clashes:
                most_clashes = num_clashes
                most_clashes_mol_index = mol_index

        # now go through each of the ones it clashes with and remove it from their lists
        for other_clasher in clash_map[most_clashes_mol_index]:
            clash_map[other_clasher].remove(most_clashes_mol_index)
            if len(clash_map[other_clasher]) == 0:
                del clash_map[other_clasher]

        # now assign its value in the lipids lists to None, to be deleted later
        try:
            lipids_to_delete[most_clashes_mol_index[1]].append(
                most_clashes_mol_index[0]
            )
        except:
            lipids_to_delete[most_clashes_mol_index[1]] = []
            lipids_to_delete[most_clashes_mol_index[1]].append(
                most_clashes_mol_index[0]
            )

        # now remove this one from clash_map as well
        del clash_map[most_clashes_mol_index]
    gc.enable()

    # now delete the lipids that have been marked for deletion by being assigned a value of None
    class remove_clashes2_multiprocessing(multiprocessing_utils.general_task):
        """A class for removing lipids that have steric clashes"""

        def value_func(self, item, results_queue):
            """Remove lipids that have steric clashes

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output

            """

            triangle_index = item[0]
            lipid_indices_to_delete = item[1]
            lipids_list = item[2]
            params = item[4]

            self.print_star_if_appropriate(item[3])

            if params["use_disk_instead_of_memory"] == "TRUE":
                lipids = file_io.load_pickle(lipids_list, params)
            else:
                lipids = lipids_list

            for lipid_index in lipid_indices_to_delete:
                lipids[lipid_index] = None
            while None in lipids:
                lipids.remove(None)

            if params["use_disk_instead_of_memory"] == "TRUE":
                file_io.save_pickle(lipids, params, lipids_list)
            else:
                gc.disable()
                self.results.append((triangle_index, lipids))
                gc.enable()

    # now actually go through and remove all the "lipids" that have been asigned None
    some_input = []
    gc.disable()
    t = 0
    for triangle_index in list(lipids_to_delete.keys()):
        t = t + 1
        lipid_indices_to_delete = lipids_to_delete[triangle_index]
        some_input.append(
            (
                triangle_index,
                lipid_indices_to_delete,
                molecules_by_triangle[triangle_index][1],
                t,
                params,
            )
        )
    gc.enable()

    tmp = multiprocessing_utils.multi_threading(
        some_input,
        params["number_of_processors"],
        remove_clashes2_multiprocessing,
        params,
        "REMARK            (step 2) ",
    )

    # update molecules_by_triangle
    if (
        params["use_disk_instead_of_memory"] != "TRUE"
    ):  # so you need to reconstruct molecules_by_triangle
        for triangle_id, lipids_list in tmp.results:
            if len(molecules_by_triangle[triangle_id]) == 2:
                molecules_by_triangle[triangle_id] = (
                    molecules_by_triangle[triangle_id][0],
                    lipids_list,
                )
            else:
                molecules_by_triangle[triangle_id] = (
                    molecules_by_triangle[triangle_id][0],
                    lipids_list,
                    molecules_by_triangle[triangle_id][2],
                )


def two_lipids_clash(
    mol1,
    mol2,
    cutoff: float,
    num_sub_partitions: int,
    params: dict,
    very_large_distance_check: bool = True,
):
    """Determine whether two lipid molecules clash

    Arguments:
    mol1 -- A Molecule object, the first lipid
    mol2 -- A Molecule object, the second lipid
    cutoff -- A float, how close the two lipids must be to constitute a "clash"
    num_sub_partitions -- An integer, the number of partitions into which the atoms of each lipid molecule are divided. Clashes are then determined pairwise on each partition, rather than comparing every atom of one lipid to every atom of the other. Important for large system to avoid memory problems, but it can usually just be set to 1.
    params -- A dictionary, the user-specified comand-line parameters
    very_large_distance_check -- An optional Boolean, to enable the option of eliminating steric clashes early by examining the distance between the first atoms of each lipid

    Returns:
    A Boolean, True if the two lipids clash, False otherwise.

    """

    # of the user provided a Molecule object, use just the coordinate numpy array
    if not type(mol1) is numpy.ndarray:
        mol1 = mol1.all_atoms_numpy
    if not type(mol2) is numpy.ndarray:
        mol2 = mol2.all_atoms_numpy

    # first, check if the bounding boxes of each lipid couldn't possiblely overlap
    margin = numpy.array([cutoff, cutoff, cutoff])
    mol1_min = numpy.min(mol1, 0) - margin
    mol2_max = numpy.max(mol2, 0) + margin

    if mol1_min[0] > mol2_max[0]:
        return False
    if mol1_min[1] > mol2_max[1]:
        return False
    if mol1_min[2] > mol2_max[2]:
        return False

    mol1_max = numpy.max(mol1, 0) + margin
    mol2_min = numpy.min(mol2, 0) - margin

    if mol2_min[0] > mol1_max[0]:
        return False
    if mol2_min[1] > mol1_max[1]:
        return False
    if mol2_min[2] > mol1_max[2]:
        return False

    # now check if the distance between their first atoms is so great, they are unlikely to overlap
    if (
        very_large_distance_check == True
        and numpy.linalg.norm(mol1[0] - mol2[0]) > params["very_distant_lipids_cutoff"]
    ):
        return False

    # now reduce the points to the ones that really need to be compared. Basically, stripping the points
    # I think this just retains points in the overlapping regions of the bounding boxes
    reduced_boundary_max = numpy.min(numpy.vstack((mol1_max, mol2_max)), 0)
    reduced_boundary_min = numpy.max(numpy.vstack((mol1_min, mol2_min)), 0)

    mol1_to_remove = mol1[:, 0] < reduced_boundary_max[0]
    if not True in mol1_to_remove:
        return False
    mol1_to_remove = numpy.logical_and(
        mol1_to_remove, (mol1[:, 0] > reduced_boundary_min[0])
    )
    if not True in mol1_to_remove:
        return False
    mol1_to_remove = numpy.logical_and(
        mol1_to_remove, (mol1[:, 1] < reduced_boundary_max[1])
    )
    if not True in mol1_to_remove:
        return False
    mol1_to_remove = numpy.logical_and(
        mol1_to_remove, (mol1[:, 1] > reduced_boundary_min[1])
    )
    if not True in mol1_to_remove:
        return False
    mol1_to_remove = numpy.logical_and(
        mol1_to_remove, (mol1[:, 2] < reduced_boundary_max[2])
    )
    if not True in mol1_to_remove:
        return False
    mol1_to_remove = numpy.logical_and(
        mol1_to_remove, (mol1[:, 2] > reduced_boundary_min[2])
    )
    if not True in mol1_to_remove:
        return False

    mol2_to_remove = mol2[:, 0] < reduced_boundary_max[0]
    if not True in mol2_to_remove:
        return False
    mol2_to_remove = numpy.logical_and(
        mol2_to_remove, (mol2[:, 0] > reduced_boundary_min[0])
    )
    if not True in mol2_to_remove:
        return False
    mol2_to_remove = numpy.logical_and(
        mol2_to_remove, (mol2[:, 1] < reduced_boundary_max[1])
    )
    if not True in mol2_to_remove:
        return False
    mol2_to_remove = numpy.logical_and(
        mol2_to_remove, (mol2[:, 1] > reduced_boundary_min[1])
    )
    if not True in mol2_to_remove:
        return False
    mol2_to_remove = numpy.logical_and(
        mol2_to_remove, (mol2[:, 2] < reduced_boundary_max[2])
    )
    if not True in mol2_to_remove:
        return False
    mol2_to_remove = numpy.logical_and(
        mol2_to_remove, (mol2[:, 2] > reduced_boundary_min[2])
    )
    if not True in mol2_to_remove:
        return False

    mol1_reduced = mol1[numpy.nonzero(mol1_to_remove)[0]]
    mol2_reduced = mol2[numpy.nonzero(mol2_to_remove)[0]]

    # now do a pairwise distance comparison between the remaining atoms
    if num_sub_partitions == 1:
        if True in (scipy.spatial.distance.cdist(mol1_reduced, mol2_reduced) < cutoff):
            return True
    else:
        a1s = numpy.array_split(
            mol1_reduced, num_sub_partitions
        )  # in my benchmarks, 35 gave good results
        a2s = numpy.array_split(mol2_reduced, num_sub_partitions)

        for mol1_reduced in a1s:
            for mol2_reduced in a2s:
                if len(mol1_reduced) > 0 and len(mol2_reduced) > 0:
                    if True in (
                        scipy.spatial.distance.cdist(mol1_reduced, mol2_reduced)
                        < cutoff
                    ):
                        return True

    return False


def indices_of_close_pts(points1, points2, cutoff: float, num_sub_partitions: int):
    """Examine two sets of points and return the indices of the points that are close to each other

    Arguments:
    points1 -- A nx3 numpy array, a set of points to be considered
    points2 -- A nx3 numpy array, a second set of points to be considered
    cutoff -- A float, how close the two lipids must be to constitute a "clash"
    num_sub_partitions -- An integer, the number of partitions into which the atoms of each lipid molecule are divided. Clashes are then determined pairwise on each partition, rather than comparing every atom of one lipid to every atom of the other. Important for large system to avoid memory problems, but it can usually just be set to 1.

    Returns:
    A tuple, containing a numpy array with indices from the first molecule and a numpy array with indices from the second molecule

    """

    if num_sub_partitions == 1:
        dists = (
            scipy.spatial.distance.cdist(points1, points2) < cutoff
        )  # which ones clash
        return numpy.nonzero(dists)  # these are indices that clash
    else:
        a1s = numpy.array_split(
            points1, num_sub_partitions
        )  # in my benchmarks, 35 gave good results
        a2s = numpy.array_split(points2, num_sub_partitions)

        points1_indices = []
        points2_indices = []

        a1s_index = 0
        for points1 in a1s:
            a2s_index = 0
            for points2 in a2s:
                if len(points1) > 0 and len(points2) > 0:
                    dists = scipy.spatial.distance.cdist(points1, points2) < cutoff
                    indices = numpy.nonzero(dists)
                    points1_indices.extend(indices[0] + a1s_index)
                    points2_indices.extend(indices[1] + a2s_index)
                a2s_index = a2s_index + len(points2)
            a1s_index = a1s_index + len(points1)

        points1_indices = numpy.array([points1_indices])
        points2_indices = numpy.array([points2_indices])

        return (points1_indices[0], points2_indices[0])
