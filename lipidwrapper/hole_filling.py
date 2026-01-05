## imports

# standard
import gc
import random
import typing

# custom
import numpy
import scipy.spatial.distance

# local
from . import multiprocessing_utils
from . import file_io
from . import clash_removal
from . import lipid_positioning


## methods


def fill_in_lipid_holes(molecules_by_triangle: list, params: dict) -> list:
    headgroup_locs = []

    index_to_try = 0
    while (
        len(headgroup_locs) < 5
    ):  # I feel like there needs to be at least 5 to get any kind of representative sample. Because the first triangel could conceivably not have any headgroups in it or only have one, keep looking until you find one that does
        # I'm not going to search through all of them to find the max because that might require loading a lot of pickles
        headgroup_locs = []

        # the first step is to get the average and minimum distance between headgroups.
        if params["use_disk_instead_of_memory"] == "TRUE":
            lipids = file_io.load_pickle(molecules_by_triangle[0][1], params)
        else:
            lipids = molecules_by_triangle[index_to_try][
                1
            ]  # Just look one of the triangles as representative

        gc.disable()
        for alipid in lipids:
            headgroup_locs.append(
                alipid.all_atoms_numpy[
                    alipid.get_headgroup_index(params["lipid_headgroup_marker"])
                ]
            )
        gc.enable()

        index_to_try = index_to_try + 1

    headgroup_locs = numpy.vstack(headgroup_locs)
    headgroup_dists = scipy.spatial.distance.squareform(
        scipy.spatial.distance.pdist(headgroup_locs)
    )
    headgroup_min_dists = numpy.empty(len(headgroup_dists))

    for indx in range(len(headgroup_dists)):
        t = headgroup_dists[indx]
        headgroup_min_dists[indx] = numpy.min(t[numpy.nonzero(t)])

    average_dist_between_headgroups = int(
        numpy.round(numpy.average(headgroup_min_dists))
    )
    min_dist_between_headgroups = numpy.min(headgroup_min_dists)
    pt_step = max(
        [1, int(average_dist_between_headgroups / 3.0)]
    )  # so grid points every third of the way between headgroups

    # now, determine which triangles are adjacent
    adjacent_triangles_map = {}
    gc.disable()
    for index1 in range(len(molecules_by_triangle) - 1):
        triangle_pts1 = molecules_by_triangle[index1][0]
        for index2 in range(index1 + 1, len(molecules_by_triangle)):
            triangle_pts2 = molecules_by_triangle[index2][0]

            if triangle_pts1.near_other_triangle(triangle_pts2, params):

                try:
                    adjacent_triangles_map[index1].append(index2)
                except:
                    adjacent_triangles_map[index1] = []
                    adjacent_triangles_map[index1].append(index2)

                try:
                    adjacent_triangles_map[index2].append(index1)
                except:
                    adjacent_triangles_map[index2] = []
                    adjacent_triangles_map[index2].append(index1)

    class lipid_inserts_multiprocessing(multiprocessing_utils.general_task):
        def value_func(
            self, item: tuple, results_queue: typing.Optional[typing.Any]
        ) -> None:
            molecules_by_triangle_index = item[0]
            triangle_pts = item[1]
            lipids = item[2]
            adjacent_lipids = item[
                3
            ]  # molecules of lipids in neighboring triangles, but NOT in this one (i.e., a triangle is not adjacent to itself)
            params = item[4]
            average_dist_between_headgroups = item[5]
            min_dist_between_headgroups = item[6]
            pt_step = item[7]

            self.print_star_if_appropriate(molecules_by_triangle_index)

            if params["use_disk_instead_of_memory"] == "TRUE":
                lipids = file_io.load_pickle(lipids, params)

            # now get the plane going between these three points

            # the order of the triangle points could potentially matter in the case of
            # right triangles. So we potentially need to make sure every order is considered,
            # though we can abort early if an acceptable solution is found.
            # basically, in the case of right triangles, the point opposite the hypotenuse
            # needs to be projected onto the hypotenuse. With other kinds of triangles,
            # it can really be any point projected onto the opposite side.

            combos = []
            combos.append((triangle_pts[0], triangle_pts[1], triangle_pts[2]))
            combos.append((triangle_pts[0], triangle_pts[2], triangle_pts[1]))
            combos.append((triangle_pts[1], triangle_pts[2], triangle_pts[0]))

            for combo in combos:

                pt1 = combo[0]
                pt2 = combo[1]
                pt3 = combo[2]

                # project pt3 onto the line segment pt1-pt2
                u = pt1 - pt2
                v = pt1 - pt3
                u = u / numpy.linalg.norm(u)
                new_pt = pt1 - numpy.dot(u, v) * u  # this is the projected point

                # make sure the project point isn't equal to one of the triangle verticies
                if not numpy.array_equal(pt3, new_pt) and not numpy.array_equal(
                    pt1, new_pt
                ):
                    break

            vec1 = pt3 - new_pt
            vec2 = pt1 - new_pt

            vec1 = vec1 / numpy.linalg.norm(
                vec1
            )  # two perpenticular vectors in the plane
            vec2 = vec2 / numpy.linalg.norm(vec2)  # and a point in the plane

            plane_normal = numpy.cross(vec1, vec2)  # a normal to the plane
            plane_normal = plane_normal / numpy.linalg.norm(plane_normal)

            # good to get a scalar equation for the plane too: ax + by + cz + d = 0
            scalar_eq_a = plane_normal[0]
            scalar_eq_b = plane_normal[1]
            scalar_eq_c = plane_normal[2]
            scalar_eq_d = -numpy.dot(triangle_pts.center(), plane_normal)

            # now that the plane has been identified, find the average distance between the plane and lipid headgroups
            # also, start adding lipids that could clash with future inserted lipids into the neighborhood_lipids_that_could_clash list. All lipids in the margin and submargin of the central triangle will be added.
            lipid_head_indices = numpy.empty(len(lipids), dtype=int)
            for indx, lipid in enumerate(lipids):
                lipid_head_indices[indx] = lipid.get_headgroup_index(
                    params["lipid_headgroup_marker"]
                )

            all_lipid_heads_loc_in_central_triangle = numpy.empty(
                (len(lipid_head_indices), 3)
            )
            neighborhood_lipids_that_could_clash = []
            headgroup_locs_of_lipids_that_could_clash = []
            for t in range(len(lipid_head_indices)):
                all_lipid_heads_loc_in_central_triangle[t] = lipids[t].all_atoms_numpy[
                    lipid_head_indices[t]
                ]

                if (
                    lipids[t].in_triangle_margin == True
                    or lipids[t].in_triangle_submargin == True
                ):
                    neighborhood_lipids_that_could_clash.append(lipids[t])
                    headgroup_locs_of_lipids_that_could_clash.append(
                        lipids[t].all_atoms_numpy[lipid_head_indices[t]]
                    )

            three_scalars = numpy.array([scalar_eq_a, scalar_eq_b, scalar_eq_c])
            dists2 = numpy.empty(len(all_lipid_heads_loc_in_central_triangle))
            for indx in range(len(all_lipid_heads_loc_in_central_triangle)):
                lipid_head_pt = all_lipid_heads_loc_in_central_triangle[indx]
                dist = numpy.fabs(
                    numpy.dot(three_scalars, lipid_head_pt) + scalar_eq_d
                ) / numpy.power(numpy.dot(three_scalars, three_scalars), 0.5)
                dists2[indx] = dist

            if (
                len(dists2) == 0
            ):  # if there are no lipid headgroups in this triangle, so you can't proceed
                positioned_molecules = []
                if params["use_disk_instead_of_memory"] == "TRUE":
                    self.results.append(
                        (
                            molecules_by_triangle_index,
                            file_io.save_pickle(positioned_molecules, params),
                        )
                    )  # here save the results for later compilation
                else:
                    self.results.append(
                        (molecules_by_triangle_index, positioned_molecules)
                    )  # here save the results for later compilation
                return

            average_headgroup_dist_to_plane = numpy.mean(dists2)

            # Find the locations of the in-margin headgroups of all adjacent triangles
            gc.disable()
            for (
                lipids2
            ) in (
                adjacent_lipids
            ):  # note that this does NOT include the lipids in the central triangle, which were identified above.
                if params["use_disk_instead_of_memory"] == "TRUE":
                    lipids2 = file_io.load_pickle(lipids2, params)

                for alipid in lipids2:

                    if (
                        alipid.in_triangle_margin == True
                    ):  # so for neighboring triangles, we only care about the lipids that are in the margin, which might clash with future inserted lipids
                        neighborhood_lipids_that_could_clash.append(alipid)
                        headgroup_locs_of_lipids_that_could_clash.append(
                            alipid.all_atoms_numpy[
                                alipid.get_headgroup_index(
                                    params["lipid_headgroup_marker"]
                                )
                            ]
                        )
            gc.enable()

            # need to numpify headgroup_locs_of_lipids_that_could_clash
            headgroup_locs_of_lipids_that_could_clash = numpy.array(
                headgroup_locs_of_lipids_that_could_clash
            )

            # now flood the surface of both bilayers with points
            # first, generate a field of points
            s = numpy.arange(
                -triangle_pts.max_radius(),
                triangle_pts.max_radius(),
                pt_step,
                dtype=int,
            )
            t = numpy.arange(
                -triangle_pts.max_radius(),
                triangle_pts.max_radius(),
                pt_step,
                dtype=int,
            )
            pts = numpy.empty((len(s) * len(t), 3))
            for s_index, s_val in enumerate(s):
                for t_index, t_val in enumerate(t):
                    pt = s_val * vec1 + t_val * vec2
                    pts[s_index * len(t) + t_index][0] = pt[0]
                    pts[s_index * len(t) + t_index][1] = pt[1]
                    pts[s_index * len(t) + t_index][2] = pt[2]
            pts = numpy.array(pts) + triangle_pts.center()

            # check which of these points are within the central triangle
            indices_of_pts_in_triangle = (
                triangle_pts.get_indices_of_points_within_triangle_boundaries(pts)
            )
            pts_in_triangle = (
                pts[indices_of_pts_in_triangle]
                if len(indices_of_pts_in_triangle) > 0
                else numpy.array([])
            )

            # now remove points that are too far in the interior. Fill only at triangle edges
            smaller_tri = triangle_pts.new_triangle_expanded_by_margin(
                -params["clashing_potential_margin"]
            )
            indices_of_pts_in_triangle = (
                smaller_tri.get_indices_of_points_within_triangle_boundaries(
                    pts_in_triangle
                )
            )
            pts_in_triangle = numpy.delete(
                pts_in_triangle, indices_of_pts_in_triangle, axis=0
            )

            # create points above and below each of these grid points
            local_pts_to_examine = numpy.empty((2 * len(pts_in_triangle), 3))
            for apt_index, apt in enumerate(pts_in_triangle):
                # now get the two points above and below the plane
                starting_pts = numpy.array(
                    [
                        apt - plane_normal * average_headgroup_dist_to_plane,
                        apt + plane_normal * average_headgroup_dist_to_plane,
                    ]
                )

                # place those two points into the local_pts_to_examine numpy array
                for starting_pt_index, starting_pt in enumerate(starting_pts):
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][0] = (
                        starting_pt[0]
                    )
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][1] = (
                        starting_pt[1]
                    )
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][2] = (
                        starting_pt[2]
                    )

            # remove all pts that are too close to the headgroups
            indices_of_clashing_pts = clash_removal.indices_of_close_pts(
                headgroup_locs_of_lipids_that_could_clash,
                local_pts_to_examine,
                average_dist_between_headgroups,
                params["memory_optimization_factor"],
            )[1]
            local_pts_to_examine = numpy.delete(
                local_pts_to_examine, indices_of_clashing_pts, 0
            )

            # remove all remaining pts that clash with other lipid atoms (headgroups first to reduce number of pair-wise distance comparisons)
            for lip in neighborhood_lipids_that_could_clash:
                indices_of_clashing_pts = clash_removal.indices_of_close_pts(
                    lip.all_atoms_numpy,
                    local_pts_to_examine,
                    min_dist_between_headgroups,
                    params["memory_optimization_factor"],
                )[1]

                # indices_of_clashing_pts could be empty, so just try
                try:
                    local_pts_to_examine = numpy.delete(
                        local_pts_to_examine, indices_of_clashing_pts, 0
                    )
                except:
                    pass

            # now position lipids
            positioned_molecules = (
                []
            )  # can't know size, so can't preallocate in numpy array
            positioned_molecules_headgroup_locs = (
                []
            )  # can't know size, so can't preallocate
            gc.disable()  # because appending complex objects to a list

            for t in range(params["fill_hole_exhaustiveness"]):
                indxs = list(range(len(local_pts_to_examine)))
                random.shuffle(indxs)  # so not examining points sequentially
                for headgroup_loc_index in indxs:
                    if headgroup_loc_index < len(
                        local_pts_to_examine
                    ):  # the point could have been deleted, in which case you should skip

                        new_head_group_loc = local_pts_to_examine[headgroup_loc_index]

                        # determine the directionality of the lipid (i.e., points "up" or "down")
                        candidates_pts = numpy.array(
                            [
                                new_head_group_loc - plane_normal,
                                new_head_group_loc + plane_normal,
                            ]
                        )
                        dists_to_center = scipy.spatial.distance.cdist(
                            candidates_pts, numpy.array([triangle_pts.center()])
                        )

                        if dists_to_center[0] < dists_to_center[1]:
                            directionality = 1
                        else:
                            directionality = -1

                        # pick a random lipid
                        lipid = random.choice(lipids)  # maybe needs to be a copy?
                        lipid_head_loc_index = lipid.get_headgroup_index(
                            params["lipid_headgroup_marker"]
                        )
                        lipid_head_loc = lipid.all_atoms_numpy[lipid_head_loc_index]
                        lipid_center_loc = numpy.mean(lipid.all_atoms_numpy, 0)
                        lipid_length = numpy.linalg.norm(
                            lipid_head_loc - lipid_center_loc
                        )

                        # you should be working with a copy of the lipid, not the original
                        lipid = lipid.copy_of()
                        lipid.in_triangle_margin = True
                        lipid.in_triangle_submargin = False

                        # get new guide (static) template. This specifies where the lipid will ultimately be moved to
                        lipid_center_guidepoint = (
                            new_head_group_loc
                            - directionality * lipid_length * plane_normal
                        )
                        guide_static_template = numpy.array(
                            [new_head_group_loc, lipid_center_guidepoint]
                        )

                        # get new dynamic template. this is the starting location of the lipid before it's moved to the new location.
                        dynamic_template = numpy.array(
                            [lipid_head_loc, lipid_center_loc]
                        )

                        # get origin template. This is a destination at the origin. You'll move it here for rotating before moving it to the new location
                        origin = numpy.array([0.0, 0.0, 0.0])
                        origin2 = origin - lipid_length * numpy.array([0.0, 0.0, 1.0])
                        guide_origin_template = numpy.array([origin, origin2])

                        # move lipid to origin.
                        transform_data = lipid_positioning.get_transformation_data(
                            guide_origin_template, dynamic_template
                        )
                        lipid_positioning.apply_transformation(lipid, transform_data)

                        # now rotate about z axis
                        theta = random.random() * numpy.pi * 2.0
                        rot_max = numpy.array(
                            [
                                [numpy.cos(theta), -numpy.sin(theta), 0.0],
                                [numpy.sin(theta), numpy.cos(theta), 0.0],
                                [0.0, 0.0, 1.0],
                            ]
                        )
                        lipid.all_atoms_numpy = numpy.dot(
                            lipid.all_atoms_numpy, rot_max
                        )

                        # now move to correct location in bilayer
                        (
                            center_dynamic_pdb,
                            rot_quat,
                            center_static_pdb,
                        ) = lipid_positioning.get_transformation_data(
                            guide_static_template, guide_origin_template
                        )
                        lipid.all_atoms_numpy = (
                            lipid.all_atoms_numpy - center_dynamic_pdb
                        )
                        lipid.rotate_mol_quat(rot_quat)
                        lipid.all_atoms_numpy = (
                            lipid.all_atoms_numpy + center_static_pdb
                        )

                        # redefine the lead group location now that things have been moved
                        lipid_head_loc = lipid.all_atoms_numpy[lipid_head_loc_index]

                        # check to see if the positioned lipid clashes with other lipids
                        some_clash = False
                        first_pt_dists = scipy.spatial.distance.cdist(
                            headgroup_locs_of_lipids_that_could_clash,
                            numpy.array([lipid_head_loc]),
                        )
                        first_pt_close_indices = numpy.nonzero(
                            first_pt_dists < params["very_distant_lipids_cutoff"]
                        )[0]
                        for indx in first_pt_close_indices:
                            if (
                                clash_removal.two_lipids_clash(
                                    lipid,
                                    neighborhood_lipids_that_could_clash[indx],
                                    params["clash_cutoff"],
                                    1,
                                    params,
                                    False,
                                )
                                == True
                            ):
                                some_clash = True
                                break

                        if some_clash == False:

                            if len(positioned_molecules_headgroup_locs) > 0:

                                positioned_pt_dists = scipy.spatial.distance.cdist(
                                    numpy.array(positioned_molecules_headgroup_locs),
                                    numpy.array([lipid_head_loc]),
                                )
                                positioned_pt_close_indices = numpy.nonzero(
                                    positioned_pt_dists
                                    < params["very_distant_lipids_cutoff"]
                                )[0]
                                for indx in positioned_pt_close_indices:
                                    if (
                                        clash_removal.two_lipids_clash(
                                            lipid,
                                            positioned_molecules[indx],
                                            params["clash_cutoff"],
                                            1,
                                            params,
                                            False,
                                        )
                                        == True
                                    ):
                                        some_clash = True
                                        break

                        if some_clash == False:  # so it doesn't clash. save it.

                            positioned_molecules.append(lipid)  # remember, a copy
                            positioned_molecules_headgroup_locs.append(lipid_head_loc)

                            # now remove surface points from local_pts_to_examine that come close to the newly positioned lipid
                            dists = (
                                scipy.spatial.distance.cdist(
                                    lipid.all_atoms_numpy, local_pts_to_examine
                                )
                                < min_dist_between_headgroups
                            )  # which ones clash
                            indices_of_clashing_pts = numpy.nonzero(dists)[
                                1
                            ]  # these are indices that clash
                            local_pts_to_examine = numpy.delete(
                                local_pts_to_examine, indices_of_clashing_pts, 0
                            )

            # now add all these positioned lipids to the molecules_by_triangle list
            if params["use_disk_instead_of_memory"] == "TRUE":
                self.results.append(
                    (
                        molecules_by_triangle_index,
                        file_io.save_pickle(positioned_molecules, params),
                    )
                )  # here save the results for later compilation
            else:
                self.results.append(
                    (molecules_by_triangle_index, positioned_molecules)
                )  # here save the results for later compilation
            gc.enable()

    # fill the lipid holes using multiple processors if possible
    some_input = []
    for molecules_by_triangle_index in range(len(molecules_by_triangle)):
        triangle_pts = molecules_by_triangle[molecules_by_triangle_index][0]
        lipids = molecules_by_triangle[molecules_by_triangle_index][1]
        adjacent_lipids = [
            molecules_by_triangle[index][1]
            for index in adjacent_triangles_map[molecules_by_triangle_index]
        ]
        some_input.append(
            (
                molecules_by_triangle_index,
                triangle_pts,
                lipids,
                adjacent_lipids,
                params,
                average_dist_between_headgroups,
                min_dist_between_headgroups,
                pt_step,
            )
        )

    gc.enable()

    positioned_lipids = multiprocessing_utils.multi_threading(
        some_input,
        params["number_of_processors"],
        lipid_inserts_multiprocessing,
        params,
        "REMARK ",
    ).results

    # now organize the positioned_lipids into the same organization as molecules_by_triangle for subsequent processing
    positioned_lipids_by_triangle = []
    gc.disable()
    for molecules_by_triangle_index, positioned_molecules in positioned_lipids:
        positioned_lipids_by_triangle.append(
            (
                molecules_by_triangle[molecules_by_triangle_index][0],
                positioned_molecules,
                molecules_by_triangle_index,
            )
        )
    gc.enable()

    return positioned_lipids_by_triangle
