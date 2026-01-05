## imports

# standard
import os
import pathlib
import tempfile

# custom
import pytest

# local
import lipidwrapper


## constants

EXAMPLES_DIR = pathlib.Path(__file__).parent.parent / "examples"
FILES_DIR = EXAMPLES_DIR / "files"
LIPID_PDB = FILES_DIR / "lipid_example.pdb"


## fixtures


@pytest.fixture
def output_directory():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir + "/"


@pytest.fixture
def base_params(output_directory):
    return {
        "lipid_pdb_filename": str(LIPID_PDB),
        "lipid_headgroup_marker": [(None, "", None, "P"), (None, "CHL1", None, "O3")],
        "show_grid_points": "FALSE",
        "create_triangle_tcl_file": "FALSE",
        "delete_clashing_lipids": "FALSE",
        "use_disk_instead_of_memory": "FALSE",
        "clash_cutoff": 2.0,
        "fill_holes": "FALSE",
        "output_directory": output_directory,
        "fill_hole_exhaustiveness": 10,
        "number_of_processors": 1,
        "clashing_potential_margin": 15.0,
        "triangle_center_proximity_cutoff_distance": 50.0,
        "memory_optimization_factor": 1,
        "very_distant_lipids_cutoff": 50.0,
        "compress_output": "FALSE",
        "memory_store_dir": output_directory + "store_in_memory.tmp/",
        "surface_filename": "",
        "surface_equation": "z = 0",
        "min_x": 0,
        "max_x": 50,
        "min_y": 0,
        "max_y": 50,
        "step_x": 25,
        "step_y": 25,
        "max_height": 0,
    }


## classes


class TestSurfaceFromEquation:
    def test_flat_surface_generates_output(self, base_params):
        base_params["surface_equation"] = "z = 0"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0

    def test_curved_surface_generates_output(self, base_params):
        base_params["surface_equation"] = "z = 25.0*x*x/10000.0 + 25.0*y*y/10000.0"
        base_params["min_x"] = -50
        base_params["max_x"] = 50
        base_params["min_y"] = -50
        base_params["max_y"] = 50
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0

    def test_with_clash_removal(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["delete_clashing_lipids"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        assert any("remove_lipids_with_clashes" in f for f in output_files)

    def test_with_hole_filling(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["delete_clashing_lipids"] = "TRUE"
        base_params["fill_holes"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        assert any("lipid_holes_plugged" in f for f in output_files)

    def test_with_grid_points(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["show_grid_points"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        grid_points_file = base_params["output_directory"] + "grid_points.pdb"
        assert os.path.exists(grid_points_file)

    def test_with_tcl_file(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["create_triangle_tcl_file"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        tcl_file = base_params["output_directory"] + "triangles.tcl"
        assert os.path.exists(tcl_file)


class TestSurfaceFromImage:
    def test_image_surface_generates_output(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "face.png")
        base_params["max_height"] = 50
        base_params["min_x"] = -25
        base_params["max_x"] = 25
        base_params["min_y"] = -25
        base_params["max_y"] = 25
        base_params["step_x"] = 25
        base_params["step_y"] = 25
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0

    @pytest.mark.slow
    def test_image_surface_with_clash_removal(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "face.png")
        base_params["max_height"] = 50
        base_params["min_x"] = -25
        base_params["max_x"] = 25
        base_params["min_y"] = -25
        base_params["max_y"] = 25
        base_params["step_x"] = 25
        base_params["step_y"] = 25
        base_params["delete_clashing_lipids"] = "TRUE"
        base_params["fill_holes"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        assert any("lipid_holes_plugged" in f for f in output_files)


@pytest.mark.slow
class TestSurfaceFromPDBPoints:
    def test_pdb_surface_generates_output(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "influenza_mesh.pdb")
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0

    def test_pdb_surface_with_clash_removal(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "influenza_mesh.pdb")
        base_params["delete_clashing_lipids"] = "TRUE"
        base_params["fill_holes"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        assert any("lipid_holes_plugged" in f for f in output_files)


@pytest.mark.slow
class TestSurfaceFromDAEFile:
    def test_dae_surface_generates_output(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "concave_model.dae")
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0

    def test_dae_surface_with_clash_removal(self, base_params):
        base_params["surface_filename"] = str(FILES_DIR / "concave_model.dae")
        base_params["delete_clashing_lipids"] = "TRUE"
        base_params["fill_holes"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        assert any("lipid_holes_plugged" in f for f in output_files)


@pytest.mark.slow
class TestDiskStorage:
    def test_disk_storage_mode(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["use_disk_instead_of_memory"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0
        assert not os.path.exists(base_params["memory_store_dir"])


@pytest.mark.slow
class TestMultiprocessing:
    def test_multiple_processors(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["number_of_processors"] = 2
        lipidwrapper.run_with_params(base_params)
        output_files = os.listdir(base_params["output_directory"])
        pdb_files = [f for f in output_files if f.endswith(".pdb")]
        assert len(pdb_files) > 0


class TestOutputValidation:
    def test_full_bilayer_pdb_contains_atoms(self, base_params):
        base_params["surface_equation"] = "z = 0"
        lipidwrapper.run_with_params(base_params)
        full_bilayer_files = [
            f
            for f in os.listdir(base_params["output_directory"])
            if "full_bilayer" in f and f.endswith(".pdb")
        ]
        assert len(full_bilayer_files) == 1
        with open(base_params["output_directory"] + full_bilayer_files[0]) as pdb_file:
            content = pdb_file.read()
            assert "ATOM" in content

    def test_tcl_file_contains_triangles(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["create_triangle_tcl_file"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        tcl_file = base_params["output_directory"] + "triangles.tcl"
        with open(tcl_file) as tcl:
            content = tcl.read()
            assert "draw triangle" in content

    def test_grid_points_pdb_contains_x_atoms(self, base_params):
        base_params["surface_equation"] = "z = 0"
        base_params["show_grid_points"] = "TRUE"
        lipidwrapper.run_with_params(base_params)
        grid_file = base_params["output_directory"] + "grid_points.pdb"
        with open(grid_file) as grid:
            content = grid.read()
            assert "ATOM" in content
            assert "X" in content
