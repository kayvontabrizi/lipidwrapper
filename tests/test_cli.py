## imports

# standard
import tempfile
import os

# custom
import pytest

# local
from lipidwrapper import cli


class TestGetCommandlineParameters:
    @pytest.fixture
    def temp_lipid_pdb(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(
                "ATOM      1  P   POPC    1       0.000   0.000   0.000  1.00  0.00           P\n"
            )
            f.write(
                "ATOM      2  C1  POPC    1       1.000   0.000   0.000  1.00  0.00           C\n"
            )
            temp_path = f.name
        yield temp_path
        os.unlink(temp_path)

    def test_default_parameters(self, temp_lipid_pdb):
        argv = ["script_name", f"--lipid_pdb_filename={temp_lipid_pdb}"]
        params = cli.get_commandline_parameters(argv)
        assert params["delete_clashing_lipids"] == "FALSE"
        assert params["fill_holes"] == "TRUE"
        assert params["clash_cutoff"] == 2.0
        assert params["number_of_processors"] == 1

    def test_surface_equation_default(self, temp_lipid_pdb):
        argv = ["script_name", f"--lipid_pdb_filename={temp_lipid_pdb}"]
        params = cli.get_commandline_parameters(argv)
        assert "numpy.sin" in params["surface_equation"]

    def test_custom_surface_equation(self, temp_lipid_pdb):
        custom_eq = "z = x + y"
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            f"--surface_equation={custom_eq}",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["surface_equation"] == custom_eq

    def test_min_max_parameters(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--min_x=100",
            "--max_x=200",
            "--min_y=150",
            "--max_y=250",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["min_x"] == 100.0
        assert params["max_x"] == 200.0
        assert params["min_y"] == 150.0
        assert params["max_y"] == 250.0

    def test_step_parameters(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--step_x=10",
            "--step_y=15",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["step_x"] == 10.0
        assert params["step_y"] == 15.0

    def test_delete_clashing_lipids_true(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--delete_clashing_lipids=TRUE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["delete_clashing_lipids"] == "TRUE"

    def test_delete_clashing_lipids_lowercase(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--delete_clashing_lipids=true",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["delete_clashing_lipids"] == "TRUE"

    def test_clash_cutoff(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--clash_cutoff=1.5",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["clash_cutoff"] == 1.5

    def test_fill_holes_false(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--fill_holes=FALSE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["fill_holes"] == "FALSE"

    def test_number_of_processors(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--number_of_processors=4",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["number_of_processors"] == 4

    def test_fill_hole_exhaustiveness(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--fill_hole_exhaustiveness=20",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["fill_hole_exhaustiveness"] == 20

    def test_output_directory(self, temp_lipid_pdb, temp_directory):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            f"--output_directory={temp_directory}",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["output_directory"].endswith("/")

    def test_output_directory_adds_slash(self, temp_lipid_pdb):
        with tempfile.TemporaryDirectory() as tmpdir:
            argv = [
                "script_name",
                f"--lipid_pdb_filename={temp_lipid_pdb}",
                f"--output_directory={tmpdir}",
            ]
            params = cli.get_commandline_parameters(argv)
            assert params["output_directory"].endswith("/")

    def test_lipid_headgroup_marker_default(self, temp_lipid_pdb):
        argv = ["script_name", f"--lipid_pdb_filename={temp_lipid_pdb}"]
        params = cli.get_commandline_parameters(argv)
        assert len(params["lipid_headgroup_marker"]) == 2
        assert params["lipid_headgroup_marker"][0] == (None, "", None, "P")
        assert params["lipid_headgroup_marker"][1] == (None, "CHL1", None, "O3")

    def test_lipid_headgroup_marker_custom(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--lipid_headgroup_marker=_N,DPPC_P",
        ]
        params = cli.get_commandline_parameters(argv)
        assert len(params["lipid_headgroup_marker"]) == 2
        assert params["lipid_headgroup_marker"][0] == (None, "", None, "N")
        assert params["lipid_headgroup_marker"][1] == (None, "DPPC", None, "P")

    def test_compress_output(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--compress_output=TRUE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["compress_output"] == "TRUE"

    def test_use_disk_instead_of_memory(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--use_disk_instead_of_memory=TRUE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["use_disk_instead_of_memory"] == "TRUE"

    def test_memory_optimization_factor(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--memory_optimization_factor=2",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["memory_optimization_factor"] == 2

    def test_very_distant_lipids_cutoff(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--very_distant_lipids_cutoff=100.0",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["very_distant_lipids_cutoff"] == 100.0

    def test_triangle_center_proximity_cutoff_distance(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--triangle_center_proximity_cutoff_distance=75.0",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["triangle_center_proximity_cutoff_distance"] == 75.0

    def test_clashing_potential_margin(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--clashing_potential_margin=30.0",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["clashing_potential_margin"] == 30.0

    def test_show_grid_points(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--show_grid_points=TRUE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["show_grid_points"] == "TRUE"

    def test_create_triangle_tcl_file(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--create_triangle_tcl_file=TRUE",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["create_triangle_tcl_file"] == "TRUE"

    def test_max_height(self, temp_lipid_pdb):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            "--max_height=50.0",
        ]
        params = cli.get_commandline_parameters(argv)
        assert params["max_height"] == 50.0

    def test_missing_lipid_file_exits(self):
        argv = ["script_name", "--lipid_pdb_filename=nonexistent_file.pdb"]
        with pytest.raises(SystemExit):
            cli.get_commandline_parameters(argv)

    def test_memory_store_dir_set(self, temp_lipid_pdb, temp_directory):
        argv = [
            "script_name",
            f"--lipid_pdb_filename={temp_lipid_pdb}",
            f"--output_directory={temp_directory}",
        ]
        params = cli.get_commandline_parameters(argv)
        assert "memory_store_dir" in params
        assert "store_in_memory.tmp" in params["memory_store_dir"]

    def test_help_exits(self):
        argv = ["script_name", "--help"]
        with pytest.raises(SystemExit):
            cli.get_commandline_parameters(argv)
