## imports

# standard
import gzip
import os

# custom

# local
from lipidwrapper import file_io


## classes


class TestSavePickle:
    def test_save_pickle_creates_file(self, temp_directory):
        params = {"memory_store_dir": temp_directory}
        data = {"key": "value", "number": 42}
        pickle_id = file_io.save_pickle(data, params)
        expected_file = os.path.join(temp_directory, f"{pickle_id}.pickle")
        assert os.path.exists(expected_file)

    def test_save_pickle_with_custom_id(self, temp_directory):
        params = {"memory_store_dir": temp_directory}
        data = [1, 2, 3, 4, 5]
        custom_id = "my_custom_id"
        returned_id = file_io.save_pickle(data, params, an_id=custom_id)
        assert returned_id == custom_id
        expected_file = os.path.join(temp_directory, f"{custom_id}.pickle")
        assert os.path.exists(expected_file)

    def test_save_pickle_generates_random_id(self, temp_directory):
        params = {"memory_store_dir": temp_directory}
        data = "test string"
        pickle_id = file_io.save_pickle(data, params)
        assert len(pickle_id) == 12

    def test_pickle_returns_valid_id(self, temp_directory):
        params = {"memory_store_dir": temp_directory}
        data = {"nested": {"data": [1, 2, 3]}}
        pickle_id = file_io.save_pickle(data, params)
        assert isinstance(pickle_id, str)
        assert len(pickle_id) > 0

    def test_multiple_pickles_unique_ids(self, temp_directory):
        params = {"memory_store_dir": temp_directory}
        ids = set()
        for i in range(5):
            pickle_id = file_io.save_pickle(f"data_{i}", params)
            ids.add(pickle_id)
        assert len(ids) == 5


class TestOpenFile:
    def test_openfile_regular_write(self, temp_directory):
        params = {"compress_output": "FALSE"}
        filepath = os.path.join(temp_directory, "test.txt")
        with file_io.openfile(filepath, "w", params) as f:
            f.write("test content")
        assert os.path.exists(filepath)
        with open(filepath, "r") as f:
            content = f.read()
        assert content == "test content"

    def test_openfile_regular_read(self, temp_directory):
        params = {"compress_output": "FALSE"}
        filepath = os.path.join(temp_directory, "test.txt")
        with open(filepath, "w") as f:
            f.write("read this content")
        with file_io.openfile(filepath, "r", params) as f:
            content = f.read()
        assert content == "read this content"

    def test_openfile_compressed_write(self, temp_directory):
        params = {"compress_output": "TRUE"}
        filepath = os.path.join(temp_directory, "test.txt")
        with file_io.openfile(filepath, "wt", params) as f:
            f.write("compressed content")
        assert os.path.exists(filepath + ".gz")

    def test_openfile_compressed_read(self, temp_directory):
        params = {"compress_output": "TRUE"}
        filepath = os.path.join(temp_directory, "test.txt")
        with gzip.open(filepath + ".gz", "wt") as f:
            f.write("compressed read test")
        with file_io.openfile(filepath, "rt", params) as f:
            content = f.read()
        assert content == "compressed read test"
