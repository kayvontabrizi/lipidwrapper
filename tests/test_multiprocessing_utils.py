## imports

# local
from lipidwrapper import multiprocessing_utils


## classes


class TestGeneralTask:
    def test_init(self):
        task = multiprocessing_utils.general_task()
        assert hasattr(task, "results")

    def test_results_list(self):
        task = multiprocessing_utils.general_task()
        task.results = []
        task.results.append("test")
        assert len(task.results) == 1


class TestMultiThreading:
    class SimpleTask(multiprocessing_utils.general_task):
        def value_func(self, item, results_queue):
            self.results.append(item * 2)

    def test_single_processor(self):
        inputs = [1, 2, 3, 4, 5]
        params = {}
        mt = multiprocessing_utils.multi_threading(
            inputs, 1, self.SimpleTask, params, ""
        )
        expected = [2, 4, 6, 8, 10]
        assert sorted(mt.results) == expected

    def test_empty_inputs(self):
        inputs = []
        params = {}
        mt = multiprocessing_utils.multi_threading(
            inputs, 1, self.SimpleTask, params, ""
        )
        assert mt.results == []

    def test_single_input(self):
        inputs = [5]
        params = {}
        mt = multiprocessing_utils.multi_threading(
            inputs, 1, self.SimpleTask, params, ""
        )
        assert mt.results == [10]


class TestMultiThreadingWithStrings:
    class StringTask(multiprocessing_utils.general_task):
        def value_func(self, item, results_queue):
            self.results.append(item.upper())

    def test_string_processing(self):
        inputs = ["hello", "world"]
        params = {}
        mt = multiprocessing_utils.multi_threading(
            inputs, 1, self.StringTask, params, ""
        )
        assert sorted(mt.results) == ["HELLO", "WORLD"]


class TestMultiThreadingWithTuples:
    class TupleTask(multiprocessing_utils.general_task):
        def value_func(self, item, results_queue):
            a, b = item
            self.results.append(a + b)

    def test_tuple_processing(self):
        inputs = [(1, 2), (3, 4), (5, 6)]
        params = {}
        mt = multiprocessing_utils.multi_threading(
            inputs, 1, self.TupleTask, params, ""
        )
        expected = [3, 7, 11]
        assert sorted(mt.results) == expected
