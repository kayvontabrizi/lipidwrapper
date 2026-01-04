## imports

# standard
import sys
import random
import multiprocessing


## classes


class multi_threading:
    """A class for running calculations on multiple processors"""

    results = []

    def __init__(
        self, inputs, num_processors, task_class, params, progress_bar_prefix=""
    ):
        """Launches a calculation on multiple processors

        Arguments:
        inputs -- A list, containing all the input required for the calculation
        num_processors -- An integer, the requested number of processors to use
        task_class -- An class, the class governing what calculations will be run on a given thread
        progress_bar_prefix -- An optional string, the prefix to append to the progress bar

        Returns:
        Nothing, though the objects self.results list is populated with the calculation results

        """

        # set up the progress bar
        self.results = []
        indices_to_star = []
        if len(inputs) < 50:
            indices_to_star = list(range(len(inputs)))
        else:
            while len(indices_to_star) < 50:
                indx_to_add = random.choice(list(range(len(inputs))))
                if not indx_to_add in indices_to_star:
                    indices_to_star.append(indx_to_add)

        if progress_bar_prefix != "":
            toadd = 78 - len(progress_bar_prefix) - 50
            progress_bar_prefix = progress_bar_prefix + (" " * toadd)
            sys.stdout.write(progress_bar_prefix)

        if num_processors == 1:  # so just running on 1 processor, perhaps under windows
            single_thread = task_class()
            single_thread.total_num_tasks = len(inputs)
            single_thread.indices_to_star = indices_to_star

            single_thread.results = []
            for item in inputs:
                single_thread.value_func(item, None)

            self.results = single_thread.results

        else:  # so it actually is running on multiple processors

            cpu_count = 1
            cpu_count = multiprocessing.cpu_count()

            # first, if num_processors <= 0, determine the number of processors to use programatically
            if num_processors <= 0:
                num_processors = cpu_count

            # reduce the number of processors if too many have been specified
            if len(inputs) < num_processors:
                num_processors = len(inputs)

            if len(inputs) == 0:  # if there are no inputs, there's nothing to do.
                self.results = []
                return

            # now, divide the inputs into the appropriate number of processors
            inputs_divided = {}
            for t in range(num_processors):
                inputs_divided[t] = []

            for t in range(0, len(inputs), num_processors):
                for t2 in range(num_processors):
                    index = t + t2
                    if index < len(inputs):
                        inputs_divided[t2].append(inputs[index])

            # now, run each division on its own processor
            running = multiprocessing.Value("i", num_processors)
            mutex = multiprocessing.Lock()

            arrays = []
            threads = []
            for i in range(num_processors):
                athread = task_class()
                athread.total_num_tasks = len(inputs)
                athread.indices_to_star = indices_to_star

                threads.append(athread)
                arrays.append(multiprocessing.Array("i", [0, 1]))

            results_queue = multiprocessing.Queue()  # to keep track of the results

            processes = []
            for i in range(num_processors):
                p = multiprocessing.Process(
                    target=threads[i].runit,
                    args=(running, mutex, results_queue, inputs_divided[i]),
                )
                p.start()
                processes.append(p)

            while running.value > 0:
                is_running = 0  # wait for everything to finish

            # compile all results into one list
            for thread in threads:
                chunk = results_queue.get()
                self.results.extend(chunk)

        if progress_bar_prefix != "":
            print()  # because the progress bar is now done


class general_task:
    """A parent class of others that governs what calculations are run on each thread"""

    results = []

    def print_star_if_appropriate(self, current_index: int):
        """If appropriate, prints an asterix as part of the progress bar

        Arguments:
        current_index -- An integer, the index of the current calculation

        """

        # if the current index is one of the ones that you should star, write out the asterix
        if current_index in self.indices_to_star:
            sys.stdout.write("*")

    def runit(self, running, mutex, results_queue, items):
        """Launches the calculations on this thread

        Arguments:
        running -- A multiprocessing.Value object
        mutex -- A multiprocessing.Lock object
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        items -- A list, the input data required for the calculation

        """

        for item in items:
            self.value_func(item, results_queue)

        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)
