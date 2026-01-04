## imports

# standard
import os
import random
import string
import gzip
import pickle as pickle


## methods


def save_pickle(item, params: dict, an_id: str = ""):
    """Save an object to a pickle file

    Arguments:
    item -- The object to be saved
    item -- A dictionary, the user-specified comand-line parameters
    an_id -- An optional string, the id of the pickled object

    Returns:
    A string, the id of the current pickle. If the user-specified an_id is '', then a random an_id is generated.

    """

    # make an initial guess at hte pickle name
    filename = params["memory_store_dir"] + an_id + ".pickle"

    # need to make up an id, change filename, if an_id is not specified
    if an_id == "":
        an_id = "".join(
            random.choice(string.ascii_uppercase + string.digits) for x in range(12)
        )
        filename = params["memory_store_dir"] + an_id + ".pickle"

        # keep trying until you come up with a unique an_id. Almost certainly on first try.
        while os.path.exists(filename):
            an_id = "".join(
                random.choice(string.ascii_uppercase + string.digits) for x in range(12)
            )
            filename = params["memory_store_dir"] + an_id + ".pickle"

    while True:  # keep trying to write until you succeed
        try:
            pickle.dump(item, gzip.open(filename, "wb"), protocol=2)
            break
        except:
            pass

    return an_id  # let user know what an_id you settled on


def load_pickle(an_id: str, params: dict):
    """Save an object to a pickle file

    Arguments:
    an_id -- A string, the id of the pickled object to load
    params -- A dictionary, the user-specified comand-line parameters

    Returns:
    A python object, loaded from the pickle file

    """

    while True:  # keep trying to return the pickle until you succeed.
        try:
            return pickle.load(
                gzip.open(params["memory_store_dir"] + an_id + ".pickle", "rb")
            )
        except:
            pass


def openfile(filename: str, mode: str, params: dict):
    if params["compress_output"] == "TRUE":  # open a gzip file
        return gzip.open(filename + ".gz", mode)
    else:  # open a regular file
        return open(filename, mode)
