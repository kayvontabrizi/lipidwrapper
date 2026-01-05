## imports

# standard
import gzip
import os
import pickle
import random
import string
import typing


## methods


def save_pickle(item: typing.Any, params: dict, an_id: str = "") -> str:
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


def load_pickle(an_id: str, params: dict) -> typing.Any:
    while True:
        try:
            return pickle.load(
                gzip.open(params["memory_store_dir"] + an_id + ".pickle", "rb")
            )
        except:
            pass


def openfile(
    filename: str, mode: str, params: dict
) -> typing.Union[gzip.GzipFile, typing.IO]:
    if params["compress_output"] == "TRUE":
        return gzip.open(filename + ".gz", mode)
    else:
        return open(filename, mode)
