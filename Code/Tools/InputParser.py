import pyBrellaSampling.Code.Tools.io as io
from pprint import pprint
import os

def ParseInputs(Defaults:dict):
    pprint(Defaults)
    path = os.path.join(Defaults["workdir"],Defaults["configfile"])
    try:
        InputFile = io.jsonRead(path)
    except FileNotFoundError:
        print("WARNING: Input file not found... Using defaults" if Defaults["verbosity"] > 0 else "")
        InputFile = Defaults
        io.jsonDump(InputFile, path)
    for key in Defaults.keys():
        try:
            Defaults[key] = InputFile[key]
        except KeyError:
            pass
    return Defaults

def check_keys(dictionary:dict, keys:list ):
    for key in keys:
        if key not in dictionary:
            dictionary[key] = None
    return dictionary