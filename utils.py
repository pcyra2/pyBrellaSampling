import numpy
import logging as log
import argparse as ap
import json


def file_read(path):
    with open(path, 'r') as file:
        data = file.readlines()
    return data

def file_write(path, lines):
    with open(path, 'w') as f:
        for i in lines:
            print(i,file=f)
def data_2d(path):
    data0 = file_read(path)
    data = []
    for i in data0:
        if "#" not in i:
           data.append(i)
        if " 4000 " in i: ### For cutting down to 2 ps of data
            break
    col1 = numpy.zeros(len(data))
    col2 = numpy.zeros(len(data))
    for i in range(len(data)):
        words = data[i].split()
        col1[i] = words[0]
        col2[i] = words[1]
    return col1, col2

def dict_read(path):
    with open(path, "r") as f:
        data = json.load(f)
    return data

def dict_write(path, dict):
    with open(path,"w") as f:
        json.dump(dict, f)
def input_parser(args):
    if args.verbose == 0:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.CRITICAL)
    elif args.verbose == 1:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
        log.info("Verbose output.")
    try:
        data = dict_read(f"{args.WorkDir}{args.input}")
        lst = vars(args)
        newlst = lst
        for i in data:
            for j in lst:
                if j == i:
                    newlst[j] = data[i]
        args = ap.Namespace(**newlst)
        dict_write(f"{args.WorkDir}{args.input}", newlst)
    except FileNotFoundError:
        lst = vars(args)
        dict_write(f"{args.WorkDir}{args.input}", lst)
    return args

