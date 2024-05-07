import sys
import time


import pyBrellaSampling.InputParser as input
import pyBrellaSampling.Umbrella as UmbrellaRun
import pyBrellaSampling.Standalone as Standalone
import pyBrellaSampling.Standalone_New as Standalone1


def main():
    starttime = time.time()
    ### Parse user inputs.
    args = input.VariableParser(sys.argv[1:], JT="Umbrella")
    if args.JobType.casefold() == "umbrella":   #JobType is case insensitive
        UmbrellaRun.main(args)
    if args.JobType.casefold() == "mm" or args.JobType.casefold() == "qmmm":
        Standalone.main(args)
    if args.JobType.casefold() == "standalone":
        Standalone1.main(args)
    endtime = time.time()
    print(f"Total time is {endtime - starttime}")

#
