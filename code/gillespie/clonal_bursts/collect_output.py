# This is a script to collect the output after many bootstraps.
# Basically, if you have a directory with a bunch of files called
# outputX[.csv], you can use this to compile them all into a single
# file called collected_output.csv.  This is useful to reduce the
# number of files hanging around.

# Note the nrows argument; if more than 1000 you may run into issues
# because it will delete the old files.

import glob
import numpy as np
import pandas as pd
import os

def main():
    output = pd.DataFrame()
    for file in glob.glob("./output*"):
        x = pd.read_csv(file,delimiter=' ',index_col=None, header=None, nrows=1000)
        output = pd.concat((output,x))
    return output

if __name__ == "__main__":
    temp = main()
    temp.to_csv('collected_output.csv',header=False,index=False)
    for file in glob.glob("./output*"):
        os.remove(file)
