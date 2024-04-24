#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os, re
import argparse
from morfeus import read_xyz


# args parsing #
################
argParser = argparse.ArgumentParser()
argParser.add_argument('-n','--mcb', type=int, nargs='+', required=True, help='list of the mcb starting from the Mo followed by the carbene number')
args = argParser.parse_args()
##############

print(args.mcb)

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key=alphanum_key) 

xyzFiles = []
for file in os.listdir(os.getcwd()):
    if file.endswith('.xyz'):
        xyzFiles.append(file)
xyzFiles = natural_sort(xyzFiles)

def getBondLenght(atom1, atom2):
    bondLenght = np.linalg.norm(coordinates[atom1]-coordinates[atom2])
    return bondLenght

geometricData = {}
for xyzFile in xyzFiles:
    if os.stat(xyzFile).st_size==0:
        continue
    else:
        elements,coordinates = read_xyz(xyzFile)
        mcbAtomNumbers = args.mcb
        mcbAtomNumbers = [x-1 for x in mcbAtomNumbers]
        atomNames = ['molybdenum', 'methylidene', 'carbon1', 'carbon2']

    geometricData[xyzFile] = {}

    for i in range(len(mcbAtomNumbers)):
        if i == len(mcbAtomNumbers)-1:
            param = atomNames[i]+'-'+atomNames[0]
            geometricData[xyzFile][param] = getBondLenght(mcbAtomNumbers[i],mcbAtomNumbers[0])
            #print(atomNames[i]+'-'+atomNames[0],'{:.3f}'.format(getBondLenght(mcbAtomNumbers[i],mcbAtomNumbers[0])))
        else:
            param = atomNames[i]+'-'+atomNames[i+1]
            geometricData[xyzFile][param] = getBondLenght(mcbAtomNumbers[i],mcbAtomNumbers[i+1])
            #print(atomNames[i]+'-'+atomNames[i+1],'{:.3f}'.format(getBondLenght(mcbAtomNumbers[i],mcbAtomNumbers[i+1])))

df = pd.DataFrame(geometricData)
df.T.to_csv('geometricData.csv')
print("Done")