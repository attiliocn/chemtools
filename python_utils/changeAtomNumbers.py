#!/usr/bin/env python3
import morfeus
import argparse
import numpy as np
from numpy.core.numeric import tensordot
import itertools

############################################################
# Parsing args                                             #
############################################################
parser = argparse.ArgumentParser()

parser.add_argument('files', nargs='+')
parser.add_argument('-n', type=int, nargs='+',help='Provide the original number followed by the desired number')
parser.add_argument('--auto-mcb', action='store_true')
args = parser.parse_args()

if not args.auto_mcb:
    atomsPairs = list()
    for atomPair in zip(args.n[::2], args.n[1::2]):
        atomsPairs.append(atomPair)

############################################################
# Functions                                                #
############################################################    
def sort_pairs(listOfPairs):
    pass

def print_header(title=None):
    print(15*'=',title,15*'=')

def getMetalCenter(listOfElements):
    metals = ['Mo','Re','W']
    for metal in metals:
        if metal in listOfElements:
            return metal

############################################################
# MAIN                                                     #
############################################################ 

for file in args.files:
    print('='*80)
    print(f"FILENAME:\t{file}")
    elements, coords = morfeus.read_xyz(file)
    atomQuantity = len(coords)
    
    if args.auto_mcb:

        #determine which atoms are part of the MCB
        conmat = morfeus.utils.get_connectivity_matrix(coords,elements, scale_factor=1.1)
        metalCenter = getMetalCenter(elements)

        atomsAndNeighours = dict()

        metalAtomNumber = int(np.where(elements == metalCenter)[0])
        metalCarbonNeighbours = [i for i in range(len(elements)) if (conmat[i][metalAtomNumber] and elements[i] == 'C')]
        atomsAndNeighours[metalAtomNumber] = {'element':elements[metalAtomNumber],'neighbours':metalCarbonNeighbours}

        for i in metalCarbonNeighbours:
            neighboursOfNeighbours = np.where(conmat[i] == 1)[0]
            atomsAndNeighours[i] = {'element':elements[i],'neighbours':neighboursOfNeighbours.tolist()}

        atomsOfInterest = [atomsAndNeighours[i]['neighbours'] for i in atomsAndNeighours.keys()] #returns a list of lists
        atomsOfInterest = set((itertools.chain.from_iterable(atomsOfInterest))) #flatten the list of lists than convert the
                                                                                #resulting list into a set to remove duplicates
        atomsOfInterest = list(atomsOfInterest) #converts back to a list

        neighboursConMat = [conmat[i,atomsOfInterest] for i in atomsAndNeighours[metalAtomNumber]['neighbours']] 
        neighboursConMat = np.array(neighboursConMat)    

        for atomIndex,atomNeighbours in enumerate(neighboursConMat.T):
            sharedNeighboursNum = np.count_nonzero(atomNeighbours == 1)
            if sharedNeighboursNum == 2:
                cycleThridMember = atomsOfInterest[atomIndex]
        
        thridMemberCarbonNeighbours = [i for i in range(len(elements)) if (conmat[i][cycleThridMember] and elements[i] == 'C')]
        atomsAndNeighours[cycleThridMember] = {'element':elements[cycleThridMember],'neighbours':thridMemberCarbonNeighbours}

        for i in atomsAndNeighours[cycleThridMember]['neighbours']:
            closeNeighbours = elements[atomsAndNeighours[i]['neighbours']]
            numOfHydrogens = np.count_nonzero(closeNeighbours == 'H')
            if numOfHydrogens == 2:
                cycleMethylidene = i
        
        cycleLastMember = [x for x in atomsAndNeighours[cycleThridMember]['neighbours'] if x != cycleMethylidene][0]
        
        mcbNumbers = [x+1 for x in [metalAtomNumber,cycleMethylidene,cycleThridMember,cycleLastMember]]
        ############################################################

        print(f"MODE:\t\tAUTOMCB")
        print(f'DETECTED MCB:\t{mcbNumbers}')

        #check if the reaction is an ethene methatesis
        methylideneNeighbours = set(atomsAndNeighours[cycleMethylidene]['neighbours'])
        lastMemberNeighbours = set(atomsAndNeighours[cycleLastMember]['neighbours'])
        unique1 = methylideneNeighbours.difference(lastMemberNeighbours) # substituents on methylidene
        unique2 = lastMemberNeighbours.difference(methylideneNeighbours) # substituents on last member
        unique = unique1.union(unique2)
        substituentsElements = elements[list(unique)]
        if np.all(substituentsElements == "H"):
            print(f'ALERT: {file} IS AN ETHENE METATHESIS!!! ')
            print('**')

        atomsPairs = [(l,k+1) for k,l in enumerate(mcbNumbers)]
        print(f'Provided exchange pairs: {atomsPairs}')
        ############################################################
    else:
        print(f"MODE:\t\tMANUAL EXCHANGE")
        print(f'Provided exchange pairs: {atomsPairs}')
        print('*'*80)

# pairwise corretion for internal exchanges
# triggers only when DEST from one pair matches SOURCE on another pair
# more testing is needed...
    for firstPairIndex in range(len(atomsPairs)):
        for secondPairIndex in range(firstPairIndex+1, len(atomsPairs)):
            intersection = set(atomsPairs[firstPairIndex]).intersection(set(atomsPairs[secondPairIndex])) 
            if intersection:
                intersection = [i for i in intersection][0]
                indexA = atomsPairs[firstPairIndex].index(intersection)
                indexB = atomsPairs[secondPairIndex].index(intersection)
                fixIndex = 0 if indexA else 1

                if indexA:
                    print("WARNING: Internal exchange detected")
                    print(f"\tFirst exchange pair: {atomsPairs[firstPairIndex]}, Second exchange pair: {atomsPairs[secondPairIndex]}")
                    print(f"\tReplacing {atomsPairs[secondPairIndex][indexB]} for {atomsPairs[firstPairIndex][fixIndex]} in the second pair\n")

                    temp = list(atomsPairs[secondPairIndex])
                    temp[indexB] = atomsPairs[firstPairIndex][fixIndex]
                    atomsPairs[secondPairIndex] = tuple(temp)
                    print(f"Corrected exchange pairs: {atomsPairs}")
    
    print('*'*80)
    print("Summary:\n")
    for atomPair in atomsPairs:
        atom0 = atomPair[0]-1
        atom1 = atomPair[1]-1

        element0 = np.copy(elements[atom0])
        coords0 = np.copy(coords[atom0])
        element1 = np.copy(elements[atom1])
        coords1 = np.copy(coords[atom1])

        print(f"\tSwapping {element0}{atom0+1} ==> {element1}{atom1+1}")

        elements[atom0], elements[atom1] = element1, element0
        coords[atom0], coords[atom1] = coords1, coords0

    split = file.rsplit('_',1)
    newFilename = split[0]+'_numbered_'+split[1]
    with open(newFilename, 'w') as newFile:
        newFile.write(str(atomQuantity)+'\n')
        newFile.write(file+'\n')

        for i in range(atomQuantity):
            line = np.append(np.array(elements[i]),coords[i].astype(str))
            newFile.write(' '.join(line)+'\n')
    print("\nAll done!")
    print('='*80)

