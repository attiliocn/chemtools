#!/usr/bin/env python

import os
import argparse
from shutil import copy2, move

parser = argparse.ArgumentParser()
parser.add_argument(
    '--ensemble',
    nargs=1,
    type=str,
    help='Folder containing crest_conformers.output'
)

args = parser.parse_args()

ensembleFolder = args.ensemble[0]
entryName = ensembleFolder.rsplit('/')[0]
ensembleName = entryName+'_crest-conformers.xyz'

MAIN_DIR = os.getcwd()
RESULTS_DIR = os.path.join(MAIN_DIR,ensembleFolder,'macromodel_trial')
SCR_DIR = os.path.join(RESULTS_DIR,'scratch')

for folder in (RESULTS_DIR,SCR_DIR):
    if not os.path.exists(folder):
        os.mkdir(folder)

print('Running '+entryName)
print('Starting clustering')

copy2(os.path.join(ensembleFolder,'crest_conformers.xyz'),RESULTS_DIR)
os.chdir(RESULTS_DIR)
os.system('/usr/bin/obabel -ixyz crest_conformers.xyz -omaegz > crest_conformers.maegz 2> /dev/null')
os.system('/home/attilio/opt/schrodinger2020-3/run conformer_cluster.py crest_conformers.maegz -a all -j cluster -l Average -n 0 -NOJOBID > cluster.output 2> /dev/null')

print('Clustering finished')

maestroClusters = list()
for file in os.listdir():
    if file.endswith('.maegz'):
        maestroClusters.append(file)
maestroClusters.remove('crest_conformers.maegz')

print('Converting Maestro files to xyz')

for file in maestroClusters:
    outputFile = entryName+'_'+file.rsplit('.')[0].rsplit('_',1)[1]+'.xyz'
    os.system('/home/attilio/opt/schrodinger2020-3/utilities/mol2convert -imae '+file+' -omol2 temp 2> /dev/null')
    os.system('/usr/bin/obabel -imol2 temp -oxyz > '+outputFile+' 2> /dev/null')

print('Moving scratch files to scratch')

for file in os.listdir():
    if not file.endswith('.xyz') and not os.path.isdir(file):
        move(file, os.path.join(SCR_DIR,file))

print('Finished '+entryName)
print('')