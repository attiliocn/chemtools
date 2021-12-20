#!/usr/bin/env python

import argparse
import re

def splitConformerEnsemble(conformerEnsemble):
    compoundIdentifier = conformerEnsemble.rsplit('_',1)[0]
    conformersFilenames = list()
    with open(conformerEnsemble) as conformers:
        data = conformers.read()
        nAtoms = int(data.split()[0])
        
        regexAtomNumber = '^\s*'+ str(nAtoms)
        dataSplitted = re.split(regexAtomNumber,data, flags=re.MULTILINE)
        
        for conformerNumber in range(1,len(dataSplitted)):
            conformerFilename = compoundIdentifier+'_conf-'+str(conformerNumber)+'.xyz'
            conformersFilenames.append(conformerFilename)
            with open(conformerFilename, "w") as f:
                f.write(str(nAtoms))
                f.write(dataSplitted[conformerNumber])
    return conformersFilenames

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ensemble',
        nargs=1,
        type=str,
        help='xyz file containing all conformers'
    )

    args = parser.parse_args()
    ensembleFile = args.ensemble[0]
    splitConformerEnsemble(ensembleFile)
