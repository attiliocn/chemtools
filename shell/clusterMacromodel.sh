#!/bin/bash

conformationalSearchFolders=$@
MAIN_FOLDER=$PWD

for folder in $conformationalSearchFolders 
do
    entry=${folder%/*}
    conformerEnsemble=""$entry"_crest-conformers.xyz"
    resultsFolder="macromodel-CL_trial"

    cd $folder
    mkdir $resultsFolder
    cp crest_conformers.xyz $resultsFolder/$conformerEnsemble
    cd $resultsFolder

    #for i in $(seq $(head -n 1 $conformerEnsemble)); do echo -n $i, >> atomlist.csv; done
    #sed -i 's/,$/\n/' atomlist.csv

    obabel -ixyz $conformerEnsemble -omaegz > ${conformerEnsemble%.*}.maegz
    /home/attilio/opt/schrodinger2020-3/run conformer_cluster.py ${conformerEnsemble%.*}.maegz -a all -j cluster -l Average -n 0
    cd ..
    cd $resultsFolder
    ls
    ls
    dir
done









exit

for i in $(seq $(head -n 1 crest_conformers.xyz)); do echo -n $i, >> atomlist.csv; done
sed -i 's/,$/\n/' atomlist.csv

obabel -ixyz crest_conformers.xyz -omaegz > inputmacromodel.maegz
run conformer_cluster.py inputmacromodel.maegz -n 0 -a_rms atomlist.csv

/home/attilio/opt/schrodinger2020-3/utilities/mol2convert -imae conf_cluster_ligand1_representatives.maegz -omol2 temp
obabel -imol2 temp -oxyz