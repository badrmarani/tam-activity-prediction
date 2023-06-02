#!/bin/bash

NUM_CORES=$(grep -c ^processor /proc/cpuinfo);

ENZ=$1;
SUB=$2;
OUTPUT=$3;
PMP="PMP.pdb";


echo "Enzyme:   $ENZ";
echo "Substrat: $SUB";
echo "Cofactor (PMP): $PMP";


rm -rf ${OUTPUT} 
mkdir ${OUTPUT};

####################################################
# ENZYME-PMP DOCKING
# ENZYME-PMP DOCKING
# ENZYME-PMP DOCKING

# obabel -:"$PMP" -O "${OUTPUT}/PMP.sdf" --gen3D
# obabel "${OUTPUT}/PMP.sdf" -O "${OUTPUT}/PMP.pdb" --gen3D

obabel PMP.pdb -O ${OUTPUT}/PMP.pdbqt --gen3D

ENZYME_NAME=$(echo $ENZ | sed "s/\.[^.]*$//");
grep ATOM $ENZ > ${OUTPUT}/${ENZYME_NAME}.pdb

obabel ${OUTPUT}/${ENZYME_NAME}.pdb -O ${OUTPUT}/${ENZYME_NAME}.pdb   -gen3d
obabel ${OUTPUT}/${ENZYME_NAME}.pdb -O ${OUTPUT}/${ENZYME_NAME}.pdbqt -xr

gnina \
    --cnn_scoring rescore \
    --cnn crossdock_default2018_4 \
    --exhaustiveness $NUM_CORES \
    --scoring vinardo \
    -r ${OUTPUT}/${ENZYME_NAME}.pdbqt \
    -l ${OUTPUT}/PMP.pdbqt \
    --autobox_ligand ${OUTPUT}/PMP.pdbqt \
    --autobox_add 20 \
    --seed 0 \
    -o ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMPs.sdf

# obabel ${OUTPUT}/PMP.pdbqt -O ${OUTPUT}/PMP.pdb

mkdir ${OUTPUT}/temp/
obabel ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMPs.sdf -O ${OUTPUT}/temp/${ENZYME_NAME}_DOCKED_PMP.pdb -m
cp ${OUTPUT}/temp/${ENZYME_NAME}_DOCKED_PMP1.pdb ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMP.pdb

python3 combine_pdbs.py \
    --pdb1 ${OUTPUT}/${ENZYME_NAME}.pdb \
    --pdb2 ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMP.pdb \
    --outp ${OUTPUT}/${ENZYME_NAME}_PMP.pdb \


####################################################
# ENZYME-PMP-SUBSTRAT DOCKING
# ENZYME-PMP-SUBSTRAT DOCKING
# ENZYME-PMP-SUBSTRAT DOCKING

SUBSTRAT_NAME=$(echo $SUB | sed "s/\.[^.]*$//");
obabel $SUB -O ${OUTPUT}/${SUBSTRAT_NAME}.pdbqt -gen3d

obabel ${OUTPUT}/${ENZYME_NAME}_PMP.pdb -O ${OUTPUT}/${ENZYME_NAME}_PMP.pdb   -gen3d
obabel ${OUTPUT}/${ENZYME_NAME}_PMP.pdb -O ${OUTPUT}/${ENZYME_NAME}_PMP.pdbqt -gen3d -xr

obabel ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMP.pdb -O ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMP.pdbqt -gen3d -xr

# --autobox_ligand ${OUTPUT}/${SUBSTRAT_NAME}.pdbqt \

# --autobox_ligand ${OUTPUT}/${ENZYME_NAME}_PMP.pdbqt \
gnina \
    --cnn_scoring rescore \
    --cnn crossdock_default2018_4 \
    --exhaustiveness $NUM_CORES \
    --scoring vinardo \
    -r ${OUTPUT}/${ENZYME_NAME}_PMP.pdbqt \
    -l ${OUTPUT}/${SUBSTRAT_NAME}.pdbqt \
    --autobox_ligand ${OUTPUT}/${SUBSTRAT_NAME}.pdbqt \
    --autobox_add 6 \
    --seed 0 \
    -o ${OUTPUT}/${ENZYME_NAME}_PMP_DOCKED_SUBs.sdf

obabel ${OUTPUT}/${ENZYME_NAME}_PMP_DOCKED_SUBs.sdf -O ${OUTPUT}/temp/${ENZYME_NAME}_PMP_DOCKED_SUB.pdb -m
cp ${OUTPUT}/temp/${ENZYME_NAME}_PMP_DOCKED_SUB1.pdb ${OUTPUT}/${ENZYME_NAME}_PMP_DOCKED_SUB.pdb

python3 combine_pdbs.py \
    --pdb1 ${OUTPUT}/${ENZYME_NAME}_PMP.pdb \
    --pdb2 ${OUTPUT}/${ENZYME_NAME}_PMP_DOCKED_SUB.pdb \
    --outp ${OUTPUT}/${ENZYME_NAME}_PMP_${SUBSTRAT_NAME}.pdb \



# CLEAN CLEAN CLEAN
rm -rf ${OUTPUT}/temp/
rm ${OUTPUT}/${ENZYME_NAME}_DOCKED_PMP.pdb
rm ${OUTPUT}/${ENZYME_NAME}_PMP_DOCKED_SUB.pdb
rm ${OUTPUT}/*.pdbqt
