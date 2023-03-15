#!/bin/bash
exhaustiveness=$1
n_cpus=$2
PROTEIN_PATH=./cox2.pdbqt
LIGANDS_PATH=./meeko_pocket_set/*
OUTPUT_PATH=./docked_pocket_"$exhaustiveness"
num_modes=10
echo 'Protein path: ' $PROTEIN_PATH
echo 'Output path: ' $OUTPUT_PATH
echo "Exhaustiveness: $exhaustiveness"
# create output directory if it doesn't exist
if [ ! -d $OUTPUT_PATH ]; then
    mkdir $OUTPUT_PATH
fi
# count number of ligands
num_ligands=$(ls -1 $LIGANDS_PATH | wc -l)
# Run docking
for ligand in $LIGANDS_PATH; do
  start_time=$(date +%s)
  filename=$(basename $ligand)
  output=$OUTPUT_PATH/$filename
  echo "Docking $(basename $ligand) to $output"
  echo "./qvina-w --receptor $PROTEIN_PATH --ligand $ligand --num_modes $num_modes --exhaustiveness $exhaustiveness --cpu $n_cpus --seed 42 --out $output --center_x 42.8405 --center_y 31.0155 --center_z 32.3135 --size_x 34.751 --size_y 43.7889 --size_z 34.821"
  ./qvina-w --receptor $PROTEIN_PATH --ligand $ligand --num_modes $num_modes --exhaustiveness $exhaustiveness --cpu $n_cpus  --seed 42 --out $output --center_x 42.8405 --center_y 31.0155 --center_z 32.3135 --size_x 34.751 --size_y 43.7889 --size_z 34.821
  end_time=$(date +%s)
  echo "Time elapsed: $((end_time-start_time)) seconds"
  num_ligands=$((num_ligands-1))
  echo "Estimated time remaining: $((num_ligands*(end_time-start_time)/60)) minutes"
done
