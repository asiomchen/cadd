#!/bin/bash
exhaustiveness=160
PROTEIN_PATH=cox.pdbqt
LIGANDS_PATH=./pocket_set/*
OUTPUT_PATH=./docked_pocket_160
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
  ./qvina-w --receptor $PROTEIN_PATH --ligand $ligand --num_modes $num_modes --exhaustiveness $exhaustiveness --seed 42 --out $output --center_x 30.43100070953369 --center_y 28.46949863433838 --center_z 21.786001205444336 --size_x 89.71600151062012 --size_y 105.20899772644043 --size_z 101.20000076293945
  end_time=$(date +%s)
  echo "Time elapsed: $((end_time-start_time)) seconds"
  num_ligands=$((num_ligands-1))
  echo "Estimated time remaining: $((num_ligands*(end_time-start_time)/60)) minutes"
done
