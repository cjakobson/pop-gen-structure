#!/bin/bash

for f in /Users/cjakobson/Dropbox/JaroszLab/pop-gen-structure-dependencies/alphafold-predictions/*.pdb; do
echo $f
filearray=(${f//\// })
echo "/Users/cjakobson/Dropbox/JaroszLab/pop-gen-structure-dependencies/dssp-output/${filearray[6]}dssp.txt"
mkdssp -i $f -o "/Users/cjakobson/Dropbox/JaroszLab/pop-gen-structure-dependencies/dssp-output/${filearray[6]}dssp.txt" -v
done