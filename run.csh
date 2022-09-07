#!/bin/bash
source /opt/phenix/core/bin/phenix_setup.csh
path=/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/HF_Analysis/start_analysis/data/
echo $path 
file="new_run_list.txt" #the file where you keep your string name

while IFS= read -r line
do
    file_id=$line
    echo $line
    root -b -q -l 'TSSA_OHF/SingleRun_Variable_histograms.C("'$path$line'","'$line'")'
done < "$file"
