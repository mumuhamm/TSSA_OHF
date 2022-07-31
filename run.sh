#!/bin/bash

input="input_files.txt"
while IFS= read -r line
do
   root -l 'Variable_histograms.C("$line","output")'
done < "$input"



