#! /bin/sh
count=0
for file in /gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/HF_Analysis/muon_output/*.root ; do
    echo $count $file
    count=$((count+1))
    purefile='basename $file'
    stem='analysis_job'
    condor_submit \
    -a "output = $stem.out" \
    -a "error = $stem.error" \
    -a "Log = $stem.log" \
    -a "Arguments = $file $count" \
    condor_v2_analysis
done

