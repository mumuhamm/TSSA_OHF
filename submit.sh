#! /bin/sh
for file in /gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/sanghwa/taxi/Run15pp200CAMUP108/18416/data/4*.root ; do
    echo $file
    purefile='run-XV $file'
    stem='run-XV $purefile.root'
    condor_submit \
    -a "output = $stem.out" \
    -a "error = $stem.error" \
    -a "Log = $stem.log" \
    -a "Arguments = $file" \
    condor_v1
done

