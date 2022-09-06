#! /bin/sh
count=0
for file in /gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/sanghwa/taxi/Run15pp200CAMUP108/18416/data/4*.root ; do
    echo $count $file
    run=${file%.*}
    run=${run##*/}
    echo $run
    echo $file
    count=$((count+1))
    tag=$count'_'$run
    purefile='basename $file'
    stem='total_job'
    condor_submit \
    -a "output = $stem.out" \
    -a "error = $stem.error" \
    -a "Log = $stem.log" \
    -a "Arguments = $file $tag" \
    condor_v2
done

