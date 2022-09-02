#!/bin/csh
source /opt/phenix/core/bin/phenix_setup.csh
root -b -q -l  '/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/HF_Analysis/start_analysis/TSSA_OHF/SingleRun_Variable_histograms.C("'${1}'" , "'output'")'
