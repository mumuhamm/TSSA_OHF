#!/bin/csh
echo ${_CONDOR_SCRATCH_DIR}
echo $0
setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
	echo $i   
	source $i
end
source $HOME/.cshrc
source /opt/phenix/core/bin/phenix_setup.csh
root.exe -b -q -l  '/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/alibordi/HF_Analysis/start_analysis/TSSA_OHF/analyzer.C("'${1}'" , "'${2}'")'
