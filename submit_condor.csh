#! /bin/csh
setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
        echo $i
        source $i
end
source $HOME/.cshrc
source /opt/phenix/core/bin/phenix_setup.csh
set count = 0
foreach file ( /gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/sanghwa/taxi/Run15pp200CAMUP108/18416/data/4*.root) 
    echo $count $file
    @ count = $count + 1
    
end
