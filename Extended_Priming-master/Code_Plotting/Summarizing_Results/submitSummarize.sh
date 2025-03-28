#!/bin/bash
#SBATCH -N 1  
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=sched_mit_arupc_long,sched_mit_arupc,sched_any,sched_mit_hill,newnodes
. /etc/profile.d/modules.sh
module add mit/matlab/2022a
matlab -nodisplay -r "run('"$1"'), exit"
