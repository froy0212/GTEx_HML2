#!/bin/bash
#SBATCH --job-name=GTEX_Pipeline		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=batch
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=1		# Number of CPU cores per task
#SBATCH --mem-per-cpu=5000
#SBATCH --time=1-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source /cluster/tufts/micr_coffin01/aburn01/Tools/Start_Variables.txt
echo ${File}

echo "Running Download of ${File}"
sh over_anvil_download.sh 

echo "Running bamtofastq of ${File}"
sbatch -W over_wrapper_bamtofastq.sh

