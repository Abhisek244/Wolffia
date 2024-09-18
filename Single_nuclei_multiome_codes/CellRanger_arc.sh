#!/bin/bash

#SBATCH -J Wmic_ARC            # Job name
#SBATCH -o Wmic_ARC.%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH -c 24                    # Total number of cores for the single task
#SBATCH --mem 220G                # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 72:00:00              # Run time (hh:mm:ss)

#SBATCH -A m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END
#SBATCH --mail-user=achakrab@uni-mainz.de

## Load all necessary modules if needed
cellranger-arc mkref --config=Wmic_Haplotype1.config --memgb=220 --nthreads=24
cellranger-arc count --id=Wmic_multiome --reference=Haplotype1_Wmic --libraries=libraries.csv --localcores=24 --localmem=220
