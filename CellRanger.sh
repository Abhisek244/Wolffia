#!/bin/bash

#SBATCH -J CR_Wmic            # Job name
#SBATCH -o CR_Wmic.%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH -c 24                    # Total number of cores for the single task
#SBATCH --mem 220G                # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 72:00:00              # Run time (hh:mm:ss)

#SBATCH -A m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END
#SBATCH --mail-user=achakrab@uni-mainz.de

## Load all necessary modules if needed
module load bio/CellRanger/7.1.0
cellranger mkref --genome=Haplotype1_Wmic --fasta=Wmic_reference_genome.fasta --genes=Wmic_reference_genome.gtf --memgb=220 --nthreads=24
cellranger count --id=Wmic_count --transcriptome=Haplotype1_Wmic --fastqs=/lustre/project/m2_jgu-evoltroph/achakrab/microscopica/snRNA_raw --sample=jgu_xu_2023_01_01_WH_230328A,jgu_xu_2023_01_02_WH_230328B --include-introns=true --chemistry=ARC-v1 --localcores=24 --localmem=220
