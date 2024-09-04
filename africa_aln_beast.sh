#!/bin/bash

#SBATCH --job-name Africa_RABV
#SBATCH -p gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH -o beast_Africa_RABV_aln_1.out -e beast_Africa_RABV_aln_1.err
#SBATCH --nodes=1
#SBATCH --mail-user=andrew.holtz@pasteur.fr --mail-type=BEGIN,END,FAIL


module purge
module load graalvm/ce-java19-22.3.1 beagle-lib/4.0.1
module load beast-mcmc/10.5.0-beta4


cd /pasteur/zeus/projets/p02/NGS_Viroscreen/aholtz/Africa/BEAST/Beast_2
xml_file=$1
beast -beagle_cuda -beagle_GPU -beagle_double -save_every 10000 -save_state beast_Africa_RABV_1_checkpoint.state Africa_RABV_sub_aln.xml


exit 0
