#!/bin/bash
#SBATCH -D /home/moritz/repos/thorsellia/
#SBATCH -J selblast
#SBATCH -o /home/moritz/repos/thorsellia/selblast.err
#SBATCH -e /home/moritz/repos/thorsellia/selblast.out
#SBATCH -A b2013086
#SBATCH -t 4-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

module load raxml
python thorsellia/__init__.py