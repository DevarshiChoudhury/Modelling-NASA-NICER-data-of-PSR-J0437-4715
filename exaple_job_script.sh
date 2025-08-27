#!/bin/bash
#SBATCH -N 30
#SBATCH --tasks-per-node=128
#SBATCH -t 1-00:00:00
#SBATCH -p thin
#SBATCH --job-name=CST-PDT_Agn_init
#SBATCH --mail-user=d.choudhury@uva.nl
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load 2021
module load intel/2021a
module load Python/2.7.18-GCCcore-10.3.0-bare

cp -r $HOME/NICER_analyses/J0437/CST_PDT/ "$TMPDIR"

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_PRELOAD=/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_def.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_avx2.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_core.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_lp64.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_thread.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so
export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python _auto_modules/main_3c50Agn.py @config_Agn.ini --multinest > out_cst_pdt_3c50Agn 2> errcst_pdt_3c50Agn

cp out_cst_pdt_3c50Agn err_cst_pdt_3c50Agn nlive4000* $HOME/NICER_analyses/J0437/CST_PDT/CST_PDT_outputs/3c50Agn/init/
#end of job file
