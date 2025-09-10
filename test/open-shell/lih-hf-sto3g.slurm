#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=faculty
##SBATCH --partition=jochena
##SBATCH --qos=jochena
##SBATCH --account=jochena
##SBATCH --mem=55000
#SBATCH --job-name=lih-hf-sto3g
#SBATCH --output=lih-hf-sto3g.log
#SBATCH --mail-type=END
#SBATCH --constraint=CPU-E5-2650v4

#export LD_LIBRARY_PATH=/projects/academic/jochena/software/orca/orca_5_0_4_linux_x86-64_shared_op    enmpi411:${LD_LIBRARY_PATH}

export PATH=/projects/academic/jochena/dumitruc/ompi-4.8/build/bin:$PATH
export LD_LIBRARY_PATH=/projects/academic/jochena/dumitruc/ompi-4.8/build/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/projects/academic/jochena/dumitruc/orca:${LD_LIBRARY_PATH}

cd $SLURMTMPDIR
pwd
export HomeDir=$SLURM_SUBMIT_DIR
export WorkDir=$SLURMTMPDIR

echo $PATH
echo $LD_LIBRARY_PATH

export | grep SLURM

#cd $WorkDir
cp $HomeDir/lih-hf-sto3g.inp .
#cp $HomeDir/lih-hf-sto3g.gbw .

/projects/academic/jochena/dumitruc/orca/orca lih-hf-sto3g.inp >  $HomeDir/lih-hf-sto3g.out

# mv test_yu.gbw $HomeDir/test_yu.gbw
 mv *qro $HomeDir
 mv *gbw $HomeDir
 mv *uno $HomeDir
 mv *loc $HomeDir
 mv *.47 $HomeDir
 mv *_property.txt $HomeDir
exit

