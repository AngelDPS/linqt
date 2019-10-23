
#! /bin/bash SIMULATION PARAMETER
GPU=0
Nx=20
Ny=$Nx
M=40
W=0.0;
RUNNAME=MomentConvNx${Nx}Ny${Ny}W${W}
RUNDIR=$(pwd)
R=10;
ID=$RANDOM
if [ $GPU == 0 ]; then
	TGPU=""
else
 if [ $GPU == 1 ]; then 
	TGPU="2"
else
	echo "There is no computer with more than 3 gpu in our cluster"
	exit
 fi
fi


#LOCAL SCRIPT
cat > $RUNDIR/scripts/local/RunLocal$ID.sh << !
#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
source /opt/intel/bin/compilervars.sh intel64
bash $RUNDIR/scripts/FullCompile.sh GrapheneMagneticField
#bash \$ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version ${BASH_VERSION}..."
PCID=\`uname -n\`;
COMMAND=$RUNDIR/bin/GrapheneMagneticField
JOBID=\`echo \$1 | sed 's/[^0-9]*//g'\`;
SEED=\$RANDOM

#for M in 50 61 75 92 113 139 171 210 257 316 387 475 583 716 878 1078
for M in 3000
 
	do
        GPU=$GPU
			 #Nx  Ny  M   W R   MachineName, GPU, seed
	\${COMMAND}\$PCID $Nx $Ny \$M $W 1 $R Magnetic \$GPU \$SEED

 	done
!
chmod +x $RUNDIR/scripts/local/RunLocal$ID.sh

#TORQUE SCRIPT
cat > $RUNDIR/scripts/TorqueRun.sh << !
#PBS -l nodes=1:ppn=1:GPU$TGPU
##PBS -l nodes=1:ppn=1:main
#PBS -N $RUNNAME
#PBS -m abe
#PBS -M jgarcia@if.ufrj.br

# Bookeeping
echo    "Executando no no  : \$HOSTNAME"
echo -n "Data              : "
date
echo    "Job ID            : \$PBS_JOBID"
echo -n "Diretorio         : "
cd $RUNDIR
pwd
$RUNDIR/scripts/local/RunLocal$ID.sh \$PBS_JOBID
rm $RUNDIR/scripts/local/RunLocal$ID.sh 
!

qsub scripts/TorqueRun.sh




