#!/bin/sh
#PBS -N rserver
#PBS -o PBS-r.log
#PBS -e PBS-r.log
#PBS -q cpu
#PBS -l nodes=1:ppn=8
#PBS -l walltime=7200:00:00
#PBS -V
#PBS -S /bin/bash

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=192.168.203.202
timestamp=$(date | awk '{print $4}')
# print tunneling instructions
echo -e "
MacOS or linux terminal command on your local computer to create your ssh tunnel
First, Create a Public/Private Key Pair with empty passphrase on your local Mac and copy the string in id_rsa.pub to server .ssh/authorized_keys

nohup ssh -N -L ${port}:${node}:${port} ${user}@${cluster} &
" >> ~/Rserver.log
# DON'T USE ADDRESS BELOW.
# DO USE TOKEN BELOW

# Configure modules to load
#module load Singularity
#module load R/3.5
#export PATH=/home/${user}/anaconda3/envs/sc/bin:$PATH

# Work Dir
cd ~
#printf "rsession-which-r=${rpath}" > rstudiodb/rserver.conf
# Start jupyter
PASSWORD='qwerty1405' singularity exec --bind /home,rstudiodb/run:/run,rstudiodb/var-lib-rstudio-server:/var/lib/rstudio-server,rstudiodb/database.conf:/etc/rstudio/database.conf,rstudiodb/rserver.conf:/etc/rstudio/rserver.conf,rstudiodb/tmp:/tmp rstudiodb/tidyverse_latest.sif rserver --auth-none=0  --auth-pam-helper-path=pam-helper --auth-stay-signed-in-days=1000 --auth-timeout-minutes=0 --server-user=${user} --www-port=${port} 2>> ~/Rserver.log
