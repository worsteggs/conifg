#!/bin/sh
#PBS -N gpu
#PBS -o PBS.log
#PBS -e PBS.log
#PBS -q gpu
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
First, Create a Public/Private Key Pair with empty passphrase on your local Mac and copy th
e string in id_rsa.pub to server .ssh/authorized_keys

nohup ssh -N -L ${port}:${node}:${port} ${user}@${cluster} &

Windows MobaXterm info
Forwarded port: same as remote port
Remote server: ${node}
Remote port: ${port}
SSH server: ${cluster}
SSH login: $user
SSH port: 22
Use a Browser on your local machine to go to:
localhost:${port} (prefix w/ https:// if using password)
You may be prompted to provide a token (a string of 48 letters + digits).
Please find '?token=...' in the 'jupyter-lab-${timestamp}.log' file under your current directory,
and copy the '...' tokens.
" >> ~/Jupyter-gpu-${timestamp}.log
# DON'T USE ADDRESS BELOW.
# DO USE TOKEN BELOW

# Configure modules to load
#module load Singularity
#module load R/3.5
export PATH=/home/${user}/anaconda3/bin:$PATH

# Work Dir
mkdir ~/github/jlab/gpu-${timestamp}
cd ~/github/jlab/gpu-${timestamp}
ln -s ~/github github
# Start jupyter
jupyter-lab --no-browser --port=${port} --ip=${node} 2>> ~/Jupyter-gpu-${timestamp}.log
