#!/bin/bash
# this is a configuration file to synchronize my laptop with my SciNet/komputer archive (using sync-client.sh)
# Andre R. Erler, 30/05/2016, GPL v3

# essential variables
export CODE_ROOT="${CODE_ROOT:-${HOME}/Code/}" # code root (makes things easier)
export DATA_ROOT="${DATA_ROOT:-/data/}" # root folder of data repository
# N.B.: these two variables need to respect defaults, otherwise the variables can not be changed from the command line
SCRIPTS="${CODE_ROOT}/WRF Tools/Scripts/SyncPP/" # folder with sync and NCO scripts
NICENESS=${NICENESS:-10}

# Environment variables used by rsync scripts are defined here

# observational datasets
export DATASETS='Unity GPCC NARR CFSR CRU PRISM PCIC EC WSC HGS' # list of observational datasets
# WRF & CESM data directories
export WRFDATA="${DATA_ROOT}/WRF/" # local WRF data root
export SUBDIR='GreatLakes' # project subfolders
#export SUBDIR='WesternCanada GreatLakes' # project subfolders
export CESMDATA="${DATA_ROOT}/CESM/" # local CESM data root

# environment variables and settings for rsync scripts
WRFCLTSREX='NONE'; WRFCLTSFT='NONE'
WRFCLIMREX='*-*'; WRFCLIMFT='wrf*_grw2_clim_*.nc'
CESMCLIMFT='cesmatm_grw2_clim_*.nc cesmlnd_grw2_clim_*.nc'

# SSH settings for SciNet (HPC cluster)
# same as ssh settings for unattended nightly update: 
# special identity/ssh key, batch mode, and connection sharing
CLUSTER='aerler@login.scinet.utoronto.ca'
CLUSTERSSH="-i ${HOME}/.ssh/rsync -o BatchMode=yes -o ControlPath=${HOME}/master-%l-%r@%h:%p -o ControlMaster=auto -o ControlPersist=1"
WRFCLSTR='/reserved1/p/peltier/aerler/'
CESMCLSTR='/reserved1/p/peltier/aerler/CESM/archive/'
OBSCLSTR='/reserved1/p/peltier/aerler/Datasets/'
# SSH settings for komputer (workstation)
WORKSTN='fskomputer' # has to be defined in .ssh/config
WORKSTNSSH='-o BatchMode=yes'
WRFWRKSTN='/data/WRF/wrfavg/'
CESMWRKSTN='/data/CESM/cesmavg/'
# N.B.: the logins can also be SSH configurations which are defined in ${HOME}/.ssh/config


