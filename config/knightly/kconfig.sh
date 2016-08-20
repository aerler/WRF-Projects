#!/bin/bash
# this is a sample configuration file for knightly.sh and associated scripts

# IMPORTANT: make sure all environment variables are exported where necessary,
#            usually that means all, except PYTHON, SCRIPTS, NICENESS 

# essential variables
export CODE_ROOT="${CODE_ROOT:-/home/data/Code/}" # code root (makes things easier)
export DATA_ROOT="${DATA_ROOT:-/data-3/}" # root folder of data repository
# N.B.: these two variables need to respect defaults, otherwise the variables can not be changed from the command line

# some sensible defaults for Linux systems
export GDAL_DATA='/usr/local/share/gdal' # for GDAL API
PYTHON='/home/data/Enthought/EPD/bin/python' # path to Python executable (do not export!)
# Python modules and other scripts
export PYTHONPATH="${CODE_ROOT}/Projects/src/:${CODE_ROOT}/GeoPy/src/:${CODE_ROOT}/WRF Tools/Python/" # required modules
SCRIPTS="${CODE_ROOT}/WRF Tools/Scripts/SyncPP/" # folder with sync and post-processing scripts
# WRF & CESM data directories
export WRFDATA="${DATA_ROOT}/WRF/" # local WRF data root
export CESMDATA="${DATA_ROOT}/CESM/" # local CESM data root
# general settings
export FIGURE_ROOT="${HOME}" # this variable is not used; only to prevent crash
NICENESS=${NICENESS:-10}
export STATIC='STATIC' # download static/const data

# location of YAML configuration files for Python scripts
PYYAML_WRFAVG="${WRFDATA}/wrfavg/wrfavg.yaml"
PYYAML_EXSTNS="${DATA_ROOT}/exstns.yaml"
PYYAML_REGRID="${DATA_ROOT}/regrid.yaml"
PYYAML_EXPORT="${DATA_ROOT}/export.yaml"
PYYAML_SHPAVG="${DATA_ROOT}/shpavg.yaml"


# Environmental variable used by sync-wrf are defined here
#export DATASETS='Unity GPCC NARR CFSR CRU PRISM PCIC EC WSC HGS' # list of observational datasets
export DATASETS='Unity GPCC NARR CFSR CRU PRISM PCIC EC WSC' # list of observational datasets (without HGS ASCII)

## user specific settings
export SSHMASTER="${SSHMASTER:-$HOME}" # for shared connection
# connection settings
if [[ "${HOST}" == 'komputer' ]]
  then
    # download from komputer instead of SciNet using sshfs connection
    export SSH="-o BatchMode=yes"
    export HOST='fskomputer' # defined in .ssh/config
    export WRFSRC="${WRFDATA}/wrfavg/" # archives with my own wrfavg files
    export SUBDIR='WesternCanada GreatLakes'
    export CESMSRC='/data/CESM/cesmavg/' # archives with my own wrfavg files
    export OBSSRC='/data/'
    export INVERT='INVERT' # invert name/folder order in source (i.e. like in target folder)
elif [[ "${HISPD}" == 'HISPD' ]]
  then
    # high-speed transfer: special identity/ssh key, batch mode, and connection sharing
    export SSH="-o BatchMode=yes -o ControlPath=${SSHMASTER}/hispd-master-%l-%r@%h:%p -o ControlMaster=auto -o ControlPersist=1"
    export HOST='datamover' # defined in .ssh/config
    export WRFSRC='/reserved1/p/peltier/aerler/'
    export SUBDIR='WesternCanada GreatLakes'
    export CESMSRC='/reserved1/p/peltier/aerler//CESM/archive/'
    export OBSSRC='/reserved1/p/peltier/aerler/Datasets/'
    export INVERT='FALSE' # source has name first then folder type (like on SciNet)
else
    # ssh settings for unattended nightly update: special identity/ssh key, batch mode, and connection sharing
    # ssh settings: special identity/ssh key, batch mode, and connection sharing  
    export SSH="-i /home/me/.ssh/rsync -o BatchMode=yes -o ControlPath=${SSHMASTER}/master-%l-%r@%h:%p -o ControlMaster=auto -o ControlPersist=1"
    export HOST='aerler@login.scinet.utoronto.ca'
    export WRFSRC='/reserved1/p/peltier/aerler/'
    export SUBDIR='WesternCanada GreatLakes'
    export CESMSRC='/reserved1/p/peltier/aerler/CESM/archive/'
    export OBSSRC='/reserved1/p/peltier/aerler/Datasets/'
    export INVERT='FALSE' # source has name first then folder type (like on SciNet)
fi # if high-speed
