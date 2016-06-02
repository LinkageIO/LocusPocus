#!/usr/bin/env bash
# Script to setup camoco in an anaconda environment
# Written by Joe Jeffers and Rob Schaefer
# Email: jeffe174@umn.edu, schae234@umn.edu


function usage(){
cat <<EOF
    Usage: $0 [flags]

    Flags
    -----
    -h | --help
        print help message
    -b | --base
        Base installation directory for camoco (default: ~/.camoco).
EOF
exit 0
}

# Configurable variables
BASE=$HOME/.locuspocus

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

function red(){
    printf "${RED}$1${NC}\n"
}
function green(){
    printf "${GREEN}$1${NC}\n"
}

while [[ $# > 0 ]]
do
key="$1"
case $key in 
    -h|--help)
    usage
    shift
    ;;
    -b|--base)
    BASE=$2
    shift
    ;;
    *)  
        #unknown options
    ;;
esac
shift
done

NAME="locuspocus"

export CWD=$(pwd)
#===================================================
#----------Setup the build Environment--------------
#===================================================
echo "Setting up the build environment"
source $HOME/.bashrc
mkdir -p $BASE
mkdir -p $BASE/conda
cd $BASE

export LD_LIBRARY_PATH=$BASE/lib:$LD_LIBRARY_PATH
export PATH=$PATH:$BASE/bin:$BASE/conda/bin

#===================================================
#----------------Install conda ---------------------
#===================================================
# if ! hash conda 2>/dev/null
if [ ! -e $BASE/conda/bin/conda ]
then
    cd $BASE
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -f -p $BASE/conda
	rm -f miniconda.sh
else
    green "Conda Already Installed" 
fi

#===================================================
#----------Build the Conda Environment--------------
#===================================================
if [ ! -d $BASE/conda/envs/camoco ]
then
    echo "Making the conda virtual environment named $NAME in $BASE"
    conda remove -y --name $NAME --all
    conda config --add envs_dirs $BASE/conda/envs
    conda create -y -n $NAME --no-update-deps python=3.5 setuptools pip
else
    green 'conda already installed'
fi

#===================================================
#----------Activate the Conda Environment-----------
#===================================================
green "activating $NAME"
source $BASE/conda/bin/activate $NAME
green 'checking python'
which python

#==================================================
#-----------------Install Camoco-------------------
#=================================================
echo "Installing LocusPocus"
cd $CWD
python setup.py install
source deactivate 

#===================================================
#------------Update the bashrc----------------------
#===================================================
if [ $(grep $BASE/conda/bin ~/.bashrc | wc -l) -eq 0 ]
then 
    red '-----------------------------------------------'
    red "Update your $HOME/.bashrc:"
    echo "export LD_LIBRARY_PATH=$BASE/lib:\$LD_LIBRARY_PATH"
    echo "export PATH=$BASE/bin:$BASE/conda/bin:\$PATH"
    red '-----------------------------------------------'
fi

#===================================================
#-------------Use Instructions----------------------
#===================================================
green '==============================================='
echo "All done, to use locuspocus, first restart the shell."
echo "Then type: source activate $NAME"
echo "When finished, type: source deactivate"
echo ""
echo ""
echo "Please report all bugs to: https://github.com/schae234/LocusPocus/issues"
green '==============================================='
