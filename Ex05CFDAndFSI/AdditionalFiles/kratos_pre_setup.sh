# appending to VALUE -> $VALUE

# pointing to the Kratos Libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos/libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos

echo "LD_LIBRARY_PATH set to"
echo $LD_LIBRARY_PATH

# pointing to the Kratos folder
export PYTHONHOME=$PYTHONHOME:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos

echo "PYTHONHOME set to"
echo $PYTHONHOME

# pointing to the Kratos and python35.zip folders
export PYTHONPATH=$PYTHONPATH:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos/python34.zip
export PYTHONPATH=$PYTHONPATH:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos/

echo "PYTHONPATH set to"
echo $PYTHONPATH

export PATH=$PATH:/home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos

echo "PATH set to"
echo $PATH

echo "###   Environment for KratosMultiphysics Precompiled 711 started   ###"

# function instead of an alias
runkratos() {
    /home/kodakkal/Software/gid14.1.7d-x64/problemtypes/kratos.gid/exec/Kratos/runkratos "$@"
}

