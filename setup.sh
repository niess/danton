#!/bin/bash

# Script root directory.
danton_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set the path for dynamic libraries.
lib_dir=$danton_dir/lib
[[ "$LD_LIBRARY_PATH" =~ "${lib_dir}" ]] || export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH

# Set the materials.
export PUMAS_MDF=$danton_dir/share/materials/materials.xml
export PUMAS_DEDX=$danton_dir/share/materials/dedx

# Set the PATH for binaries.
bin_dir=$danton_dir/bin
[[ "$PATH" =~ "${bin_dir}" ]] || export PATH=${bin_dir}:$PATH

# Set the PYTHONPATH For python module.
python_dir=$danton_dir/lib/python
[[ "$PYTHONPATH" =~ "${python_dir}" ]] || export PYTHONPATH=${python_dir}:$PYTHONPATH
