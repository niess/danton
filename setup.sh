#!/bin/bash

# Script root directory.
danton_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set the path for dynamic libraries.
lib_dir=$danton_dir/ent/lib
[[ "$LD_LIBRARY_PATH" =~ "${lib_dir}" ]] || export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH
lib_dir=$danton_dir/pumas/lib
[[ "$LD_LIBRARY_PATH" =~ "${lib_dir}" ]] || export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH
lib_dir=$danton_dir/alouette/lib
[[ "$LD_LIBRARY_PATH" =~ "${lib_dir}" ]] || export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH

# Set the materials.
export PUMAS_MDF=$danton_dir/materials/materials.xml
export PUMAS_DEDX=$danton_dir/materials/dedx

# Set the PATH for binaries.
bin_dir=$danton_dir/bin
[[ "$PATH" =~ "${bin_dir}" ]] || export PATH=${bin_dir}:$PATH
