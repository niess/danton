#!/bin/bash

# Script root directory.
danton_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set the environment for the submodules.
. $danton_dir/ent/setup.sh
. $danton_dir/pumas/setup.sh
. $danton_dir/tauola-c/setup.sh

# Set the materials.
export PUMAS_MDF=$danton_dir/materials/materials.xml
export PUMAS_DEDX=$danton_dir/materials/dedx

# Set the PATH for binaries.
bin_dir=$danton_dir/bin
[[ "$PATH" =~ "${bin_dir}" ]] || export PATH=${bin_dir}:$PATH
