#! /bin/bash

# Script base directory.
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Generate the binaries
SCRIPT=$(cat <<-END
echo "building danton with \$(gcc --version | head -n 1)"
cd /work
make PREFIX=AppDir
chown --recursive $(id -u):$(id -g) AppDir build
END
)

mkdir -p AppDir
make clean
docker run --mount type=bind,source=$(pwd),target=/work                        \
                   quay.io/pypa/manylinux1_x86_64 /bin/bash -c "${SCRIPT}"

# Copy extra application data
mkdir -p AppDir/share
cp -rL share/danton AppDir/share
cp -rL include AppDir

# Force the generation of the materials dump
./AppDir/bin/danton examples/cards/backward-tau-decays.json > /dev/null
rm -f steps.json

# Bundle AppImage meta and AppRun
cp ${basedir}/danton.desktop AppDir
cp ${basedir}/danton.png AppDir
cp ${basedir}/danton.appdata.xml AppDir
cp ${basedir}/AppRun AppDir

# Build the AppImage
appimagetool=appimagetool-$(arch).AppImage
[ -e ${appimagetool} ] || wget https://github.com/AppImage/AppImageKit/releases/download/continuous/${appimagetool} && chmod u+x ${appimagetool}
./${appimagetool} AppDir
