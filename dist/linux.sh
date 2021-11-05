#! /bin/bash

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

# Bundle AppImage meta and entry point
cat << END > AppDir/danton.desktop
[Desktop Entry]
Type=Application
Name=danton
Exec=danton
Comment=Simulate the coupled transport of ultra high energy taus and neutrinos through the Earth, by Monte-Carlo
Icon=danton
Categories=Science;
Terminal=true
END

cat << END > AppDir/AppRun
#! /bin/bash
DANTON_PREFIX=\$APPDIR \$APPDIR/bin/danton \$*
END
chmod u+x AppDir/AppRun
