mkdir build
meson setup --wipe build
cd build
ninja
ninja python
cd ..
