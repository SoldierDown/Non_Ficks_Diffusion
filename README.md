# MPM Snow Simulation

## Compilation
The following instruction is tested on Ubuntu 18.04 LTS with g++7.5

### Dependencies
#### Install Nova library and OpenGL viewer
The Nova library can be downloaded from GitHub using the following command:
```
git clone https://github.com/OrionQuest/Nova.git
```
Example projects that use the Nova core library live in a separate repository, and can be downloaded using the following commands:
```
cd Nova/Projects
git clone https://github.com/OrionQuest/Nova_Examples.git
```
Nova depends on several libraries such as `GLM`,`FreeType`,`GLFW`,`GLEW`,`Boost`,`Assimp`, etc. Run the following command to install all the dependencies:
```
sudo apt-get install libboost1.62-dev libboost-program-options1.62-dev libboost-filesystem1.62-dev libboost-regex1.62-dev libglfw3-dev libglew-dev libglm-dev libassimp-dev libfreetype6-dev
```
An alternative:
```
sudo apt-get install libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-regex-dev libglfw3-dev libglew-dev libglm-dev libassimp-dev libfreetype6-dev
```
The build system uses cmake (version 3.0 or higher). We recommend the use of the graphical version ccmake for easy configuration of the environment variables. The best way to install the latest version of cmake and ccmake is to first run the following command:
```
sudo apt-get install cmake
sudo apt-get install cmake-qt-gui
sudo apt-get install libncurses5-dev
```
#### Clone Project
Nova/Projects. The partial directory should look like this:
```
cd Nova/Projects
git clone https://github.com/SoldierDown/Non_Ficks_Diffusion.git non_ficks_diffusion
```
The directory should look like this:

```
.
+-- Projects 
|   +-- non_ficks_diffusion
|   +-- Nova_Examples
|   +-- CMakeLists.txt
...
```
### Build
Run the following commands in order:
```
cd Nova
mkdir build
cd build
ccmake ..
```
Set ``CMAKE_BUILD_TYPE`` to ``Release``, and turn ``OFF`` the following flags: ``ENABLE_OPENIMAGEIO_PLUGIN``, ``ENABLE_SOIL_PLUGIN``, 
``ENABLE_VARIATIONAL_FLUIDS``, ``USE_C11_REGEX``, ``USE_DOUBLES``.Turn all the other flags ``ON``. 
Press ``c`` to configure, and then ``g`` to generate the ``Makefile``. Finally, run the following command:
```
make -j 8
```
### Configure OpenGL
Run the following commands in order:
```
cd build/
ln -s ../Projects/Nova_Examples/opengl/example/nova.conf .
ln -s ../Projects/Nova_Examples/opengl/example/fonts/ .
ln -s ../Projects/Nova_Examples/opengl/plugins/Grid/shaders/ .
```
Choose the correct plugin for visualization in nova.conf in *Nova/build/*

```
...
#Plugin=libplugin_Embedded_Deformables
Plugin=libplugin_MPM
...
```
Here are some of the command options for both plugins.
- ALT-P: planar camera
- P: play-pause
- S: next frame
- SHIFT-S: previous frame
- R: reset frame
- F2: show/hide background grid
- Scroll: zoom in/out
- Left-mouse: rotate
- Right-mouse: translate

### Run
Go to the build directory (*Nova/build/*), run the following command
```
./bin/mpm_2d -test_number 18 -size 64 64 -threads 8 -last_frame 300
```
### View Result
```
./bin/opengl mpm_2d_18
```
