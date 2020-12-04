#! /bin/bash  
# a simple script to build 
#_______________________________________________________________________________
function cleanup {
   echo "Cleaning up previous build..."
   rm -r CMakeFiles *.cmake *.txt *.cxx *.pcm *.rootmap Makefile *.so
   echo "--> Done!"
}
#_______________________________________________________________________________
function configure {
   echo "Configuring..."
   cmake -DCMAKE_INSTALL_PREFIX=../install \
   ../UMass-ZC
   echo "--> Done!"
}
#_______________________________________________________________________________
function build_lib {
   echo "Building..."
   make && make install
   # for macOS
   cd ../install/lib && rm *.so && ln -s libUMassZC.dylib libUMassZC.so
   echo "--> Done!"
}
#_______________________________________________________________________________
# main 
cleanup
configure
build_lib
