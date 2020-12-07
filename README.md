# UMass-ZC
Zero crossing algorithm for NMR signal analysis. 

# Compilation and Installation

## Configuration and Building   
Create a `build` and `install` directory in parallel to the main source directory.  Then change to 
the build directory.  From here, use CMake to configure the installation: 
`cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation /path/to/source`.  Alternately, this is encapsulated in 
the `config.sh` script in the source directory.  Then do: `make -jN && make install`, where N is the 
number of cores in your machine.  

## Configuring Your Environment
To be able to use the library in your application, you have to set up your environment with specific 
variables to successfully link this library to your application.  These will be defined soon.  
