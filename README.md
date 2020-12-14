# UMass-ZC
A C++ library containing a zero crossing counting algorithm to extract the frequency of an input NMR signal.

# System Requirements 
Your system needs to have the following: 
- C++11 or higher   
- CMake version 3.9 or higher
- GNU Scientific Library installed 

# Configuration and Building   
Create a `build` and `install` directory in parallel to the main source directory.  Then change to 
the build directory.  From here, use CMake to configure the installation: 
`cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation /path/to/source`.  Alternately, this is encapsulated in 
the `config.sh` script in the source directory.  Then do: `make -jN && make install`, where N is the 
number of cores in your machine.  

# Configuring Your Environment
To be able to use the library in your application, you have to set up your environment with specific 
variables to successfully link this library to your application.  In your `bashrc`, define: 

```
export UMZCSYS=/path/to/UMass-ZC/installation
export LD_LIBRARY_PATH=$UMZCSYS/lib:$LD_LIBRARY_PATH
```
