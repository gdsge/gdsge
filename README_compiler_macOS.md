# Instruction for setting up the Intel C++ compiler on a macOS

- Install Xcode from the Mac App Store

- Acquire the Intel C++ Compiler for macOS from the official website: [Intel® oneAPI standalone component installation files](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#dpcpp-cpp)

- Set the compiler environment by typing in the macOS terminal

  ```bash
  source /opt/intel/oneapi/setvars.sh
  ```

  or exporting the path that contains the compiler executable

  ```bash
  export PATH=/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/:$PATH
  ```

  Note that you need to modify the path-to-intel-compiler accordingly

- *In the same macOS terminal*, run MATLAB to make sure it can find the intel compiler