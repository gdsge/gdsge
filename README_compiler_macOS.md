# Instruction for setting up he g++ as the Mex compiler on macOS

- Why g++ but not Xcode/Intel Compiler? (1) Xcode does not officially support openmp. (2) Support from Intel is conditional on the open software license which may not be renewed. (3) Supporting g++ automatically supports MinGW under Windows.

- Install [homebrew](https://brew.sh/) in a macOS terminal

  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

- Install g++ version 8. Mex files compiled with newer version of g++ crash MATLAB for known compatibility issues.
  ```bash
  brew install gcc@8
  ```

- Launch MATLAB, change directory to "gdsge/plugins". Run
  ```matlab
  mex -setup -f g++_maci64.xml
  ```

  to  setup g++ as the mex compiler.