# Instruction for setting up g++ as the Mex compiler on macOS Silicon

- Warning: need MacOS silicon processor (i.e., M1+)

- Install [homebrew](https://brew.sh/) in a macOS terminal

  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

- Install g++ version 15.
  ```bash
  brew install gcc@15
  ```

- Launch MATLAB, change directory to "gdsge/plugins". Run
  ```matlab
  mex -setup -f g++_maca64.xml
  ```

  to  setup g++ as the mex compiler.