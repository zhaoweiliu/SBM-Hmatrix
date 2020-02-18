
# Singular Boundary Method with H-matrix

### Developed based on from a CMake C++ Project Template

Thank you for your interest in this project!

## Usage

### Prerequisites

You will need:

 * A modern C/C++ compiler
 * CMake 3.1+ installed (on a Mac, run `brew install cmake`)
 * If you prefer to code in a great IDE, I highly recommend [Jetbrains CLion](https://www.jetbrains.com/clion/). It is fully compatible with this project.

### Building The Project

#### Git Clone

First we need to check out the git repo:

```bash
❯ mkdir ~/workspace
❯ cd ~/workspace
❯ git clone \
    https://github.com/zhaoweiliu/SBM-Hmatrix.git \
    my-project
❯ cd my-project
❯ bash build-and-run
```

The output of this script is rather long and is shown [on this screenshot](doc/build-and-run.png).

The script `build-and-run` is a short-cut — you shouldn't really be using this script to build your project, but see how to do it properly below.

#### Project Structure

There are three empty folders: `lib`, `bin`, and `include`. Those are populated by `make install`.

The rest should be obvious: `src` is the sources, and `test` is where we put our unit tests.

Now we can build this project, and below we show three separate ways to do so.

#### Building Manually

```bash
❯ rm -rf build && mkdir build
❯ git submodule init && git submodule update
❯ cd build
❯ cmake ..
❯ make && make install
❯ cd ..
```

### Using it as a C++ Library

We build a static library that, given a simple fraction will return the integer result of the division, and the remainder.

We can use it from C++ like so:

```cpp
#include <iostream>
#include <division>

Fraction       f = Fraction{25, 7};
DivisionResult r = Division(f).divide();

std::cout << "Result of the division is " << r.division;
std::cout << "Remainder of the division is " << r.remainder;
```

## File Locations

 * `src/*` — C++ code that ultimately compiles into a library
 * `bin/`, `lib`, `include` are all empty directories, until the `make install` install the project artifacts there.
