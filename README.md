# FVLib

**A Fortran library of advanced computational algorithms and numerical methods to solve partial differential equations within the finite volume philosophy.**

## Description

The project aims to bring high-fidelity and high-performance simulations of a wide range of physics and engineering problems in relevant industrial, environmental, and biomedical applications. In that regard, the development of the FVLib code is based on the following principles:

- *Modern object-oriented Fortran (2003/2008 standards)* - for better code reuse and organisation
- *High-scalability with multi-threading execution* - to take advantage of modern HPC systems
- *Unstructured meshes with general element shapes* - for complex geometries in real applications
- *Very high-order of convergence both in space and time* - for highly-accurate and efficient solutions

Currently, the following problems can be solved in the FVLib code:

- Convection-diffusion problems for heat and species transfer
- Conjugate heat transfer with solid/solid and solid/fluid interfaces
- Incompressible (non-)isothermal (non-)Newtonian fluid flows with the Euler/Stokes/Oseen/Navier-Stokes equations

## Building

### Requirements

The FVLib code can be installed on any Linux and other Unix-based operating systems and requires a modern compiler that implements Fortran 2003 standard or above. Additionally, the compiler should provide OpenMP version 3 or above for multi-threading execution. The following compilers are currently supported:

- GNU Compiler Collection (GCC), version 5.0 or above.
- Fujitsu Technical Computing Suite (TCSDS), version 1.2.26 or above.
- Arm Compiler for Linux (ACFL), version 20.3 or above.

Additionally, for a parallel installation (recommended), a build automation tool is also required. The following make utilities are supported:

- GNU Make (Make), version 3.8 or above.

> **NOTE**: The choice of the compiler version strongly determines the build time, as well as the processor of the computer where the FVLib code is being installed. Newer versions require significantly more time to compile, but more aggressive optimisations are usually performed. It is recommended to choose the compiler version according to the purpose: for developing and testing in a personal computer, GCC version 5.0 or 6.0 is recommended for short build times; for intensive testing in high-performance computing systems, GCC version 7.0 or above is preferable.

### Quick installation

Provided that a compatible GCC compiler and make utility are available on the computer, and that the default configuration is appropriate (otherwise, follow the step-by-step installation in the detailed instructions), the following commands will install the FVLib code in the current directory of a bash session:
```
git clone https://github.com/ricardodpcosta/FVLib.git
cd FVLib
source etc/bashrc.sh
fvl-make --make
cd build
make all -s -j
cd ..
echo "alias fvlib='source $(pwd)/etc/bashrc.sh' >> ~/.bashrc"
```

### Detailed instructions

Open a bash session on the computer and follow the instructions below.

**Step 1.** Make sure the compiler and the other required utilities are correctly configured on the computer and check whether the installed version is compatible as provided in the requirements section:

- GCC compiler:
```
gfortran --version
```

- Fujitsu compiler:
```
frt --version
```

- Arm compiler:
```
armflang --version
```

- Make:
```
make --version
```

**Step 2.** Download the FVLib code package to any directory of your choice. For instance, `$HOME` can be used to install only for the current user, whereas `/opt` or `/share` is preferable to install for all the users. Uncompress the FVLib code package if necessary and enter the folder:

- Using Git:
```
git clone https://gitlab.com/fvlib/fvlib-release.git
cd fvlib-release
```

- Manual download in tar.gz format:
```
tar -xzvf fvlib-release.tar.gz
cd fvlib-release
```

- Manual download in .zip format:
```
unzip fvlib-release.zip
cd fvlib-release
```

**Step 3.** The default configuration of the FVLib code is to compile with the GCC compiler (gfortran) on an x86_64 architecture. For other combinations of compiler-architecture, several configurations files are provided inside the `etc/` folder. Rename the most appropriate configuration file to `compiler.sh`:

- GNU compiler in AArch64 architecture:
```
cp etc/gfortran-aarch64.sh etc/compiler.sh
```

- Fujitsu compiler in AArch64 architecture:
```
cp etc/fuji-aarch64.sh etc/compiler.sh
```

- Arm compiler in AArch64 architecture:
```
cp etc/arm-aarch64.sh etc/compiler.sh
```

**Step 4.** Source the global configurations file:
```
source etc/bashrc.sh
```

**Step 5.** The FVLib code provides a library management utility to compile, update, clean and create dependencies of the library and applications. Make sure the global configurations file is correctly sourced and the library management utility is available on the computer:
```
fvl-make --version
```

The installation of the FVLib code can be performed in serial or parallel. Since the installation might take several hours to finish depending on the compiler and the processor, the second option is usually preferable.

- For a serial installation of the library and applications run:
```
fvl-make --all
```

- For a parallel installation of the library and applications, generate the Makefile with targets and dependencies for the make utility:
```
fvl-make --make
```

A file is created inside the `build/` folder, from where the make utility can be run:
```
cd build
make all -s -j
```

> **NOTE**: The compilation of the FVLib code and applications can take several hours to finish depending on the compiler and the processor. In the serial installation, a percentage is shown on the terminal to indicate the progress. In case the installation needs to be stopped and later resumed, simply run the `fvl-make --all` (serial installation) or `make all -s -j` (parallel installation) commands again and the installation will proceed from the previous state.

> **NOTE**: Library source and applications can also be installed in several steps. Run the `fvl-make --all` (serial installation) or `make all -s -j` (parallel installation) commands replacing the `all` option with `library` to install the library, then with `tools` to install the tools, and finally with `solvers` to install the solvers.

**Step 6.** An alias (shortcut) to source the global configurations file can be created in the `~/.bashrc` file (or the `~/.profile` file):
```
echo "alias fvlib-release='source $HOME/etc/bashrc.sh' >> ~/.bashrc"
```

Replace `$HOME` with the installation path of the FVLib code if necessary. The alias will be available after sourcing the edited file or after opening a new Bash session (command windows or tab). Sourcing the global configurations file is required for development and for running the applications (including pre- and post-processing). With the alias created, source the global configurations file automatically whenever needed with the command:
```
fvlib-release
```

### Debugging mode

The default configuration installs the FVLib code in production mode, meaning that the compiled library and applications are optimised for performance rather than for debugging. On the other side, for developing and testing, the debugging mode is recommended to enable verification routines and specific compiler debugging options. For that purpose, all the previous steps involving the `fvl-make` command should be run with the `--debugging` option before any other option.

> **NOTE**: Only either production or debugging mode should be used for the current installation. In order to have both modes available on the computer, create two installations in separate folders. For instance, `fvlib-release/` for the production mode and `fvlib-release-dbg/` for the debugging mode.

Provided that a compatible GCC compiler and make utility are available on the computer, and that the default configuration is appropriate, the following commands will install the FVLib code in debugging mode in the current directory of the bash session:
```
git clone https://gitlab.com/fvlib/fvlib-release.git
mv fvlib-release fvlib-release-dbg
cd fvlib-release-dbg
source etc/bashrc.sh
fvl-make --debugging --make
cd build
make all -s -j
cd ..
echo "alias fvlib-release-dbg='source $(pwd)/etc/bashrc.sh' >> ~/.bashrc"
```

## Third-party software

The FVLib code has minimal dependencies on third-party software. At the current version of the release, only [Gmsh](https://gmsh.info/) is required to be installed. Version 4.8.4 or above is recommended for compatibility.

Make sure Gmsh is available on the computer:
```
gmsh --version
```

> **NOTE**: In Ubuntu and other APT-based Linux distributions, Gmsh can be easily installed with `sudo apt-get install gmsh`. For other distributions, follow the instructions provided on the [Gmsh](https://gmsh.info/#Download) webpage.

---

## Contributing

Contributions are encouraged in all forms, including the addition of functionalities, refinement of source code, enhancement of documentation, and resolution of issues.
For detailed guidelines, refer to [CONTRIBUTING.md](CONTRIBUTING.md).

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

📧 [Ricardo Costa](mailto:rcosta\@dep.uminho.pt)  
🌐 [Academic page](https://ricardodpcosta.github.io/)
