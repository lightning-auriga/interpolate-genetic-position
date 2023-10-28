# interpolate-genetic-position

## Brief Summary

Use reference genetic map to interpolate genetic position for a query set of variants

## Overview

This README is an automated stub generated from a `cookiecutter` template.
Documentation below reflects the state of the templated project immediately
after creation and may not reflect the current state of the project after
development updates.

## Requirements

  - g++ >= 8.2.0
  - automake/autoconf
  - make >= 4.2
  - git >= 2.28.0
  - nodejs (for commitizen)
  - pre-commit
  - associated linting tools for C++: cppcheck, clang-format
  - [boost headers](https://www.boost.org)
  - [boost program_options](https://www.boost.org/doc/libs/1_75_0/doc/html/program_options.html)
  - [boost filesystem/system](https://www.boost.org/doc/libs/1_75_0/libs/filesystem/doc/index.htm)
  - [boost iostreams](https://www.boost.org/doc/libs/1_74_0/libs/iostreams/doc/index.html)

## Build

By default, a build process involving a [conda/mamba](https://mamba.readthedocs.io/en/latest/installation.html) environment is supported.

  - if you wish to use `conda` and it's not currently available, you can install it with the instructions [here](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)
  - navigate into your project directory (interpolate-genetic-position)
  - create the `conda` environment for installation as follows:
  
     `mamba env create -f environment.yaml`
  - activate the conda environment:
  
     `mamba activate interpolate-genetic-position-env`
  - (one time only per environment) install `commitizen`:
  
     `npm install -g commitizen cz-conventional-changelog`
  - (one time only per environment) install `pre-commit` linters:
  
     `pre-commit install`

  - update (create) the necessary `configure` scripts with `autoreconf`:
  
     `autoreconf --force --install`
	 
     - note that this can also be run with `./generate.bash` inside the repo
  - run `configure`:
  
	 `CC=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc CXX=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-g++ ./configure --with-boost=${CONDA_PREFIX} --with-boost-libdir=${CONDA_PREFIX}/lib --with-yaml-cpp=${CONDA_PREFIX}`

	 - if you are planning on installing software to a local directory, run instead `./configure --prefix=/install/dir [...]`
  - run `make -j{ncores}`

  - if desired, run `make install`. if permissions issues are reported, see above for reconfiguring with `./configure --prefix`.
  
## Usage

By default, the final compiled program can be run with

`./interpolate-genetic-position.out`

## Version History

28 10 2023: project generated from cookiecutter template
