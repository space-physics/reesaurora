[![image](https://zenodo.org/badge/36744637.svg)](https://zenodo.org/badge/latestdoi/36744637)
[![image](https://travis-ci.org/scivision/reesaurora.svg)](https://travis-ci.org/scivision/reesaurora)
[![image](https://coveralls.io/repos/scivision/reesaurora/badge.svg?branch=master)](https://coveralls.io/github/scivision/reesaurora?branch=master)
[![Maintainability](https://api.codeclimate.com/v1/badges/fae4ee1dfb20a766ebce/maintainability)](https://codeclimate.com/github/scivision/reesaurora/maintainability)
[![PyPi Download stats](http://pepy.tech/badge/reesaurora)](http://pepy.tech/project/reesaurora)


# Rees-Sergienko-Ivanov Auroral model

Rees-Sergienko-Ivanov model of excitation rates, relevant to auroral
optical emissions inspired/based upon Gustavsson / Brandstrom et al
[AIDA_Tools](https://github.com/scivision/AIDA-tools)

Model designed for **100 - 10,000 eV**, and is essentially a *parameter
fit* to more advanced models, making for convenient computation in this
energy range with the PCs of the early 1990s. Today, much more advanced
physics-based models are tractable on a PC.

Uses MSISE-00 to generate O, O~2~, N~2~ densities, and models outcome of
primary electron precipitation on this neutral background.

![volume production rate](tests/demo.png)


## Install

Requires a Fortran compiler, such as `gfortran`. 
Any Fortran compiler should work:

    pip install -e .

## Example

    ReesSergienkoIvanov -t 2011-03-15T12:34:56Z -c 65 -148

* `-o` specify output file (HDF5) 
* `-c` specify geographic lat,lon

