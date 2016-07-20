.. image:: https://travis-ci.org/scienceopen/reesaurora.svg
    :target: https://travis-ci.org/scienceopen/reesaurora

.. image:: https://coveralls.io/repos/scienceopen/reesaurora/badge.svg?branch=master&service=github 
    :target: https://coveralls.io/github/scienceopen/reesaurora?branch=master

.. image:: https://codeclimate.com/github/scienceopen/reesaurora/badges/gpa.svg
   :target: https://codeclimate.com/github/scienceopen/reesaurora
   :alt: Code Climate    
==========
ReesAurora
==========

Rees-Sergienko-Ivanov model of excitation rates, relevant to auroral optical emissions
inspired/based upon Gustavsson / Brandstrom et al `AIDA_Tools <https://github.com/scienceopen/AIDA-tools>`_
 
Model designed for **100 - 10,000 eV**, and is essentially a *parameter fit* to more advanced
models, making for convenient computation in this energy range with the PCs of the early 1990s. Today, much more advanced physics-based models are tractable on a PC.

Uses MSISE-00 to generate O, O\ :sub:`2`, N\ :sub:`2` densities, and models outcome of primary electron precipitation on this neutral background. 

.. image:: test/demo.png
   :alt: volume production rate

.. contents::

Installation
============
::

  python setup.py develop

Example
==================
::

  python RunRees.py -t 2011-03-15T12:34:56Z -c 65 -148

-o    specify output file (HDF5)
-c    specify geographic lat,lon
