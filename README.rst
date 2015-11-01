==========
ReesAurora
==========

Rees-Sergienko-Ivanov module of excitation rates, relevant to auroral optical emissions
inspired/based upon Gustavsson / Brandstrom et al `AIDA_Tools <https://github.com/scienceopen/AIDA-tools>`_

.. contents::

Installation
============
::

  git clone --depth 1 https://github.com/scienceopen/reesaurora
  conda install --file requirements.txt
  python setup.py develop

Example
==================
::

  python RunRees.py 2011-03-15T12:34:56Z -c 65 -148

-o    specify output file (HDF5)
-c    specify geographic lat,lon
