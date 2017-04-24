#!/usr/bin/env python
req = ['python-dateutil','nose','pytz','numpy','scipy','xarray','h5py','astropy','matplotlib','seaborn','pathlib2']
pipreq=['msise00','gridaurora','glowaurora']
# %%
import pip
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    pip.main(['install'] +req)
pip.main(['install']+pipreq)
# %%
from setuptools import setup  # enables develop

setup(name='reesaurora',
      packages=['reesaurora'],
      author='Michael Hirsch, Ph.D',
      description='Model of Earth ionosphere.',
      version='1.0.1',
      url = 'https://github.com/scivision/reesaurora',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',],
	  )
