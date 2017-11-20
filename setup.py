#!/usr/bin/env python
req = ['python-dateutil','nose','pytz','numpy','scipy','xarray',
       'msise00','gridaurora','glowaurora']
# %%
from setuptools import setup,find_packages

setup(name='reesaurora',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D',
      description='Model of Earth ionosphere.',
      version='1.0.1',
      url = 'https://github.com/scivision/reesaurora',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python',],
      install_requires=req,
      python_requires='>=3.6',
      extras_require={'plot':['matplotlib','seaborn'],},
	  )
