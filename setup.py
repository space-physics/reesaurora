#!/usr/bin/env python3

from setuptools import setup

with open('README.rst','r') as f:
  	long_description = f.read()

setup(name='reesaurora',
      version='0.1',
	  description='Rees-Sergienko-Ivanov auroral model',
	  long_description=long_description,
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/reesaurora',
	  install_requires=['h5py','six','nose','pytz'],
      packages=['reesaurora'],
	  )

