#!/usr/bin/env python3

from setuptools import setup
import subprocess

with open('README.rst','r') as f:
  	long_description = f.read()

try:
    subprocess.run(['conda','install','--yes','--quiet','--file','requirements.txt'])
except Exception as e:
    print('you will need to install packages in requirements.txt  {}'.format(e))

setup(name='reesaurora',
      version='0.1',
	  description='Rees-Sergienko-Ivanov auroral model',
	  long_description=long_description,
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/reesaurora',
   dependency_links = ['https://github.com/scienceopen/msise00/tarball/master#egg=msise00',
                      'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora'],
	  install_requires=['msise00','gridaurora'],
      packages=['reesaurora'],
	  )
