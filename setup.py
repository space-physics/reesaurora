#!/usr/bin/env python
from setuptools import setup
import subprocess

try:
    subprocess.call(['conda','install','--quiet','--file','requirements.txt'])
except Exception as e:
    pass

setup(name='reesaurora',
	  description='Rees-Sergienko-Ivanov auroral model',
	  url='https://github.com/scienceopen/reesaurora',
   dependency_links = ['https://github.com/scienceopen/msise00/tarball/master#egg=msise00',
                      'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora'],
	  install_requires=['msise00','gridaurora',
                    'pathlib2'],
      packages=['reesaurora'],
	  )
