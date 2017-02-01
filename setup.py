#!/usr/bin/env python
from setuptools import setup
try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)

setup(name='reesaurora',
      dependency_links = ['https://github.com/scienceopen/msise00/tarball/master#egg=msise00',
                      'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora'],
	  install_requires=['msise00','gridaurora'],
      packages=['reesaurora'],
	  )
