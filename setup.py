#!/usr/bin/env python
from setuptools import setup
try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)

setup(name='reesaurora',
   dependency_links = [
      'https://github.com/scienceopen/msise00/tarball/master#egg=msise00',
      'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora'
      'https://github.com/scienceopen/glowaurora/tarball/master#egg=glowaurora'
        ],
	  install_requires=['msise00','gridaurora'],
      extras_require = {'glowaurora':['glowaurora']},
      packages=['reesaurora'],
	  )
