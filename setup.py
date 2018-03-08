#!/usr/bin/env python
install_requires = ['python-dateutil','pytz','numpy','scipy','xarray',
       'msise00','gridaurora']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import setup,find_packages

setup(name='reesaurora',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D',
      description='Model of Earth ionosphere.',
      long_description=open('README.rst').read(),
      version='1.0.2',
      url = 'https://github.com/scivision/reesaurora',
      classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3.6',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      install_requires=install_requires,
      python_requires='>=3.6',
      extras_require={'plot':['matplotlib','seaborn'],
                      'tests':tests_require},
      tests_require=tests_require,
	  )
