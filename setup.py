#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
setup(
  name = 'goodvibes',
  packages = ['goodvibes'],
  version = '2.0.1',
  description = 'Calculates quasi-harmonic free energies from Gaussian output files with temperature and haptic corrections',
  author = 'Robert Paton',
  author_email = 'robert.paton@chem.ox.ac.uk',
  url = 'https://github.com/bobbypaton/goodvibes',
  download_url = 'https://github.com/bobbypaton/GoodVibes/archive/2.0.1.zip',
  keywords = ['compchem', 'thermochemistry', 'gaussian', 'vibrational-entropies', 'temperature'],
  classifiers=[
    'Programming Language :: Python',
    'Natural Language :: English',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License:: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering :: Chemistry',
  ],
  install_requires=["numpy", ],
  python_requires='>=2.6',
  include_package_data=True,
  entry_points='''
    [console_scripts]
    goodvibes=goodvibes.GoodVibes:main
    '''
)
