from __future__ import print_function
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import glob
import sys
PY2 = sys.version_info[0] == 2
PYX6 = sys.version_info[1] <= 6

install_requires = [
    'matplotlib',
    'scipy',
    'numpy',
    'astropy'
    ]

if PY2 and PYX6:
    install_requires += ['unittest2']

print(install_requires)
setup(name='srttools',
      version='0.1',
      description="SRT Single dish tools",
      packages=['srttools', 'srttools.core'],
      package_data={'': ['README.md'], 'TEST_DATASET': ['*']},
      include_package_data=True,
      author='OAC high-energygroup',
      author_email="blablabla@bla.net",
      license='3-clause BSD',
      url='https://bitbucket.org/srttools/srt-single-dish-tools',
      keywords='Radio astronomy maps',
      scripts=glob.glob('scripts/*'),
      platforms='all',
      classifiers=[
          'Intended Audience :: Science/Research, Education, Developers',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Astronomy'
          ],
      install_requires=install_requires
      )
