#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError as e:
    from distutils.core import setup

requirements = [
    'numpy>=1.14',
    'scipy',
    'astropy',
    'omegaconf>=2.2.2',
    'psutil',
    'matplotlib',
    'future-fstrings',
]
__casa_max__ = 3.8

if 3.6 <= float(sys.version[:3]) <= __casa_max__:
    requirements.append('casatasks')
else:
    print(f'You are running a python version ({sys.version[:3]}) for which modular casa is not avaliable.')

if float(sys.version[:3]) < 3.7:
    requirements.append('importlib_resources>=3.3.0')

PACKAGE_NAME = 'pyHIARD'
__version__ = '1.1.9'


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name=PACKAGE_NAME,
      version=__version__,
      description="Development Status :: 1 - Beta",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="P. Kamphuis",
      author_email="peterkamphuisastronomy@gmail.com",
      url="https://github.com/PeterKamphuis/pyHIARD",
      packages=[PACKAGE_NAME],
      python_requires='>=3.6',
      install_requires=requirements,
      include_package_data=True,
      # package_data - any binary or meta data files should go into MANIFEST.in
      scripts=["bin/" + j for j in os.listdir("bin")],
      license="GNU GPL v3",
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python :: 3",
          "Topic :: Scientific/Engineering :: Astronomy"
      ]
      )
