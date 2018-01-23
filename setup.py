#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import os

setup(
	cmdclass = {'build_ext': build_ext},
	ext_modules = [Extension("dtw",
							sources=["dtw.pyx", "dtw_sub.c"],
							include_dirs=[numpy.get_include()])],
)
