from __future__ import with_statement, print_function, absolute_import

from setuptools import setup, find_packages, Extension
from distutils.version import LooseVersion

import os
from glob import glob
from os.path import join

from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
	def finalize_options(self):
		_build_ext.finalize_options(self)
		__builtins__.__NUMPY_SETUP__ = False
		import numpy
		self.include_dirs.append(numpy.get_include())

dtw_c_src = join("src")
dtw_c_sources = glob(join(dtw_c_src, "*.c"))

ext_modules = [
	Extension(
		name="dtw_c.dtw_c",
		include_dirs=[dtw_c_src],
		sources=[join("dtw_c","dtw_c.pyx")] + dtw_c_sources,
		language="c")]

setup(
	name="dtw_c",
	ext_modules=ext_modules,
	cmdclass={'build_ext': build_ext},
	version='0.0.9.dev1',
	packages=find_packages(),
	setup_requires=[
		'numpy',
		'cython',
	],
	install_requires=[
		'numpy',
		'cython',
	],
	extras_require={
		'test': ['nose'],
		'sdist': ['numpy', 'cython'],
	},
	author="Patrick Lumban Tobing",
	author_email="patrickltobing@gmail.com",
	url="https://github.com/patrickltobing/dtw_c",
	description="Dynamic Time Warping (DTW) implementation in C for Python",
	keywords=['dynamic time warping','dtw'],
	classifiers=[],
)
