#!/usr/bin/env python
import builtins
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from pybind11.setup_helpers import Pybind11Extension, build_ext
import numpy

class CustomBuildExt(build_ext):
    def finalize_options(self):
        super().finalize_options()
        # Prevent numpy from thinking it is in setup process
        builtins.__NUMPY_SETUP__ = False
        self.include_dirs.append(numpy.get_include())

# Define the pybind11 extension
cpp_extension = [
    Pybind11Extension(
        'pylandau',
        sources=['pyLandau/cpp/pylandau_src.cpp'],  # Path to the pybind11-based C++ file
        language="c++"
    ),
]

# Requirements
install_requires = ['pybind11>=2.6', 'numpy>=1.21']

setup(
    cmdclass={'build_ext': CustomBuildExt},
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
    package_data={'': ['README.*', 'VERSION'], 'docs': ['*'], 'examples': ['*']},
    ext_modules=cpp_extension,
    keywords=['Landau', 'Langau', 'PDF'],
    platforms='any'
)
