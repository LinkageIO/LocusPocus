#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext

import io
import os
import re
import numpy

from setuptools.command.develop import develop
from setuptools.command.install import install
from subprocess import check_call

class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        print('Running post-installation setup')
        check_call('''\
	    pip install -r requirements.txt
        '''.split())
        develop.run(self)

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        print('Running post-installation setup')
        check_call('''\
        pip install -r requirements.txt
        '''.split())
        install.run(self)

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

locusdist = Extension(
    'locuspocus.LocusDist',
    sources=['locuspocus/LocusDist.pyx'],
#    extra_compile_args=['-ffast-math'],
    include_dirs=[numpy.get_include()]
)

setup(
    name = 'locuspocus',
    version = find_version('locuspocus','__init__.py'),
    packages = find_packages(),
    scripts = [],
    ext_modules = [locusdist],
    cmdclass = {
        'build_ext' : build_ext,
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    package_data = {
        '':['*.cyx']
    },
    install_requires = [
        'cython>=0.16',
        'numpy>=1.9.1'
    ],
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'rob@linkage.io',
    description = 'Genetic coordinates so easy, it seems like MAGIC!',
    license = "Copyright Linkage Analytics 2017. MIT License",
    url = 'linkage.io',
)
