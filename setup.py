#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

import io
import os
import re

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
)

setup(
    name = 'locuspocus',
    version = find_version('locuspocus','__init__.py'),

    description = 'Genetic coordinates made easy!',
    url = 'http://linkage.io',
    author = 'Rob Schaefer',
    author_email = 'rob@linkage.io',
    license = "Copyright Linkage Analytics 2016. Available under the MIT License",

    classifiers=[
	# How mature is this project? Common values are
	#   3 - Alpha
	#   4 - Beta
	#   5 - Production/Stable
	'Development Status :: 4 - Beta',

	# Indicate who your project is intended for
	'Intended Audience :: Developers',
	'Topic :: Software Development :: Build Tools',

	# Pick your license as you wish (should match "license" above)
	 'License :: OSI Approved :: MIT License',

	# Specify the Python versions you support here. In particular, ensure
	# that you indicate whether you support Python 2, Python 3 or both.
	'Programming Language :: Python :: 3',
	'Programming Language :: Python :: 3.6',
    ],
    keywords='data genetics biology coordinates', 
    project_urls={
        'Documentation' : 'http://linkage.io',
        'Source' : 'https://github.com/LinkageIO/LocusPocus',
        'Tracker' : 'https://github.com/LinkageIO/LocusPocus/issues'
    },


    packages = find_packages(),
    scripts = [],
    ext_modules = [locusdist],
    cmdclass = {
    },
    package_data = {
        '':['*.cyx']
    },
    setup_requires = [
        # Setuptools 18.0 properly handles Cython extensions.
        'setuptools>=18.0',
        'cython',
    ],
    install_requires = [
        'minus80>=0.1.2',
        'Cython>=0.16.0',
        'numpy>=1.14.3',
        'scipy>=0.19.0'
    ],
    include_package_data=True,

)
