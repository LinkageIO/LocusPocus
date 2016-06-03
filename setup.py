#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

import io
import os
import re

from setuptools import setup, find_packages, Extension

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


setup(
    name = 'locuspocus',
    version = find_version('locuspocus','__init__.py'),
    packages = find_packages(),
    scripts = [],
    install_requires = [],
    package_data = {},
    include_package_data=True,
    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Genetic coordinates so easy, it seems like MAGIC!',
    license = "MIT License",
    url = 'https://github.com/schae234/LocusPocus',


)
