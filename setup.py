#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

setup(
    name = 'locuspocus',
    version = '0.1.0',
    packages = find_packages(),
    scripts = [],

    install_requires = [
    ],

    package_data = {
    },
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Genetic coordinates so easy, it seems like MAGIC!',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/LocusPocus',


)
