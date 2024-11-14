#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (C) 2023 Intel Corporation.  All rights reserved.

from setuptools import setup, find_packages

def read_requirements(file):
    with open(file) as f:
        return f.read().splitlines()

def read_file(file):
   with open(file) as f:
        return f.read()
    
long_description = read_file("README.md")
requirements = read_requirements("requirements.txt")

setup(
    name = 'EAS_angles_sphere',
    version = '1.0.,
    author = 'Clemence Azra',
    author_email = 'clemenceanastasia@gmail.com',
    url = 'https://github.com/ClemenceAzra/Master_scientific_work/',
    description = 'Python package for the MSU SPHERE experiment',
    long_description_content_type = "text/x-rst",  # If this causes a warning, upgrade your setuptools package
    long_description = long_description,
    license = "MIT license",
    packages = find_packages(exclude=["test"]),  # Don't include test directory in binary distribution
    install_requires = requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]  # Update these accordingly
)
