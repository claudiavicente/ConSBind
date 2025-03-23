#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="ConSBind",
    version="1.0.0",
    author="Noelia Gil, Xavier Vílchez, Clàudia Vicente",
    description="Consensus Structural Binding site predictor",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/claudiavicente/ConSBind",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.19.0",
        "scipy>=1.5.0",
        "biopython>=1.78",
        "scikit-learn>=0.23.0",
    ],
    entry_points={
        "console_scripts": [
            "consbind=main:main",
        ],
    },
) 