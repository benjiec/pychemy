#!/usr/bin/env python

from setuptools import setup

setup(name='pychemy',
      version='0.2',
      description='Helpers for handling chemical formulas in Python',
      author='Mostly adopted from work of Christoph Gohlke, by Benjie Chen',
      author_email='benjie@alum.mit.edu',
      packages=["pychemy"],
      package_dir={"pychemy": "."},
     )
