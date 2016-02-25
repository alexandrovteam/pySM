from setuptools import setup, find_packages

from pySM import __version__

setup(name='pySM',
      version=__version__,
      description='Scientific demonstration library for spatial metabolomics',
      url='https://github.com/alexandrovteam/pySM',
      author='Andrew Palmer, Alexandrov Team, EMBL',
      packages=find_packages())
