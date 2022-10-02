from gettext import install
from setuptools import setup

setup(
    name='Coordinate',
    version='0.0.1',
    description='Utilities for coordination chemistry and molecular dft automation with TURBOMOLE',
    py_modules=["core"],
    package_dir={'':"core"},
    install_requires=["matplotlib==3.5.1", "networkx==2.8.6", "numpy==1.21.5", "pexpect==4.8.0", "PyYAML==6.0", "scipy==1.8.0", "setuptools==59.6.0", "turbomoleio"],
    author='William Laderer',
    author_email='wtladerer@gmail.com',
    url='https://github.com/wladerer/Coordinate'
)

