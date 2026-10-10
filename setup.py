from setuptools import setup, find_packages
import os

def read_version():
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, "brightify", "__init__.py"), "r") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip("'\"")
    return "0.0.0" # (in case nothing is found)

setup(
    name='brightify',
    version=read_version(), # This reads __init__.py and its __version__ variable. If you want to update the version number, please do so in the __init__.py file.
    author='Mina Akhyani',
    description='A python package for brightness calculation based on MCPL files',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/BrightnessTools/Brightify',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.22',
        'pandas>=1.4',
        'matplotlib>=3.4',
        'tqdm>=4.66',
        'mcpl',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # adjust as needed
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)