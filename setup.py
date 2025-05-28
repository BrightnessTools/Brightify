from setuptools import setup, find_packages

setup(
    name='Brightify',
    version='0.1.0',
    author='Mina Akhyani',
    description='A python package for brightness calculation based on MCPL files',
    long_description=open('README.md').read(),
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
