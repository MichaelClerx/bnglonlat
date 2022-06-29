#
# SetupTools script for bnglonlat
#
from setuptools import setup, find_packages


# Load text for description and license
with open('README.md') as f:
    readme = f.read()


# Go!
setup(
    # Module name (lowercase) and version (string)
    name='bnglonlat',
    version='0.0.1',

    # Description
    description='Pure python port of convertbng.util.convert_lonlat.',
    long_description=readme,
    long_description_content_type='text/markdown',

    # Author and license
    license='MIT license',
    author='Michael Clerx',
    author_email='work@michaelclerx.com',

    # Project URLs (only first is required, rest is for PyPI)
    url='https://github.com/MichaelClerx/bnglonlat',
    project_urls={
        'Bug Tracker': 'https://github.com/MichaelClerx/bnglonlat/issues',
        'Source Code': 'https://github.com/MichaelClerx/bnglonlat',
    },

    # Classifiers for pypi
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],

    # Packages to include
    packages=find_packages(include=('bnglonlat', 'bnglonlat.*')),

    # List of dependencies
    install_requires=[
        'numpy',
    ],

    # Optional extras
    extras_require={
        'dev': [
            'flake8>=3',
            'convertbng>=0.6.32',
        ],
    },

)
