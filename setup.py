#
# SetupTools script for bnglonglat
#
from setuptools import setup, find_packages


# Load text for description and license
with open('README.md') as f:
    readme = f.read()


# Go!
setup(
    name='bnglonglat',

    version='0.0.1',

    description='Pure python port of convertbng.convert_longlat.',
    long_description=readme,
    long_description_content_type='text/markdown',

    license='MIT license',

    author='Michael Clerx',
    author_email='mail@michaelclerx.com',

    packages=find_packages(include=('bnglonglat', 'bnglonglat.*')),

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

    # Classifiers for pypi
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
)
