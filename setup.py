#!/usr/bin/env python


from setuptools import setup, find_packages

__version__ = "1.7.0"

# PyQt has not been added here are it needs to be installed via apt-get which is a Tigger requirement.
# Versions below are set to astLib 0.11.6 tested and compatible versions found
requirements = ['astro_kittens==1.4.3',
                'numpy==1.18.1',  # set to astLib recommended
                'scipy==1.5.2',  # recommends 1.3.1, this fails, next available version
                'astlib==0.11.6',  # latest version that uses astropy WCS at the backend
                'astropy==3.2.3',  # recommends 3.2.1, this fails, next available version (last of 3.x)
                'future']

scripts = [
    'Tigger/bin/tigger-convert',
    'Tigger/bin/tigger-make-brick',
    'Tigger/bin/tigger-restore',
    'Tigger/bin/tigger-tag',
]

package_data = {
}

extras_require = {
}


setup(
    name="astro-tigger-lsm",
    version=__version__,
    packages=find_packages(),
    extras_require=extras_require,
    scripts=scripts,
    package_data=package_data,
    description="Python libraries and command-line tools for manipulating Tigger-format LSMs",
    author="Oleg Smirnov",
    author_email="osmirnov@gmail.com",
    url="https://github.com/ska-sa/tigger-lsm",
    install_requires=requirements,
)

