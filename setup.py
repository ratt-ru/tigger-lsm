#!/usr/bin/env python3


from setuptools import setup, find_packages

__version__ = "1.7.2"

# PyQt 5 has not been added here are. It needs to be installed via apt-get which is a Tigger v1.6.0 requirement.
requirements = ['astro_kittens', 'numpy', 'scipy', 'astlib', 'astropy', 'future', 'python-casacore']

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
    python_requires='>=3.6',
    description="Python libraries and command-line tools for manipulating Tigger-format LSMs",
    author="Oleg Smirnov",
    author_email="osmirnov@gmail.com",
    url="https://github.com/ska-sa/tigger-lsm",
    install_requires=requirements,
)

