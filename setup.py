#!/usr/bin/env python


from setuptools import setup, find_packages

__version__ = "1.7.0"

# PyQt 5 has not been added here are. It needs to be installed via apt-get which is a Tigger v1.6.0 requirement.
requirements = ['astro_kittens==1.4.3', 'numpy==1.18.1', 'scipy==1.6.2', 'astlib==0.11.6', 'astropy==4.2', 'future', 'python-casacore==3.4.0']

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
    python_requires='>=2.7.0',
    description="Python libraries and command-line tools for manipulating Tigger-format LSMs",
    author="Oleg Smirnov",
    author_email="osmirnov@gmail.com",
    url="https://github.com/ska-sa/tigger-lsm",
    install_requires=requirements,
)

