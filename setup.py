#!/usr/bin/env python


from setuptools import setup, find_packages

__version__ = "1.7.0"

requirements = ['astro_kittens==1.4.3',
                'numpy>=1.17',
                'scipy==1.5.2',
                'astlib==0.10.2',
                'astropy==4.1',
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

