==========================
Tigger-LSM: LSM Libs/utils
==========================

Installing Tigger-LSM
=====================

Ubuntu package
--------------

Enable the KERN suite and install the tigger-lsm package.


from pypi or from source
------------------------

requirements:

 * Assorted python packages:
    + astro-kittens==v1.4.3,
    + astropy==v4.1,
    + numpy==v1.18.1,
    + scipy==v1.5.2,
    + astLib==v0.11.6,
    + python-casacore==v3.4.0,
    + future==v0.18.2.

 With the exception of astLib, these are already present in most Linux
 distros.  astLib may be downloaded here: http://astlib.sourceforge.net/

 * Purr/Kittens. Available from pip as astro-kittens. Else, install the purr package from a MeqTrees binary
 distribution (see http://www.astron.nl/meqwiki/Downloading). Alternatively, check it out from svn (see below),
 and make sure the parent of the Kittens directory is in your PYTHONPATH.

To obtain on ubuntu you can run::

  $ sudo apt-get install python-kittens python-pyfits python-astlib python-scipy python-numpy

now from pip::

    $ pip install astro-tigger-lsm

or from source::

    $ git clone https://github.com/ska-sa/tigger-lsm
    $ cd tigger-lsm
    $ python setup.py install


Using Tigger-LSM
================

In python:

    $ import Tigger
    $ model = Tigger.load("foo.lsm.html")

In the shell

    $ tigger-convert foo.txt foo.lsm.html


Questions or problems
=====================

Open an issue on github

https://github.com/ska-sa/tigger-lsm


Travis
======

.. image:: https://travis-ci.org/ska-sa/tigger-lsm.svg?branch=master
    :target: https://travis-ci.org/ska-sa/tigger-lsm
