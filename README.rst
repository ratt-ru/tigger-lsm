==========================
Tigger-LSM: LSM Libs/utils
==========================

N.B. THIS README IS THE BETA TESTER VERSION
===========================================

Installing Tigger-LSM
=====================

Python dependencies
-------------------

Automatically installed dependencies:

* astro_kittens == v1.4.3
* numpy >= v1.17
* scipy == v1.5.2
* astlib == v0.10.2
* astropy == v4.1
* future

from source with Ubuntu 20.04
-----------------------------

Build the source with the following::

    git clone https://github.com/razman786/tigger_lsm_pyqt5
    cd tigger_lsm_pyqt5
    python setup.py install --user

Using Tigger-LSM
================

In python:

    import Tigger
    model = Tigger.load("foo.lsm.html")

In the shell

    tigger-convert foo.txt foo.lsm.html


Beta Tester Questions or problems
=================================

Open an issue on github

https://github.com/razman786/tigger_lsm_pyqt5/issues


Travis
======

.. image:: https://travis-ci.org/ska-sa/tigger-lsm.svg?branch=master
    :target: https://travis-ci.org/ska-sa/tigger-lsm
