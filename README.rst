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

from source::

    $ git clone https://github.com/razman786/tigger_lsm_pyqt5
    $ cd tigger_lsm_pyqt5
    $ git checkout develop
    $ python setup.py install --user


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
