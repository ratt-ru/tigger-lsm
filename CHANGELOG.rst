================
Tigger Changelog
================

1.7.0
=====

* Ported from PyQt 4 to PyQt 5
* Works with Tigger 1.6.0 and Python 3
* Various bug fixes

1.3.9
=====

Changes since 1.3.8:
 
* Improve freq0 parsing logic (#84)
* Install vext.pyqt4 if in virtualenv (#86)
* Use KERN-2 in Dockerfile
* Define f0 outside the "if" block #87

1.3.8
=====

changes since 1.3.7:
 * hack around problem with pkg_resources package bug


1.3.5
=====

 * Provide MS list to tigger-convert --app-to-int operation (#69)
 * Tigger incompatible with pyfits>=3.4 (#71)
 * It's given the correct name (install_requires) so that it will actually have an
   effect on the package manager.
 * PyQt4 is removed from install_requires, since it is not a PyPA-installable package.
   Instead, a check is added to fail setup if it is not already installed.
 * Added missing scipy and pyfits dependencies.



1.3.3
=====

 * renamed package to astro-tigger to resolve name conflict on pypi


