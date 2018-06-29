roboptim-core-plugin-ipopt
==========================

[![License LGPL 3][badge-license]](http://www.gnu.org/licenses/lgpl-3.0.txt)
[![Build Status](https://travis-ci.org/roboptim/roboptim-core-plugin-ipopt.png?branch=master)](https://travis-ci.org/roboptim/roboptim-core-plugin-ipopt)
[![codecov.io](https://codecov.io/github/roboptim/roboptim-core-plugin-ipopt/coverage.svg?branch=master)](https://codecov.io/github/roboptim/roboptim-core-plugin-ipopt?branch=master)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/877/badge.svg)](https://scan.coverity.com/projects/877)
[![Join the chat at https://gitter.im/roboptim/development](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/roboptim/development?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


[![Zenodo](https://zenodo.org/badge/doi/10.5281/zenodo.10332.png)](http://zenodo.org/record/10332)

This package is the IPOPT plug-in of the RobOptim framework. It is
released under the [LGPL-3](COPYING.LESSER) license.

**Warning:** this repository contains [Git
submodules][git-submodules]. Please clone this repository using the
`git clone --recursive` command. If you already have cloned the
repository, you can run `git submodule init && git submodule update`
to retrieve the submodules.


For general information about the project, please refer to its
homepage: http://www.roboptim.net/


Documentation
-------------

To get started with this library, please read the [online Doxygen
documentation][doxygen-documentation].

It can also be generated locally by running the `make doc`
command. After the package is installed, the documentation will be
located in the `$prefix/share/doc/roboptim-core` directory where
`$prefix` is your installation prefix (`/usr/local` by default).


Getting Help
------------

Support is provided through:
 * the RobOptim mailing-list: roboptim@googlegroups.com
 * the following Gitter room: https://gitter.im/roboptim/development


Contributing
------------

If you want to contribute, please refer to the
[CONTRIBUTING.md](CONTRIBUTING.md) file.


Credits
-------

This package authors are credited in the [AUTHORS](AUTHORS) file.


Available packages
------------------

 * Arch Linux (Git master branch):
   https://aur.archlinux.org/packages/roboptim-core-plugin-ipopt-git/
 * Mac OS X Homebrew Formula (Git HEAD):
   https://www.github.com/roboptim/homebrew-roboptim

[badge-license]: https://img.shields.io/badge/license-LGPL_3-green.svg

[doxygen-documentation]: http://www.roboptim.net/roboptim-core-plugin-ipopt/doxygen/HEAD/

[git-submodules]: http://git-scm.com/book/en/Git-Tools-Submodules

[Boost]: http://www.boost.org/
[CMake]: htttp://www.cmake.org/
[Doxygen]: http://www.stack.nl/~dimitri/doxygen/
[Eigen]: http://eigen.tuxfamily.org/
[Git]: http://git-scm.com/
[Libtool]: https://www.gnu.org/software/libtool/
[pkg-config]: http://www.freedesktop.org/wiki/Software/pkg-config/
[RobotPkg]: http://robotpkg.openrobots.org/
