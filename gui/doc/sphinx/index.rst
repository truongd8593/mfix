MFIX-UI developer documentation
===============================

This documentation describes the MFIX graphical user interface API.

(FIXME: The following introductory and license text was taken from the SRS and
doc strings in the MFIX code. Please correct or ammend gui/doc/sphinx/index.rst)

MFIX [Multiphase Flow with Interphase eXchanges] is a general-purpose computer
code developed at the National Energy Technology Laboratory [NETL] for
describing the hydrodynamics, heat transfer and chemical reactions in
fluid-solid systems.

The purpose of this application is to provide a user-friendly interface to:

* Set up the MFIX simulation
* Control and interact with the MFIX code during execution
* Provide post-processing capabilities for analyzing output data

Please visit: https://mfix.netl.doe.gov/

License
-------

As a work of the United States Government, this project is in the public domain
within the United States. As such, this code is licensed under
CC0 1.0 Universal public domain.

Developing
----------

This application is written primarily in Python. The Anaconda Python
distribution is the officially supported environment. All core Python
libraries used are provided by Anaconda. As of version {{ version }} of this
application, development is done with Anaconda 4 Python 2.7 and 3.5
distributions.

Additional packages and libraries are required for development, testing and
building the application.

VTK, numpy, F2PY and Flask are build requirements. To execute the included
tests, nose and Xvfb are required. The documentation makefile target requires
Sphinx.

See {{ link['Getting Started'] }} for development environment set up examples.


Contents
--------

.. toctree::
   :maxdepth: 4

   modules



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

