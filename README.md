# splitri
Splitri is a Python package for Splines functions on triangulations. It is based on Bernstein representation with respect to the domain points (also called knots to refeer to the tensorial B-splines case).

Splitri's philosophy is the following:

1. first construct a triangulation (you may actually give your own triangulation, generated for example by gmsh, or any other mesh tool)
2. construct the refined triangulation with domain points (knots). For the moment, only **uniform** points are used (with arbitrary degree)
3. construct the Splines functions space given for each spline its local representation on the refined triangulation (B-net)
4. finally, you can export all needed data for a Finite Elements code

For more details, please read [**SPLITRI**](http://ratnani.org/splitri_doc/)

Requierements
=============

**numpy**
---------

[**NumPy**](http://www.numpy.org/) is the fundamental package for scientific computing with Python

Installation can be done using

   `sudo apt-get install python-numpy`

**scipy**
---------

[**SciPy**](http://www.scipy.org/) is a Python-based ecosystem of open-source software for mathematics, science, and engineering.

Installation can be done using

   `sudo apt-get install python-scipy`

You can install both **numpy** and **scipy** using 

    sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

**matplotlib**
--------------

[**matplotlib**](http://www.matplotlib.org/)  is a python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. matplotlib can be used in python scripts, the python and ipython shell

**igakit**
----------

[**igakit**](http://bitbucket.org/dalcinl/igakit) is a package that implements many of the NURBS routines from Piegl's book using Fortran and Python.

**caid**
------------
[**caid**](https://github.com/ratnania/caid) the GUI interface is not needed, only the caid package is needed.

Installation
============

Installation can be done by runing the following command, giving **PATH_FOR_INSTALLATION**

    python setup.py install --prefix=PATH_FOR_INSTALLATION 
