# -*- coding: UTF-8 -*-
#! /usr/bin/python
import sys
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

NAME    = 'splitri'
VERSION = '0.0.1'
AUTHOR  = 'Ahmed Ratnani'
EMAIL   = 'ratnaniahmed@gmail.com'
URL     = 'http://www.ratnani.org/splitri/'
DESCR   = 'Splines on Triangulations.'
KEYWORDS = ['CAD', 'FEM', 'IGA', 'OpenGL']
LICENSE = "LICENSE"

setup_args = dict(
    name             = NAME,
    version          = VERSION,
    description      = DESCR,
    long_description = open('README.md').read(),
    author           = AUTHOR,
    author_email     = EMAIL,
    license          = LICENSE,
    keywords         = KEYWORDS,
    url              = URL,
#    download_url     = URL+'/get/default.tar.gz',
)

packages=[  'splitri' \
          , 'splitri.core' \
          , 'splitri.gallery' \
          , 'splitri.utils' \
         ]
package_dir={  'splitri': 'splitri'\
              ,'caid.core':  'caid/core' \
              ,'caid.gallery':  'caid/gallery' \
              ,'caid.utils':  'caid/utils' \
              ,}

ext_modules  = [ \
                # ... bsplines extension
                 Extension('caid.core.bezier', \
                           sources = ['caid/core/bezier.pyf', \
                                      'caid/core/bezier.F90'], \
                           f2py_options = ['--quiet'], \
                           define_macros = [ \
                                            #('F2PY_REPORT_ATEXIT', 0),
                                            ('F2PY_REPORT_ON_ARRAY_COPY', 0)] \
                          ) \
                # ...
                ,]

def setup_package():
    if 'setuptools' in sys.modules:
        setup_args['install_requires'] = ['numpy']
    setup(  packages = packages \
          , package_dir=package_dir \
#          , ext_modules=ext_modules \
          , **setup_args)

if __name__ == "__main__":
    setup_package()
