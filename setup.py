#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

try:
    import tables
except ImportError:
    print "ERROR: Need to install PyTables"
    sys.exit(1)

try:
    from Pyrex.Distutils import build_ext
except ImportError:
    print "ERROR: Need to install pyrex to compile extensions"
    sys.exit(1)

try:
    import numpy
except ImportError:
    print "ERROR: Need to install Numpy"
    sys.exit(1)

import sys

setup(name='BICPBS',
      version='0.2.1',
      description='Binary Indexed Chaining Parallel Biclustering System',
      author='Luke Imhoff',
      author_email='imho0030@umn.edu',
      url='http://tewfik.dtc.umn.edu/',
      package_dir = {'': 'src/python'},
      packages=['Biclustering'],
      ext_modules=[ 
          Extension("Biclustering.BitSet",
                    ["src/pyrex/BitSet.pyx"], #"src/pyrex/c_numpy.pxd", "src/pyrex/c_python.pxd"],
                    include_dirs = [sys.prefix + '/lib/python' + 
                                    sys.version[:3] + 
                                    '/site-packages/numpy/core/include/'])
                  ],
      cmdclass = {'build_ext': build_ext}
     )