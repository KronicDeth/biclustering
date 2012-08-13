# Parallel Biclustering Algorithm - Fast Algorithm for finding all biclusters in a GEM
# Copyright (C) 2006  Luke Imhoff
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Contact Info:
#   Luke Imhoff (imho0030@umn.edu)
#   220 Delaware St. SE
#   Minneapolis, MN 55455
"""Helper routines for sizing structures

@author Luke Imhoff
@license GPLv2
"""

import numpy
import tables

def sizeAtom(dimSize):
    """Returns Atom class with minimum number of bits to store this dimension
    
    @param dimSize max number of elements in dimension
    """
    
    atomClass = None
    for i in (8, 16, 32, 64):
        if dimSize <= 2 ** i:
            atomClass = eval("tables.UInt" + str(i) + "Atom")
            break
    
    if atomClass is None:
        raise ValueError("Dimension too large to represent in atom")
    
    return atomClass

def sizeCol(dimSize):
    """Returns Col class with minimum number of bits to store this dimension
    
    @param dimSize max number of elements in dimension
    """
    
    colClass = None
    for i in (8, 16, 32, 64):
        if dimSize <= 2 ** i:
            colClass = eval("tables.UInt" + str(i) + "Col")
            break
    
    if colClass is None:
        raise ValueError("Dimension too large to represent in atom")
    
    return colClass

def sizeArray(dimSize):
    """Returns numpy array class with minimum number of bits to store this
    dimension
    
    @param dimSize max number of elements in dimension
    """
    
    arrayClass = None
    for i in (8, 16, 32, 64):
        if dimSize <= 2 ** i:
            arrayClass = eval("numpy.uint" + str(i))
            break
    
    if arrayClass is None:
        raise ValueError("Dimension too large to represent in array")
    
    return arrayClass