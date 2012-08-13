# Parallel Biclustering Algorithm - Fast Algorithm for finding all biclusters
# in a GEM
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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
#
# Contact Info:
#   Luke Imhoff (imho0030@umn.edu)
#   220 Delaware St. SE
#   Minneapolis, MN 55455
"""Array utilities missing from scipy/numpy, but needed for porting matlab code

@author Luke Imhoff
@license GPLv2
"""

# use full name for import so pylint doesn't  complain
import numpy.core.multiarray

def sortRows(matrix, columns=None):
    """Sorts rows of matrix based on columns
    
    Should work similar to sortrows in Matlab.
    @param matrix matrix whose rows to sort
    @type matrix 2D array
    @return (sortedMatrix, indices)
            sortedMatrix - matrix of sorted rows
            indices - index matrix with number indicating matrix row in
                      sortedMatrix
    """
    
    # use all columns if None given
    if columns == None:
        columns = numpy.core.multiarray.arange(matrix.shape[1])
    
    # reverse indices
    columns = columns[::-1]
    
    lexsortable = tuple(matrix[:, columns].transpose())
    indices = numpy.core.multiarray.lexsort(lexsortable)
    
    return (matrix[indices], indices)
