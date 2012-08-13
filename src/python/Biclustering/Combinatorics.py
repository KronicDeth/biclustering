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
"""Combinatoric utilities not included in scipy or faster implementation

@author Luke Imhoff
@license GPLv2
"""

from scipy import zeros
from scipy.misc import factorial

import operator

def nChooseK(n, k):
    """n choose k
    
    Returns number of combinations of k elements in n elements
    @param n total number of elements in set
    @param k number of elements in subset to choose
    @return C(n, k)
    """
    numerator = 1
    for i in xrange(max(n - k, k) + 1, n + 1):
        numerator *= i
    
    return numerator / factorial(min(n - k, k), exact = 1)

class xcombinations(object):
    """Returns all combinations of subsetSize number from [0, setSize)
    
    Supports len()
    @param setSize number of elements in set to select from
    @param subsetSize number of elements chosen from set
    @return generator of combinations.  Combination is stored in a 1D array of subsetSize length.
            None if subsetSize is 0.
    """
    
    def __init__(self, setSize, subsetSize):
        if setSize < subsetSize:
            raise ValueError("setSize (%d) must be >= subsetSize (%d)" % 
                             (setSize, subsetSize))
        
        self.setSize = setSize
        self.subsetSize = subsetSize
    
    def len(self):
        return nChooseK(self.setSize, self.subsetSize)
    
    def __iter__(self):
        a = zeros(self.subsetSize)
        # using a closure faster than passing setSize
        def recursiveXCombinations(setStart, recursiveSubsetSize):
            if recursiveSubsetSize == 0:
                yield None
            else:
                index = self.subsetSize - recursiveSubsetSize
                for i in xrange(setStart, self.setSize - recursiveSubsetSize  + 1):
                    a[index] = i
                    for j in recursiveXCombinations(i + 1, recursiveSubsetSize - 1):
                        yield a
        
        return recursiveXCombinations(0, self.subsetSize)

def permutations(setSize, subsetSize):
    """Return 'setSize permute subsetSize' or nPk
    
    @param setSize number of elements in set to select from
    @param subsetSize number for elements chosen from set
    @return number of permutations of subsetSize of setSize
    """
    count = 1
    for i in xrange(setSize - subsetSize + 1, setSize + 1):
        count *= i
    
    return count

class xpermutations(object):
    """Returns all permutations of subsetSize number from [0, setSize)
    
    Supports len()
    @param setSize number of elements in set to select from
    @param subsetSize number for elements chosen from set
    @return generator of permutations.  Permutations stored in a 1D array of
            subsetSize length. None if subsetsize is 0
    """
    
    def __init__(self, setSize, subsetSize):
        if setSize < subsetSize:
            raise ValueError("setSize (%d) must be >= subsetSize (%d)" % 
                             (setSize, subsetSize))
        
        self.setSize = setSize
        self.subsetSize = subsetSize
    
    def len(self):
        return permutations(self.setSize, self.subsetSize)
    
    def __iter__(self):
        a = zeros(self.subsetSize)
        items = range(self.setSize)
        # using a closure faster than passing setSize
        def recursiveXPermutations(recursiveItems, recursiveSubsetSize):
            if recursiveSubsetSize == 0:
                yield None
            else:
                index = self.subsetSize -recursiveSubsetSize
                for i in xrange(len(recursiveItems)):
                    a[index] = recursiveItems[i]
                    for cc in recursiveXPermutations(recursiveItems[:i] + recursiveItems[i + 1:], recursiveSubsetSize - 1):
                        yield a
        
        return recursiveXPermutations(items, self.subsetSize)

def selections(setSizes):
    length = 1
    for setSize in setSizes:
        length *= setSize
    
    return length

class xselections(object):
    """Returns all selections of one item from each set of setSizes
    
    Supports len()
    @param collection of set sizes
    @return generator of selections.  Selection is stored in a 1D array with
            each entry corresponding to the selected index from the respective
            set.
            None if no sets
    """
    
    def __init__(self, setSizes):
        self.setSizes = setSizes
    
    def __len__(self):
        return selections(self.setSizes)
    
    def len(self):
        return selections(self.setSizes)
    
    def __iter__(self):
        a = zeros(len(self.setSizes))
        setCount = len(self.setSizes)
        
        def recursiveXSelections(setIndex):
            if setIndex == setCount:
                yield None
            else:
                for i in xrange(0, self.setSizes[setIndex]):
                    a[setIndex] = i
                    for j in recursiveXSelections(setIndex + 1):
                        yield a
        
        return recursiveXSelections(0)
