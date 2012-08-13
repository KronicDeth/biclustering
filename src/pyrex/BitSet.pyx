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
"""Collections implemented as using bit vectors

Bit vectors can only represent binary values for collection elements, but have
lower memory requirements and faster speed for operations that can be 
represented as low-level bit ops

@author Luke Imhoff
@license GPLv2
"""

cimport c_numpy
cimport c_python
import numpy

c_numpy.import_array()

import scipy
import sys

cdef enum:
    cBITS = 32
    cSHIFT = 5
    cMASK = 0x1F

# python accessible 
BITS = cBITS

# 32 bit population count from AMD Athlon optimization guide
cdef unsigned int populationCount(unsigned long v):
    cdef unsigned long w
    cdef unsigned long x
    cdef unsigned long y
    cdef unsigned long z
    
    w = v - ((v >> 1) & 0x55555555)
    x = (w & 0x33333333) + ((w >> 2) & 0x33333333)
    y = (x + (x >> 4)) & 0x0f0f0f0f
    z = (y * 0x01010101) >> 24
    
    return z

cdef unsigned int mask(unsigned char bits):
    """Returns mask for masking lower bits
    
    @param bits number of bits to mask
    @return mask
    """
    cdef unsigned int bitMask
    cdef unsigned char i
    cdef unsigned int bitsMask
    
    bitsMask = 0
    if bits == 0:
        return bitsMask
    
    bitMask = 0x1
    for i from 1 <= i <= bits:
        bitsMask = bitsMask | bitMask
        bitMask = bitMask << 1
    
    return bitsMask

cdef unsigned long vectorSize(unsigned long universe):
    size = (universe + cBITS - 1) / cBITS
    # min vector size is 2 for storage reasons as singleton arrays have the
    # habit of being converted to scalars
    if size == 1:
        size = 2
    
    return size

def arraySize(universe):
    """Returns size of array needed to hold a BitSet with the given universe
    
    @param universe size of universe of set
    """
    return int(vectorSize(universe))

cdef class BitSet:
    """Set stored in a bit vector"""
    
    cdef unsigned long _universe
    cdef c_numpy.ndarray _vector
    cdef unsigned int _cachedPopCount
    cdef char _cached
    
    def __init__(self, universe, initial=None, formatted=False):
        """
        BitSet(universe [, initial])
        
        @param universe size of universe.  max size of set
        @param initial initial data set
        @param formatted True if initial is already formatted to internal format
                         Only used for quick reconstructions from persistent
                         storage.  Array returned from asArray() is formatted.
        """
        
        self._universe = universe
        size = vectorSize(self._universe)
        
        if formatted:
            if (initial.size != size or
                initial.dtype != numpy.dtype(numpy.uint32)):
                raise ValueError("initial is not properly formatted")
            self._vector = initial
        else:
            self._vector = scipy.zeros(size, dtype = numpy.uint32)
            if initial is not None:
                self.initialize(initial)
        
        self._cached = 0
    
    cdef initialize(self, initial):
        cdef unsigned int index
        cdef unsigned int bit
        cdef unsigned int cElement
        cdef unsigned int *data
        
        data = <unsigned int *> self._vector.data
        
        for element in initial:
            # convert to c variable
            cElement = element;
            
            if cElement >= self._universe:
                raise IndexError("element not in universe")    
            
            index = cElement >> cSHIFT
            bit = cElement & cMASK
            data[index] = data[index] | (0x1 << bit)
    
    def __hash__(self):
        return sum(self._vector) % sys.maxint
    
    def __eq__(self, object obj):
        if type(obj) is not BitSet:
            return False
        
        cdef BitSet bitSet
        # cast to type
        bitSet = obj
        
        if self._universe != bitSet._universe:
            return False
        
        cdef unsigned int *selfData
        cdef unsigned int *bitSetData
        
        selfData = <unsigned int *> self._vector.data
        bitSetData = <unsigned int *> bitSet._vector.data
        
        cdef unsigned int i
        
        for i from 0 <= i < self._vector.size:
            if selfData[i] != bitSetData[i]:
                return False
        
        return True
    
    def __ne__(self, object obj):
        return not self.__eq__(obj)
    
    def __contains__(self, element):
        """Returns whether the element is in the set
        
        @param element element to check for in set
        """
        
        if element >= self._universe:
            return False
        
        cdef unsigned int cElement
        cdef unsigned int index
        cdef unsigned int bit
        cdef unsigned int *data
        
        cElement = element
        index = cElement >> cSHIFT
        bit = cElement & cMASK
        
        data = <unsigned int *> self._vector.data
        
        if data[index] & (0x1 << bit) != 0:
            return True
        else:
            return False
    
    def universe(self):
        """Returns size of universe"""
        return self._universe
    
    def __len__(self):
        """Returns number of elements in set"""
        
        if self._cached == 1:
            return self._cachedPopCount
        
        cdef unsigned long count
        cdef unsigned int *data
        cdef unsigned int i
        
        data = <unsigned int *> self._vector.data
        count = 0
        for i from 0 <= i < self._vector.size:
            count = count + populationCount(data[i])
        
        self._cachedPopCount = count
        self._cached = 1
        
        return count
    
    def intersection(BitSet self, object obj):
        """Returns intersection of this and bitSet
        
        @param bitSet bit set to intersection with self
        @return intersection of self and bitSet
        """
        
        if not isinstance(obj, BitSet):
            raise TypeError("Can only produce intersection with another BitSet")
        
        cdef BitSet bitSet
        bitSet = obj
        
        if self._universe != bitSet._universe:
            raise ValueError("BitSet Universe sizes do not match")
        
        cdef BitSet intersection
        intersection = BitSet(self._universe)
        
        # set intersection w/o constructor as data already in BitSet format
        intersection._vector = self._vector & bitSet._vector
        
        return intersection
    
    def __and__(BitSet self, object obj):
        """Returns intersection of this and bitSet
        
        @param bitSet bit set to intersection with self
        @return intersection of self and bitSet
        """
        return self.intersection(obj)
    
    def isSingletonIntersection(BitSet self, object obj, singleton):
        """Return whether intersection of this and bitSet has only 1 element
        
        Equivalent to
        def isSingletonIntersection(self, obj):
            return len(self & obj) == 1
        But faster as no intermediate BitSet is constructed
        @param
        @param
        @param singleton only element in intersection
        """
        cdef BitSet other
        
        if not isinstance(obj, BitSet):
            raise TypeError("Can only produce intersection with another BitSet")
        
        other = obj
        
        cdef unsigned int cElement
        cdef unsigned int index
        cdef unsigned int bit
        
        cSingleton = singleton
        index = cSingleton >> cSHIFT
        bit = cSingleton & cMASK
        
        cdef unsigned int *selfData
        cdef unsigned int *otherData
        cdef unsigned int i

        selfData = <unsigned int *> self._vector.data
        otherData = <unsigned int *> other._vector.data
        
        if selfData[index] & otherData[index] != 1 << bit:
            return False
        
        for i from 0 <= i < index:
            if selfData[i] & otherData[i] != 0:
                return False
        
        for i from index < i < self._vector.size:
            if selfData[i] & otherData[i] != 0:
                return False
        
        return True
    
    def union(BitSet self, object obj):
        """Returns union of this and obj
        
        @param obj bit set to union with self
        @return union of self and obj
        """
        
        if not isinstance(obj, BitSet):
            raise TypeError("Can only produce union with another BitSet")
        
        cdef BitSet bitSet
        bitSet = obj
        
        if self._universe != bitSet._universe:
            raise ValueError("BitSet Universer size do not match")
        
        cdef BitSet union
        union = BitSet(self._universe)
        
        # set union w/o constructor as data already in BitSet format
        union._vector = self._vector | bitSet._vector
        
        return union
    
    def __or__(BitSet self, object obj):
        """Returns union of this and obj
        
        @param obj bit set to union with self
        @return union of self and obj
        """
        return self.union(obj)
    
    def issubset(BitSet self, object obj):
        """Test whether every element in self is in obj"""
        if (not isinstance(obj, BitSet)):
            return False
        
        cdef BitSet bitSet
        bitSet = obj
        
        if (self._universe != bitSet._universe):
            return False
        
        cdef unsigned int *selfData
        cdef unsigned int *bitSetData
        
        selfData = <unsigned int *> self._vector.data
        bitSetData = <unsigned int *> bitSet._vector.data
        
        cdef unsigned int i
        
        for i from 0 <= i < self._vector.size:
            if selfData[i] & bitSetData[i] != selfData[i]:
                return False
        
        return True
    
    def __le__(BitSet self, object obj):
        """Test whether every element in self is in obj"""
        return self.issubset(obj)
    
    def issuperset(BitSet self, object obj):
        """Test whether every element in t is in s"""
        if (not isinstance(obj, BitSet)):
            return False
        
        cdef BitSet bitSet
        bitSet = obj
        
        if (self._universe != bitSet._universe):
            return False
        
        cdef unsigned int *selfData
        cdef unsigned int *bitSetData
        
        selfData = <unsigned int *> self._vector.data
        bitSetData = <unsigned int *> bitSet._vector.data
        
        cdef unsigned int i
        
        for i from 0 <= i < self._vector.size:
            if selfData[i] & bitSetData[i] != bitSetData[i]:
                return False
        
        return True
    
    def __ge__(BitSet self, object obj):
        return self.issuperset(obj)
    
    def complement(BitSet self):
        """Return U - self or the complement of the set"""
        cdef c_numpy.ndarray complementVector
        complementVector = ~self._vector
        
        # clean word that is not completely filled by universe
        cdef unsigned int *data
        data = <unsigned int *> complementVector.data
        
        cdef unsigned long index
        index = self._universe >> cSHIFT
        data[index] = data[index] & mask(self._universe & cMASK)
        
        # make sure pad word is rezero'd
        if self._universe < cBITS:
            data[index + 1] = 0
        
        return BitSet(self._universe, complementVector, True)
        
    def __invert__(BitSet self):
        """Return U - self or the complement of the set"""
        return self.complement()
    
    def __str__(self):
        cdef unsigned long base
        cdef unsigned int wordIndex
        cdef unsigned int word
        cdef unsigned char bit
        cdef unsigned int *data
        
        strList = list()
        strList.append('{')
        
        first = True
        base = 0
        data = <unsigned int *> self._vector.data
        for wordIndex from 0 <= wordIndex < self._vector.size:
            word = data[wordIndex]
            bit = 0
            while word != 0:
                if word & 0x1 == 0x1:
                    if first:
                        first = False
                    else:
                        strList.append(", ")
                    strList.append(str(base + bit))
                bit = bit + 1
                word = word >> 1
            base = base + cBITS
        
        strList.append("}")
        return ''.join(strList)
    
    def __repr__(self):
        return "BitSet(%d, %s)" % (self._universe, list(self))
    
    def asArray(self):
        return self._vector.copy()
    
    def __iter__(self):
        return BitSetIterator(self)

cdef class BitSetIterator:
    
    # kept reference so peek at data vector is legal
    cdef BitSet _bitSet
    
    # iterator state
    cdef unsigned int _size
    cdef unsigned int *_data
    cdef unsigned long _base
    cdef unsigned int _wordIndex
    cdef unsigned int _word
    cdef unsigned char _bit
    
    def __init__(self, BitSet bitSet):
        self._bitSet = bitSet
        self._size = bitSet._vector.size
        self._data = <unsigned int *> bitSet._vector.data
        self._base = 0
        self._wordIndex = 0
        self._word = self._data[0]
        self._bit = 0
    
    def __iter__(self):
        return self
        
    def __next__(self):
        if self._wordIndex >= self._size:
            raise StopIteration
        
        retVal = -1
        for wordIndex from self._wordIndex <= wordIndex < self._size:
            while self._word != 0:
                if self._word & 0x1 == 0x1:
                    retVal = self._base + self._bit
                self._bit = self._bit + 1
                self._word = self._word >> 1
                if retVal >= 0:
                    return int(retVal)
            self._base = self._base + cBITS
            self._word = self._data[wordIndex + 1]
            self._wordIndex = self._wordIndex + 1
            self._bit = 0

        raise StopIteration
