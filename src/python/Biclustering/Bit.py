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
"""Bit-based data structures

@author Luke Imhoff
@license GPLv2
"""

import itertools

import numpy.core.ma
import numpy.core.multiarray

import tables

import Biclustering.BitSet

class OrderedBitSet(object):
    """An ordered set (a collection of numbers where order matters but no
    repeated elements) that uses a BitSet for faster membership calculations"""
    
    def __init__(self, order, universe=None, bitSet=None):
        """OrderedBitSet(order, universe)
            OR
        OrderedBitSet(order, set)
        """
        self.order = order
        
        if bitSet is None:
            bitSet = Biclustering.BitSet.BitSet(universe, order)
        self.set = bitSet
    
    def __contains__(self, element):
        return element in self.set
    
    def __len__(self):
        return self.order.size
    
    def __str__(self):
        return self.order.__str__()
    
    def __getitem__(self, index):
        return self.order[index]
    
    def setIntersection(self, other):
        """Returns the intersection of this set and the other
        
        @param other other set to intersect
        @return set (non-ordered) that is the intersecton of this set and the other set
        """
        return self.set & other.set
    
    def isSingletonIntersection(self, other, singleton):
        """Returns whether this set and the other set only have the singleton
        member in common
        
        @param other other OrderedBitSet
        @param singleton number of member
        """
        return self.set.isSingletonIntersection(other.set, singleton)
    
    def reverse(self):
        """Return new OrderedBitSet with order reverse of this one"""
        return OrderedBitSet(self.order[::-1], set = self.set)
    
    def chain(self, tail):
        """Creates new OrderedBitSet with the union of self and tail sets and
        tail[-1] appended to self's order
        
        @param tail OrderedBitSet to append
        """
        order = numpy.core.multiarray.concatenate((self.order, tail.order[1:]))
        bitSet = self.set | tail.set
        
        return OrderedBitSet(order, set = bitSet)
    
    def isOrderedSubset(self, orderedBitSet):
        """Returns whether orderedBitSet is an ordered subset of self
        
        @param orderedBitSet bit set to test for ordered subset
        """
        
        # subset test
        if not self.set.issubset(orderedBitSet.set):
            return False
        
        # order test
        orderedSubset = True
        # ordered subset test for nested conditions
        lastMatch = -1
        superSize = orderedBitSet.order.size
        for condition in self.order:
            for i in xrange(lastMatch + 1, superSize):
                if orderedBitSet.order[i] == condition:
                    lastMatch = i
                    break
            else:
                orderedSubset = False
                break
        
        return orderedSubset

class OrderedBitSetAccessor(object):
    """Accessor to access rows of a Table as OrderedBitSets"""
    
    def __init__(self, name, width, universe):
        self.name = name
        self.width = width
        self.universe = universe
        self.set = BitSetAccessor(name + '/set', universe)
        
        orderCol = Biclustering.Sizing.sizeCol(universe)
    
        class OrderedBitSetTable(tables.IsDescription):
            """Table of OrderedBitSets of a single width"""
            
            order = orderCol(shape = width)
            set = self.set.description
        
        OrderedBitSetTable._v_flavor = 'numpy'
        
        self.description = OrderedBitSetTable()
    
    def unpack(self, row):
        """Retrieves data from row and presents it as an OrderedBitSet
        
        @param row row holding packed OrderedBitSet
        """
        order = row[self.name + '/order']
        bitSet = self.set.unpack(row)
        return OrderedBitSet(order, set = bitSet)
    
    def pack(self, row, orderedBitSet):
        """Stores OrderedBitSet in row
        
        @param row row in Table to store OrderedBitSet
        @param orderedBitSet Ordered BitSet to store in the row
        """
        row[self.name + '/order'] = orderedBitSet.order
        self.set.pack(row, orderedBitSet.set)
    
class OrderedSetArray(object):
    """Array of OrderedBitSets of a single width"""
    
    def __init__(self, nodeFile, where, name, width=None, universe=None):
        """Creates or loads OrderedSetArray nodeFile
        
        @param nodeFile file OrderedSetArray is in
        @param where path to parent of OrderedSetArray in nodeFile
        @param name name of node OrderedSetArray is
        @param width width of OrderedBitSets (only needed for creation)
        @param universe universe for each set (only needed for creation)
        """
        try:
            group = nodeFile.getNode(where, name)
        except tables.NoSuchNodeError:
            group = nodeFile.createGroup(where, name)
        
        try:
            self.orders = group.orders
        except tables.NoSuchNodeError:
            ordersClass = Biclustering.Sizing.sizeAtom(universe)
            shape = (0, width)
            atom = ordersClass(shape = shape, flavor = 'numpy')
            self.orders = nodeFile.createEArray(group, "orders", atom)
        
        self.sets = Biclustering.Bit.SetArray(nodeFile, group, "sets",
                                              universe)
    
    def append(self, orderedBitSet):
        """Appends data in OrderBitSet to end of Array
        
        @param orderedBitSet Ordered BitSet to append to array
        """
        order = orderedBitSet.order.copy()
        order.shape = (1, order.size)
        self.orders.append(order)
        self.sets.append(orderedBitSet.set)
    
    def __iter__(self):
        for order, bitSet in itertools.izip(self.orders, self.sets):
            yield OrderedBitSet(order, set = bitSet)
    
    def __getitem__(self, index):
        bitSet = self.sets[index]
        return OrderedBitSet(self.orders[index], set = bitSet)
    
    def where(self, position, value):
        """Returns rows where position has value
        
        @param position column
        @param value value for which to search
        """
        return numpy.core.multiarray.where(self.orders[:, position] == value)
    
    def whereNot(self, value):
        """Returns rows where value is not part of the sets for those rows
        
        @param value member not in set
        """
        return self.sets.whereNot(value)
    
class BitSetAccessor(object):
    """Accessor to access rows of data as BitSets"""
    
    def __init__(self, name, universe):
        self.name = name
        self.universe = universe
        
        bitSet = Biclustering.Sizing.sizeCol(2 ** Biclustering.BitSet.BITS)
        bitSetShape = Biclustering.BitSet.arraySize(universe)
        
        self.description = bitSet(shape = bitSetShape)
    
    def unpack(self, row):
        """Presents data in row as a BitSet
        
        @param row row in which data is stored
        """
        numpyFormatArray = numpy.core.ma.asarray(row[self.name])
        return Biclustering.BitSet.BitSet(self.universe, numpyFormatArray,
                                          True)
    
    def pack(self, row, bitSet):
        """Stored bitSet into row
        
        @param row row in which to store data in bitSet
        @param bitSet bitSet with data
        """
        row[self.name] = bitSet.asArray()

class SetArray(object):
    """Array of BitSets"""
    
    def __init__(self, nodeFile, group, name, universe=None):
        """
        BitSetArray(file, group, universe)
            OR
        BitSetArray(file, group)
        
        @param file file to create/load array in/from
        @param group parent group of array
        @param name name of array
        @param universe universe size of BitSets in array.  Must be specified
               when creating.  If not given, then assume array is to be loaded
        """
        
        self.file = nodeFile
        
        try:
            self.bitSets = self.file.getNode(group, name)
            self.universe = self.file.getNodeAttr(self.bitSets, "universe")
        except tables.NoSuchNodeError:
            self.universe = universe
            
            shape = (0, Biclustering.BitSet.arraySize(universe))
            bitRange = 2 ** Biclustering.BitSet.BITS
            atomClass = Biclustering.Sizing.sizeAtom(bitRange)
            atom = atomClass(shape = shape, flavor = 'numpy')
            
            self.bitSets = self.file.createEArray(group, name, atom)
            
            self.file.setNodeAttr(self.bitSets, "universe", self.universe)
    
    def append(self, bitSet):
        """Appends bitSet to array
        
        @param bitSet BitSet to append
        """
        
        # reshape to match rank of EArray
        bitSetArray = bitSet.asArray()
        bitSetArray.shape = (1, bitSetArray.size)
        
        self.bitSets.append(bitSetArray)
    
    def __iter__(self):
        for row in self.bitSets:
            yield Biclustering.BitSet.BitSet(self.universe, row, True)
    
    def __getitem__(self, index):
        return Biclustering.BitSet.BitSet(self.universe, self.bitSets[index],
                                          True)
    
    def whereNot(self, value):
        """Returns array of indexes where value is not a member of the set
        
        @param value value not in set
        """
        # TDDO replace with BitSet function
        index, mask = divmod(value, Biclustering.BitSet.BITS)
        return numpy.core.multiarray.where(self.bitSets[:, index] & 
                           numpy.uint32(1 << mask) == 0)
        