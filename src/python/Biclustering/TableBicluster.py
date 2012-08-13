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
#   3480 Golfview Dr Apt 1208
#   Eagan, MN 55123
"""Bicluster class and Group (for use with Pytables)

@author Luke Imhoff
@license GPLv2
"""

import itertools
import numarray
import numpy
import sys
import tables

import Biclustering.Array
import Biclustering.Bit
import Biclustering.Sizing

NESTED_LIST = ['nonnested', 'nested', 'unknown']
NESTED = tables.Enum(NESTED_LIST)
NESTED_ATOM = tables.EnumAtom(NESTED, dtype = 'UInt8', shape = (0,),
                              flavor = 'numpy')

class Group(object):
    """Group node for holding bilcusters"""
    
    def __init__(self, nodeFile, group, maxConditions, maxGenes, create=True,
                 minGenes=2):
        """Creates Bicluster Table
        
        @param nodeFile file to create Table on
        @param group group to which to attach Group
        @param maxConditions max condition indexes in any one bicluster
        @param maxGenes max gene indexes in any one bicluster
        @param minGenes min genes index for a valid bicluster
        @param create create bilcuster Table in as child of group in file.
                      Default = True.  If False checks for existence of table
                      in file
        """
        self.file = nodeFile
        
        self.maxConditions = maxConditions
        self.maxGenes = maxGenes
        self.minGenes = minGenes
        
        if create:
            self.biclusters = self.file.createGroup(group, "biclusters")
            self.file.setNodeAttr(self.biclusters, "minGenes",
                                  numarray.array((self.minGenes,)))
        else:
            self.biclusters = self.file.getNode(group, "biclusters")
            self.minGenes = self.file.getNodeAttr(self.biclusters,
                                                  "minGenes")[0]
        
        self.cache = WidthGroupCache(self.file, self.biclusters,
                                     self.maxConditions, self.maxGenes)
        
        # Chain Performance Monitors
        self.noHeadWidth = 0
        self.noHeadLink = 0
        self.noTailLink = 0
        self.redundantCondition = 0
        self.insufficientGenes = 0
    
    def pool(self, conditions, genes):
        """Pool biclusters
        
        @param conditions condition indexes of bicluster
        @param genes dependent indexes of the bicluster
        @param returns true if biclusters valid
        """
        
        if len(genes) >= self.minGenes:
            self.cache[len(conditions)].pool(conditions, genes)
            return True
        
        return False
    
    def index(self, width):
        """Returns index of width biclusters
        
        @param width width of biclusters to index
        """
        if width not in self.cache:
            return
        
        self.cache[width].index()
    
    def chain(self, headWidth, link):
        """Chains biclusters
        
        @param headWidth number of conditions in first bicluster
        @param link condition linking chain
        @return number of valid biclusters chained
        """
        
        if headWidth not in self.cache:
            self.noHeadWidth += 1
            return 0
        headGroup = self.cache[headWidth]
        headIndexes = headGroup[link]
        
        if len(headIndexes) == 0:
            self.noHeadLink += 1
        
        tailGroup = self.cache[2]
        tailIndexes = tailGroup['tails', link]
        
        if len(tailIndexes) == 0:
            self.noTailLink += 1
        
        progressBar = \
            Biclustering.Timing.ProgressBar(headIndexes.size,
                                            "  Link %d" % link)
        
        count = 0
        for headIndex in headIndexes:
            progressBar.update()
            
            # BUG FIX pytables doesn't understand numpy integer types
            headIndex = int(headIndex)
            headConditions = headGroup.conditions[headIndex]
            headGenes = headGroup.genes[headIndex]
            for tailIndex in tailIndexes:
                # BUG FIX pytables doesn't understand numpy integer types
                tailIndex = int(tailIndex)
                
                tailConditions = tailGroup.conditions[tailIndex]
                
                # if reduntant conditions besides linking condition
                if tailConditions.order[-1] in headConditions:
                    self.redundantCondition += 1
                    continue
                
                tailGenes = tailGroup.genes[tailIndex]
                genes = headGenes & tailGenes
                
                # if not enough common genes for valid bicluster
                geneCount = len(genes)
                if (geneCount < self.minGenes):
                    self.insufficientGenes += 1
                    continue
                
                conditions = headConditions.chain(tailConditions)
                
                self.pool(conditions, genes)
                count += 1
                # under special conditions merged biclusters can be pruned
                if geneCount == len(headGenes) == len(tailGenes):
                    headGroup.nested[headIndex] = NESTED.nested
                    tailGroup.nested[tailIndex] = NESTED.nested
        
        progressBar.finish()
        self.file.flush()
        
        return count
    
    def chainStats(self, reset=False):
        """Prints performance stats for chain
        
        @param reset [False] True to reset stats to 0.
        """
        print "No Head Width:", self.noHeadWidth
        print "No Head Link:", self.noHeadLink
        print "No Tail Link:", self.noTailLink
        print "Redundant Condition", self.redundantCondition
        print "Insufficient Genes", self.insufficientGenes
        
        if reset:
            self.noHeadWidth = 0
            self.noHeadLink = 0
            self.noTailLink = 0
            self.redundantCondition = 0
            self.insufficientGenes = 0
    
    def isNested(self, width, index):
        """Marks the bicluster at index of width conditions if it is nested in
        another bicluster
        
        A Bicluster is considered nested if its genes are a subset of the
        enclosing Bicluster's genes and the Bicluster's conditions are an
        ordered subset of the enclosing Biclusters conditions.
        @param bicluster Bicluster to test for nested-ness
        @return True if bicluster is nested; False otherwise.
        """
        
        if width not in self.cache:
            return False
        innerGroup = self.cache[width]
        
        # if already marked
        if innerGroup.nested[index] == NESTED.nested:
            return True
        elif innerGroup.nested[index] == NESTED.nonnested:
            return False
        
        if width + 1 not in self.cache:
            innerGroup.nested[index] = NESTED.nonnested
            return False
        outerGroup = self.cache[width]
        
        genes = innerGroup.genes[index]
        conditions = innerGroup.conditions[index]
        
        for outer in xrange(outerGroup.depth()):
            # if nested genes are a subset
            if (genes.issubset(outerGroup.genes[outer]) and
                conditions.isOrderedSubset(outerGroup.conditions[outer])):
                # nested-ness is a short-circuited 'or' attribute, so as so soon
                # as one enclosing bicluster is found function can exit
                innerGroup.nested[index] = NESTED.nested
                return True
        
        # bicluster can only be marked as nonnested after all possible
        # enclosing biclusters are checked
        innerGroup.nested[index] = NESTED.nonnested
        return False
    
    def depth(self, width, includeNested=True):
        """Returns number of biclusters of width conditions
        
        @param width number of conditions in biclusters
        @param includeNested True to include all biclusters; False to only
               include biclusters not marked as nested
        """
        if width not in self.cache:
            return 0
        
        return self.cache[width].depth(includeNested)
    
    def __str__(self):
        rows = list()
        
        for width in xrange(2, self.maxConditions + 1):
            if width not in self.cache:
                break
            
            rows.append(str(self.cache[width]))
        
        return ''.join(rows)

class WidthGroupCache(object):
    
    SLOTS = 3
    
    def __init__(self, file, parent, maxConditions, maxGenes):
        """Creates a WidthGroup cache with 2 slots
        
        Cache is fully associative
        @param maxConditions maxConditions in biclusters in width groups held in cache
        """
        self.file = file
        self.parent = parent
        self.maxConditions = maxConditions
        self.maxGenes = maxGenes
        
        widthClass = Biclustering.Sizing.sizeArray(self.maxConditions)
        self.widths = numpy.zeros(self.SLOTS,
                                  dtype = widthClass)
        self.groups = numpy.zeros(self.SLOTS, dtype = WidthGroup)
        ageClass = Biclustering.Sizing.sizeArray(self.SLOTS)
        self.ages = numpy.arange(self.SLOTS, dtype = ageClass)
    
    def __contains__(self, width):
        return widthGroupName(width) in self.parent
    
    def updateAges(self, slot):
        if self.ages[slot] != 0:
            agedSlots = numpy.where(self.ages < self.ages[slot])
            self.ages[agedSlots] += 1
            self.ages[slot] = 0
    
    def __getitem__(self, width):
        slot = numpy.where(self.widths == width)
        # if single slot returned
        if len(slot) == 1:
            slot = slot[0]
            group = self.groups[slot]
            
            self.updateAges(slot)
        else:
            if width == 2:
                group = SeedGroup(self.file, self.parent, self.maxConditions,
                                  self.maxGenes)
            else:
                group = WidthGroup(self.file, self.parent, self.maxConditions,
                                   self.maxGenes, width)
            
            slot = numpy.where(self.ages == self.SLOTS - 1)
            self.groups[slot] = group
            self.widths[slot] = width
            
            self.updateAges(slot)
        
        return group

class BiclusterTableAccessor(object):
    
    def __init__(self, maxGenes, maxConditions, width):
        self.conditions = Biclustering.Bit.OrderedBitSetAccessor('conditions',
                                                                 width,
                                                                 maxConditions)
        self.genes = Biclustering.Bit.BitSetAccessor('genes', maxGenes)
        
        class BiclusterTable(tables.IsDescription):
            conditions = self.conditions.description
            genes = self.genes.description
            nested = tables.EnumCol(NESTED, 'unknown')
        
        BiclusterTable._v_flavor = 'numpy'
        self.description = BiclusterTable
    
    def unpack(self, row, name):
        """Unpacks data in Row row in col name
        
        @param row Row packed data is in
        @param name name of column packed data is in
        """
        try:
            unpacked = getattr(self, name).unpack(row)
        except AttributeError:
            unpacked = row[name]
        
        return unpacked
    
    def pack(self, row, name, data):
        """Packs data into row in col name
        
        @param row Row to pack data into
        @param name name of column data should be packed into
        @param data data to be packed
        """
        try:
            getattr(self, name).pack(row, data)
        except AttributeError:
            row[name] = data

def widthGroupName(width):
    return "width" + str(width)

class WidthGroup(object):

    def __init__(self, file, parent, maxConditions, maxGenes, width):
        """Returns group for storing bicluster of width conditions."""
        self.file = file
        self.maxConditions = maxConditions
        self.maxGenes = maxGenes
        self.width = width
        
        name = widthGroupName(width)
        
        try:
            self.group = file.getNode(parent, name)
        except tables.NoSuchNodeError:
            self.group = file.createGroup(parent, name)
        
        self.biclusterPoolAccessor = BiclusterTableAccessor(maxGenes,
                                                            maxConditions,
                                                            width)
        self.biclusterPool = \
            file.createTable(self.group, "biclusterPool",
                             self.biclusterPoolAccessor.dictDescription)
    
    def __getitem__(self, link):
        return self.group.heads[link]
    
    def pool(self, conditions, genes):
        row = self.biclusterPool.row
        
        self.biclusterPoolAccessor.pack(row, 'conditions', conditions)
        self.biclusterPoolAccessor.pack(row, 'genes', genes)
        
        row.append()
    
    def indexPosition(self, name, position):
        atomClass = Biclustering.Sizing.sizeAtom(self.maxConditions)
        atom = atomClass(flavor = 'numpy')
        expectedSizeInMB = atom.atomsize() * self.nested.nrows / float(1 << 20)
        indexArray = \
            self.file.createVLArray(self.group, name, atom,
                                    expectedsizeinMB = expectedSizeInMB)
        self.file.setNodeAttr(indexArray, "position", -1)
        
        progressBar = \
            Biclustering.Timing.ProgressBar(self.maxConditions, name)
        
        for i in xrange(self.maxConditions):
            progressBar.update()
            
            entry = self.conditions.where(position, i)
            indexArray.append(entry)
        
        progressBar.finish()
        indexArray.flush()
    
    def index(self):
        self.indexPosition("heads", -1)
    
    def __str__(self):
        rows = list()
        first = True
        for row in self.biclusterPool.iterrows():
            if first:
                first = False
            else:
                rows.append("\n")
            
            conditions = self.biclusterPoolAccessor.unpack(row, 'conditions')
            genes = self.biclusterPoolAccessor.unpack(row, 'genes')
            rows.append("Conditions: %s Genes: %s" % (conditions, genes))
        
        return ''.join(rows)
    
    def depth(self, includeNested=True):
        """Returns number of biclusters of width conditions
        
        @param width number of conditions in biclusters
        @param includeNested True to include all biclusters; False to only
               include biclusters not marked as nested
        """
        
        if includeNested:
            count = self.nested.nrows
        else:
            count = len(numpy.where(self.nested[:] != NESTED.nested))
        
        return count

class SeedGroup(WidthGroup):
    
    def __init__(self, file, parent, maxConditions, maxGenes):
        super(SeedGroup, self).__init__(file, parent, maxConditions, maxGenes, 2)
    
    def __getitem__(self, link):
        if isinstance(link, tuple):
            if link[0] == "tails":
                return self.group.tails[link[1]]
            else:
                link = link[1]
        
        return super(SeedGroup, self).__getitem__(link)
    
    def index(self):
        super(SeedGroup, self).index()
        self.indexPosition("tails", 0)
