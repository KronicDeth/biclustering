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
"""Gene Expression Matrix

@author Luke Imhoff
@license GPLv2
"""

import datetime
import numpy
import tables
import time

import Biclustering.Array
import Biclustering.Bicluster
import Biclustering.Bit
import Biclustering.BitSet
import Biclustering.Combinatorics
import Biclustering.Sizing
import Biclustering.Timing

def fullCoverageData(conditions, minGenes=2):
    """Returns a GEM containing all possible biclusters of [2,conditions]
    conditions with minGenes genes
    
    @param conditions number of conditions in GEM
    @param minGenes genes in each bicluster of conditions conditions in GEM
    @return data for a GEM
    """
    patterns = Biclustering.Combinatorics.xpermutations(conditions, conditions)
    
    data = numpy.zeros((2 * patterns.len(), conditions))
    
    offset = 0
    for p in patterns:
        data[offset : offset + 2] = p
        offset += 2
    
    return data

def packData(data):
    """Compacts data by converting it to the minimum int size
    
    @param data data to compact
    @return compacted data
    """
    
    compact = Biclustering.Sizing.sizeArray(data.shape[1])
    
    return data.argsort().argsort().astype(compact)

class GeneExpressionMatrix(object):
    
    FILE_EXTENSION = "gem"
    FILTERS = tables.Filters(complevel = 1, complib= 'lzo')
    
    def __init__(self, name, data=None, path=None, minGenes=2, filters=None):
        self.name = name
        
        if path is None:
            # default to PWD
            path = "./"
            
        fileName = path + name + "." + GeneExpressionMatrix.FILE_EXTENSION
        
        if filters is None:
            filters = self.FILTERS
        
        # if creating this GEM
        if data is not None:
            self.file = tables.openFile(fileName, mode = "w",
                                        title = name,
                                        filters = filters)
            
            group = self.file.createGroup("/", "gem")
            # save raw version in case it's needed for algorithm enhancements
            self.file.createArray(group, "raw", data)
            self.data = packData(data)
            
            createBiclusters = True
        else:
            self.file = tables.openFile(fileName, mode = "r+",
                                        filters = filters)
            
            group = self.file.getNode("/", "gem")
            raw = self.file.getNode(group, "raw")
            self.data = packData(raw[:])
            
            createBiclusters = False
        
        self.maxConditions = self.data.shape[1]
        self.maxGenes = self.data.shape[0]
        
        self.biclusters = Biclustering.Bicluster.Group(self.file, "/",
                                                       self.maxConditions,
                                                       self.maxGenes,
                                                       createBiclusters,
                                                       minGenes)
    
    def splitSubset(self, conditions):
        """Identifies all biclusters with a given subset of 2 conditions
        
        @param conditions subset of conditions in matrix to search for biclusters
        @return number of valid biclusters found
        """
        
        if conditions.size != 2:
            raise ValueError("conditions subset can only have 2 conditions for split")
        
        increasing = numpy.where(self.data[:, conditions[0]] < self.data[:, conditions[1]])
        increasingGenes = Biclustering.BitSet.BitSet(self.data.shape[0], increasing)
        
        count = 0
        # increasing set
        increasingConditions = Biclustering.Bit.OrderedBitSet(conditions, self.data.shape[1])
        if self.biclusters.pool(increasingConditions, increasingGenes):
            count += 1
        
        # decreasing set
        
        decreasingConditions = increasingConditions.reverse()
        if self.biclusters.pool(decreasingConditions, ~increasingGenes):
            count += 1
        
        return count
    
    def splitBiclusters(self):
        """Finds all biclusters with 2 conditions"""
        
        combinations = Biclustering.Combinatorics.xcombinations(self.maxConditions, 2)
        
        progressBar = \
            Biclustering.Timing.ProgressBar(combinations.len(),
                                                "Splitting")
        
        count = 0
        for conditions in combinations:
            progressBar.update()
            
            count += self.splitSubset(conditions)
        
        progressBar.finish()
        self.file.flush()
        
        return count
    
    def indexBiclusters(self, width):
        """Indexes all biclusters
        
        Indexes are used to speed up chain*() methods
        @param width width of biclusters to index.  chain*(width) cannot be
                     called before calling indexBiclusters(width)
        """
        self.biclusters.index(width)
    
    def chainBiclusters(self, tailWidth):
        """Chains biclusters into larger biclusters
        
        Chained biclusters are formed by chaining one bicluster of tailWidth
        with a bicluster of width 2.  
        @param tailWidth number of conditions in first array of biclusters
        @return number of biclusters found
        """
        
        title = "(%d %d) => (%d)" % (2, tailWidth,
                                     tailWidth + 1)
        progressBar = \
            Biclustering.Timing.ProgressBar(self.maxConditions, title)
        
        count = 0
        for link in xrange(self.maxConditions):
            progressBar.update()
            
            count += self.biclusters.chain(tailWidth, link)
        
        progressBar.finish()
        self.file.flush()
        
        return count
    
    
    def chainBiclustersPreCrest(self, headWidth, doubling = False):
        """Chains biclusters into larger biclusters
        
        Chained biclusters are formed by chaining one bicluster of headWidth
        with a bicluster of width 2.  
        @param headWidth number of conditions in first array of biclusters
        @param doubling 
        @return number of biclusters found
        """
        
        if doubling:
            tailWidth = headWidth
        else:
            tailWidth = 2
        
        title = "(%d %d) => (%d)" % (headWidth, tailWidth,
                                     headWidth + tailWidth - 1)
        progressBar = \
            Biclustering.Timing.ProgressBar(self.maxConditions, title)
        
        count = 0
        for link in xrange(self.maxConditions):
            progressBar.update()
            
            count += self.biclusters.chainPreCrest(headWidth, link, doubling)
        
        progressBar.finish()
        self.file.flush()
        
        return count
    
    def pruneBiclusters(self, width):
        """Prunes biclusters of width conditions that are nested
        
        @param width number of conditions in biclusters to prune
        """
        
        indexes = xrange(self.biclusters.depth(width))
        
        title = "(%d not in %d)" % (width, width + 1)
        progressBar = Biclustering.Timing.ProgressBar(len(indexes), title)
        
        for index in indexes:
            progressBar.update()
            
            self.biclusters.isNested(width, index)
        
        progressBar.finish()
        self.file.flush()
    
    def biclusterCount(self, includeNested=True):
        """Returns total number of biclusters
        
        @param includeNested True to include nested biclusters.
                             False otherwise.
        @return total number of biclusters
        """
        
        count = 0
        for i in xrange(2, self.maxConditions + 1):
            count += self.biclusters.depth(i, includeNested)
        
        return count
    
    def allBiclusters(self):
        """Finds all biclusters in the GEM"""
        
        totalStartTime = time.time()
        
        # seed clusters need 2 conditions so biclusters
        # can be grown by 1 condition if needed
        # 0 biclusters is unlikely, but may occur to too high of minGenes
        if self.splitBiclusters() == 0:
            logging.error("No seed biclusters found.  "
                          "Perhaps minimum genes (%d) is too high?",
                          self.minGenes)
        
        logging.info("Chaining")
        
        # search for valid bicluster with most conditions
        maxConditions = self.maxConditions
        # only look for holes above minimum valid bicluster conditions and
        # smaller than the known maxConditions that may still yield genes
        progressBar = \
            Biclustering.Timing.ProgressBar(maxConditions - 2, "Chaining")
        for i in xrange(2, maxConditions + 1):
            progressBar.update()
            
            self.indexBiclusters(i)
            if self.chainBiclusters(i) == 0:
                maxConditions = i
                break
        
        progressBar.finish()
        
        logging.info("Chains constructed. Biclusters: %d Max Conditions: %d",
                     self.biclusterCount(), maxConditions)
        logging.info("Pruning nested Biclusters")
        
        progressBar = \
            Biclustering.Timing.ProgressBar(maxConditions - 2, "Pruning")
        
        for i in xrange(2, maxConditions):
            progressBar.update()
            
            self.pruneBiclusters(i)
        
        progressBar.finish()
        
        logging.info("Nested Biclusters pruned.  Biclusters: %s ",
                     self.biclusterCount(False)) 
        logging.info("Total Time: %s",
                    datetime.timedelta(seconds = time.time() - totalStartTime))
    
    def stats(self):
        """Prints stats on GEM
        
        @param includeNested True to count nested (pruned) bicluster.  False to
               ignore them
        """
        
        lines = list()
        for i in xrange(2, self.maxConditions + 1):
            if i != 2:
                lines.append("\n")
            lines.append("(%d): %d T %d NSUB" %
                         (i, self.biclusters.depth(i),
                          self.biclusters.depth(i, False)))
        
        lines.append("\n")
        lines.append("Total: %d NSUB: %d" % 
                      (self.biclusterCount(), self.biclusterCount(False)))
        
        return ''.join(lines)
    
    def annotate(self, annotation, axis="genes"):
        """
        """
        pass

class Annotation(dict):
    
    def __init__(self, universe=None, file=None, where=None):
        """
        
        @param universe number of elements being annotated
        @param file file to load Annotation from
        @param where Node in file annotation is on
        """
        super(Annotation, self).__init__()
        
        if file is not None and where is not None:
            pass
        elif universe is not None:
            self.universe = universe
        else:
            raise ValueError("Could not load or create Annotation with given arguments")
            
    
    def __setitem__(self, key, value):
        if not isinstance(value, Biclustering.BitSet.BitSet):
            value = Biclustering.BitSet.BitSet(self.universe, value)
        
        super(Annotation, self).__setitem__(key, value)
    
    def categories(self):
        return self.keys()

    def saveTo(self, file, where, name):
        """Creates SetArray which is a pickled version of this annotation
        """
        
        atom = tables.StringAtom(flavor = 'numpy')
        node = file.createVLArray(where, name, atom)
        
