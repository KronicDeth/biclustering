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
"""Timing utils.

Times functions and statements for average execution time
Adapted from public domain code released on ASPN

@author Luke Imhoff
@license GPLv2
"""

import datetime
import time

class ProgressBar(object): 
    """Text Progress bar for console apps with ETF
    
    @author Luke Imhoff
    """
    
    def __init__(self, maxCount, title="", updateInterval=10):
        self.title = title
        self.maxCount = maxCount
        self.count = 0
        self.start = None
        
        self.updateInterval = updateInterval
        self.nextUpdate = 0
        
    def update(self):
        if self.start is None:
            self.start = time.time()
            logging.info("  Starting at", time.asctime())
            return False
        
        self.count += 1
        
        if self.count > self.maxCount:
            self.count = self.maxCount
        
        elapsed = time.time() - self.start
        
        # if enough time has elapsed
        if elapsed >= self.nextUpdate:
            
            self.nextUpdate = elapsed + self.updateInterval
            
            estimatedTotalTime = self.maxCount * elapsed / self.count
            estimatedTimeToFinish = estimatedTotalTime - elapsed
            etfText = datetime.timedelta(seconds = estimatedTimeToFinish)
            
            logging.info("%s: %d / %d %s", self.title, self.count,
                         self.maxCount,  etfText)
            return True
        else:
            return False
    
    def finish(self):
        if self.start is None:
            elapsed = 0
        else:
            elapsed = time.time() - self.start
        
        logging.info("  Ended at", time.asctime())
        logging.info("%s: %d Completed %s",
                     self.title, self.count + 1,
                     datetime.timedelta(seconds = elapsed))
