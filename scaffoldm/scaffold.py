#!/usr/bin/env python
###############################################################################
#                                                                             #
#    sacffold  .py                                                            #
#                                                                             #
#    Class for storing and printing scaffolds and assoicated summary info     #
#                                                                             #
#    Copyright (C) Alexander Baker                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Alexander Baker"
__copyright__ = "Copyright 2015"
__credits__ = ["Alexander Baker"]
__license__ = "GPLv3"
__maintainer__ = "Alexander Baker"
__email__ = "Alexander.baker@uqconnect.edu.au"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
#Overall scaffold structre
#Include how to print scaffold
#Will have dict detailing necessary info on scaffolds
#involved contig and scaffold names will be extracted
# Key - Tuple - contig and scaffold
# item - start point and end point in contig, orientation, downstream contig
###############################################################################
class Scaffold(object):
    """Utilities wrapper"""
    def __init__(self,
                 cov=False,
                 links=True,
                 bammloc,
                 contigloc):

        
	self.Scaffold={}
    def extractcontigs(self,contigname,contigloc):
    ###Just assigns contigs file via contigloc.
	###Temporary just for use when making scaffold
    ###Will extract the text for that contig
        import sys
        try:
            with open(contigfile,'r') as Contigs:
                head=Contigs.readline()
                if not head.startswith('>'):
                    raise TypeError("Not a FASTA file:")
                Contigs.seek(0)
                title=head[1:].rstrip() ##Strips whitespace and >
                record=0
                contigseq=[]
                for line in Contigs:
                    if line.startswith('>') and line.find(contigname)>=0:
                        record=1
                        contigseq.append(line)
                    elif line.startswith('>') and contigname not in line:
                        record=0
                    elif record==1:
                        contigseq.append(line)
                    else:
                        pass
                seq=''.join(contigseq)
                return seq
        except:
            print "Error opening file:", contigfile,sys.exc_info()[0]
            raise
	
	self.contignames=self.scaffold.keys()
	#all_keys = set().union(*(d.keys() for d in mylist)) gets all
	# from list of dicts
    def extname(pathname):
		
		
		return nam
        
        
	def printscaffold(self):
		
		
		
		
		
		
    def runCommand(self, cmd):
        """Run a command and take care of stdout

        expects 'cmd' to be a string like "foo -b ar"

        returns (stdout, stderr)
        """
        from multiprocessing import Pool
        from subprocess import Popen, PIPE

        p = Popen(cmd.split(' '), stdout=PIPE)
        return p.communicate()

    def parseFile(self, filename):
        """parse a file"""
        import sys
        try:
            with open(filename, "r") as fh:
                for line in fh:
                    print line
        except:
            print "Error opening file:", filename, sys.exc_info()[0]
            raise 
