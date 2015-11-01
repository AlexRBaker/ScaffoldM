#!/usr/bin/env python
###############################################################################
#                                                                             #
#    dataloader.py                                                            #
#                                                                             #
#    A class for loading bamm files and contigs (FASTA)                       #
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
import os
import sys
import numpy as np
import datetime
import re #Regular expressions for name matching in contigs
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Few goals in this file:
# Load bamms - use BamM to give cov and link information
# be able to extract links and coverage
# Check at least some contigs in contig file are present in bam
# if not abort - indicate they don't match
# 
###############################################################################
class DataLoader(object):
    """Utilities wrapper"""
    def __init__(self,
                 bamnames,
                 contigloc='Dupes.fna',
                 cov=True,
                 links=True,
                 useBamm=True,
                 linksname='links.tsv',
                 covname='covs.tsv',
                 insertname='inserts.tsv',
                 libno='1'):
    ####
    # Variables:
    # coverage - binary - store coverage info
    # links - binary - store links info
    # bammloc - location of bamm file/s
    # contig loc - location of contig files
    ####
    # Attributes
    # One for all var (same name)
    # Storage of Bamm processing
    # Extraction of useful Bamm bits
    # Some functions for getting info
    ####
        if cov:
            self.cov=1
        else:
            self.cov=0
        if links:
            self.links=1
        else:
            self.links=0
        if os.path.isfile(covname) and os.path.isfile(covname.split(".tsv")[0]+"var.tsv"):
            self.cov=0
        if os.path.isfile(linksname):
            self.links=0
        self.bamnames=bamnames
        self.contigloc=contigloc
        self.contigdict=None #Initialise as None
        try:
            if len(libno)==len(bamnames):
                pass
            elif len(libno)==1:
                libno=libno*len(bamnames)
            else:
                pass
        except TypeError or NameError:
            try:
                libno=[str(libno)]*len(bamnames)
            except ValueError:
                print "libno must be integer, a string of an integer or a list of such"
        if useBamm:
            #self.bammparse=BP.parseBams(bammloc,doLinks=self.links, /
            #doCovs=self.cov,threads=min(len(bammloc),CompThreads)
            #BP.printLinks('links.tsv')
            sbamnames=" ".join(bamnames) #assumes that it ends in .bam
            #Alternative - continue using as CLI tool since python lib is bugged
            if self.cov==1 and self.links==1:
                os.system("bamm parse -n {0} -b {1} -i {2} -l {3} -c {4} -m pvariance".format(' '.join(libno),sbamnames,insertname,linksname,covname.split(".tsv")[0]+"var.tsv"))
                os.system("bamm parse -b {0} -c {1} -m pmean".format(sbamnames,covname))
            elif self.links==1:
                os.system("bamm parse -n {0} -b {1} -i {2} -l {3}".format(' '.join(libno),sbamnames,insertname,linksname))
            elif self.cov==1:
                os.system("bamm parse -n {0} -b {1} -i {2} -c {3} -m pvariance".format(' '.join(libno),sbamnames,insertname,covname.split(".tsv")[0]+"var.tsv"))
                os.system("bamm parse -b {0} -c {1} -m pmean".format(sbamnames,covname))
            else:
                print "Please Ensure you have provided appropiate input"
        else:
            print 'Error using bamm'
        #[1:] to remove header
        #print os.getcwd()
        #np.loadtxt autoignores the header
        print "step 1"
        print datetime.datetime.now().time()
        self.contigNames=list(set(self.parsetsv(covname)[:,0])) #First column of coverage file has all contigs
        print "step2"
        print datetime.datetime.now().time()
        self.links=self.parsetsv(linksname) #Get all links
        print "step3"
        print datetime.datetime.now().time()
        self.inserts=self.parsetsv(insertname)
        print "step4"
        print datetime.datetime.now().time()
        self.getcovs(covname) #Extract coverage mean and variance - defines self.coverages as
        #dict[bamname][tigname]=[coveragemean,coveragesd]
        print "step5"
        print datetime.datetime.now().time()
        self.extractcontigs(self.contigNames,self.contigloc)
        print "Contigs loaderd"
        print datetime.datetime.now().time()

    def parsetsv(self,textfile="links.tsv",delim='\t',comment='#'):
        ###Need to extend to attempt to detect delimiter
        ##Should be existing package for this
        import os
        import csv
        import numpy as np
        try:
            return np.loadtxt(textfile, dtype=str, comments=comment, delimiter=delim)
        except:
            print "This is the current directory", os.getcwd()
            
    def getcovs(self,covname):
        '''Opens variance and mean coverage file. Returns a dictionary with the mean and std of each
        contig in each bam file'''
        tig=(0,0)
        bam=(0,0)
        try:
            means=self.parsetsv(covname,comment='|#$|') #comment is just
            #to stop deafult removal of header since header starts with #
            N_bams=len(means[1,:])-2 #Excludes contig name and length columns
            N_tigs=len(means[:,0])-1 #Number of contigs
            variance=self.parsetsv(covname.split(".tsv")[0]+"var.tsv",comment='|#$|')
            Nam_Bam1=means[2:,0] #The names of the bam
            Tigs= means[1:,0]#The names of the contigs
            coverage={} #coverage value
            #Assuming bamnames are in same order - check later - looks to always be case
            #Create a dictionary for each bam file with a mean,sd entry for each 
            coverages={(Nam_Bam1[bam].split(".bam")[0]):\
            {(Tigs[tig]):[means[tig+1,bam+2].astype(float),np.sqrt(variance[tig+1,bam+2].astype(float))] \
            for tig in range(N_tigs)} for bam in range(N_bams)}
            self.coverages=coverages
            tiglen={Tigs[tig]:int(means[tig+1,1]) for tig in range(N_tigs)}
            self.tiglen=tiglen
        except IndexError:
            print "The index was out of bounds"
            print tig, "Contig index"
            print bam, "Bam Index"
            print coverage
        except:
            print "Unidentified Error occurred"
            raise

    def extractcontigs(self,contignames,contigloc,header=True):
        '''Just assigns contigs file via contigloc.
        Temporary just for use when making scaffold
        Will extract the text for that contig'''
        #~ print "Starting contig extraction"
        if isinstance(self.contigdict,type(None)):
            curname=None
            try:
                if isinstance(contignames,list):
                    re_contignames=[re.compile('{0}[^\d]'.format(contigname)) for contigname in contignames] #Excludes decimal follow on in name
                    #I.E contig1 will not match contig11 in the search for the contig name in the header
                    #Any extra char after the contig a int (since other contig names might conflict)
                    self.contigdict={}
                    with open(contigloc,'r') as Contigs:
                        head=Contigs.readline()
                        if not head.startswith('>'):
                            print "An error in file format"
                            raise TypeError("Not a FASTA file:")
                        Contigs.seek(0)
                        Go=False
                        for line in Contigs:
                            if line.startswith('>'):
                                if any(not isinstance(regexp.search(line),type(None)) for regexp in re_contignames):
                                    curname=self.findname(contignames,line)
                                    self.contigdict[curname]=[]
                                    Go=True
                                else:
                                    Go=False
                            elif not isinstance(curname,type(None)) and not line.startswith('>') and Go:
                                self.contigdict[curname].append(line.translate(None,'\n')) #Adds sequence line to current contig
                                #Also removes all line breaks
                            else:
                                pass
                    for contig in self.contigdict.keys():
                        self.contigdict[contig]=''.join(self.contigdict[contig]) #Turns into large string
                else:
                    print "List of contigs error"
                    raise TypeError("Need list of contigs to form contigdict")
            except:
                print "Another error"
                print self.contigNames, "These are the names of the contigs"
                raise
        else:
            try:
                return self.contigdict[contignames]
            except:
                print "Probably a key error"
                print contignames
                print self.contigdict
                raise
                
    def findname(self,contignames,line):
        reg_exp=[re.compile('{0}[^\d]'.format(contigname)) for contigname in contignames]
        finalname=None
        for i,contigname in enumerate(reg_exp):
            if not isinstance(contigname.search(line),type(None)):
                #print contigname
                finalname=contignames[i]
            else:
                pass
        return finalname
        
