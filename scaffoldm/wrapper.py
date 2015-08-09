#!/usr/bin/env python
###############################################################################
#                                                                             #
#    wrapper.py                                                               #
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
#This files is the overall wrapper for the project
#It will call all the classes, dataloader, dataparser and scaffold  
#(Will maybe make a contig class)
#Progess is expected as follows
#dataloader is called - passes data to data parsers
#data parser calls gapest and algo to create the scaffold dicitonary
#Scaffold class is called and created for each scaffold
#Handles commandline input and python package organisation
###############################################################################
import os
import argparse
###############################################################################
if __name__ == "__main__": ###Check if arguments coming in from command line
    
    parser = argparse.ArgumentParser(description='CreateNewContigs.')

    parser.add_argument('-N','--name', type=str, nargs='?', \
		help='The name of the file',default='MG1655ref.fasta')

    parser.add_argument('-P','--path', type=str, nargs='?', \
	default='ReferenceFasta/', \
		help='The name of absolute path') #old abspath /home/baker/Documents/Geneslab/TestDataset/

    parser.add_argument('-L','--list', type=int, nargs='*', help='Test',default=[1000,5000,1,4000,20000,1])


    args = parser.parse_args()
    print(args.path+args.name,end="\n",sep='')
    print(args.list,end="\n",sep='')

    data=DataLoader(bamnames,
                 contigloc,
                 libno=something,
                 linksname=linkfile,
                 covname=covfile,
                 insertname=insertfile,
):)
    ###Note dataparser will have inside it all the scaffolds made from the data.
    parser=DataParser(  data.links,
                            data.coverages,
                            data.inserts,
                            data.bamnames,
                            data.contigNames)
    parser.parse()
