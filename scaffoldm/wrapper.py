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
dataloader=DataLoader()
dataparser=DataParser()
