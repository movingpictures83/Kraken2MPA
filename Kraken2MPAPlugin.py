#!/usr/bin/env python
####################################################################
#kreport2mpa.py converts a Kraken-style report into mpa [MetaPhlAn) format
#Copyright (C) 2017-2020 Jennifer Lu, jennifer.lu717@gmail.com

#This file is part of KrakenTools.
#KrakenTools is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

####################################################################
#Jennifer Lu, jlu26@jhmi.edu
#11/06/2017
#Updated: 07/12/2020
#
#This program reads in a Kraken report file and generates
#an mpa-format (MetaPhlAn) style report. Each line represents
#a possible taxon classification. The first column is lists the 
#domain, kingdom, phyla, etc, leading up to each taxon.
#The levels are separated by the | delimiter, with the type of 
#level specified before each name with a single letter and underscore
#(d_ for domain, k_ for kingdom, etc). 
#The second column is the number of reads classified within 
#that taxon's subtree.
#
#Input file:
#   - Kraken report file generates from the kraken raw output file
#Input Parameters to Specify [OPTIONAL]:
#   - header_line = prints a header line in mpa-report 
#       [Default: no header]
#   - intermediate-ranks = includes non-traditional taxon levels
#       (traditional levels: domain, kingdom, phylum, class, order, 
#       family, genus, species)
#       [Default: no intermediate ranks]
#Output file format (tab-delimited)
#   - Taxonomy tree levels |-delimited, with level type [d,k,p,c,o,f,g,s,x]
#   - Number of reads within subtree of the specified level
#
#Methods
#   - main
#   - process_kraken_report
#
import os, sys, argparse
import PyIO
import PyPluMA

#process_kraken_report
#usage: parses a single line in the kraken report and extracts relevant information
#input: kraken report file with the following tab delimited lines
#   - percent of total reads
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level
#       (U, D, P, C, O, F, G, S, -)
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria,...etc)
#   - spaces + name
#returns:
#   - classification/genome name
#   - level name (U, -, D, P, C, O, F, G, S)
#   - reads classified at this level and below in the tree
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    if len(split_str) < 4:
        return []
    try:
        int(split_str[1])
    except ValueError:
        return []
    percents = float(split_str[0])
    all_reads = int(split_str[1])
    #Extract relevant information
    try:
        taxid = int(split_str[-3]) 
        level_type = split_str[-2]
        map_kuniq = {'species':'S', 'genus':'G','family':'F',
            'order':'O','class':'C','phylum':'P','superkingdom':'D',
            'kingdom':'K'}
        if level_type not in map_kuniq:
            level_type = '-'
        else:
            level_type = map_kuniq[level_type]
    except ValueError:
        taxid = int(split_str[-2])
        level_type = split_str[-3]
    #Get name and spaces 
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break
    name = name.replace(' ','_')
    #Determine level based on number of spaces
    level_num = spaces/2
    return [name, level_num, level_type, all_reads, percents]

class Kraken2MPAPlugin:
   def input(self, infile):
        self.parameters = PyIO.readParameters(infile)
        if ("report" in self.parameters):
            self.r_file = PyPluMA.prefix()+"/"+self.parameters["report"]
        if ("display_header" in self.parameters):
            self.add_header = self.parameters["display_header"]
        else:
            self.add_header = "False"
        if ("read_count" in self.parameters):
            self.use_reads = self.parameters["read_count"]
        else:
            self.use_reads = "True"
        if ("intermediate_ranks" in self.parameters):
            self.x_include = self.parameters["intermediate_ranks"]
        else:
            self.x_include = "False"


   def run(self):
       pass

   def output(self, outfile):

    #Process report file and output 
    curr_path = [] 
    prev_lvl_num = -1
    r_file = open(self.r_file, 'r')
    o_file = open(outfile, 'w')
    #Print header
    if self.add_header == "True":
        o_file.write("#Classification\t" + os.path.basename(self.r_file) + "\n")
    
    #Read through report file 
    main_lvls = ['R','K','D','P','C','O','F','G','S']
    for line in r_file:
        report_vals = process_kraken_report(line)
        #If header line, skip
        if len(report_vals) < 5: 
            continue
        #Get relevant information from the line 
        [name, level_num, level_type, all_reads, percents] = report_vals
        if level_type == 'U':
            continue
        #Create level name 
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "K":
            level_type = "k"
        elif level_type == "D":
            level_type = "k"
        level_str = level_type.lower() + "__" + name
        #Determine full string to add
        if prev_lvl_num == -1:
            #First level
            prev_lvl_num = level_num
            curr_path.append(level_str)
        else:
            #Move back if needed
            while level_num != (prev_lvl_num + 1):
                prev_lvl_num -= 1
                curr_path.pop()
            #Print if at non-traditional level and that is requested
            if (level_type == "x" and self.x_include == "True") or level_type != "x":
                #Print all ancestors of current level followed by |
                for string in curr_path:
                    if (string[0] == "x" and self.x_include == "True") or string[0] != "x":
                        if string[0] != "r": 
                            o_file.write(string + "|")
                #Print final level and then number of reads
                if self.use_reads == "True":
                    o_file.write(level_str + "\t" + str(all_reads) + "\n")
                else:
                    o_file.write(level_str + "\t" + str(percents) + "\n")
            #Update
            curr_path.append(level_str)
            prev_lvl_num = level_num
    o_file.close()
    r_file.close()

