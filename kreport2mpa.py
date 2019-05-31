#! /usr/bin/env python
####################################################################
#kreport2mpa.py converts a Kraken-style report into mpa [MetaPhlAn) format
#Copyright (C) 2019 Jennifer Lu, jlu26@jhmi.edu

#This file is part of Kraken-Tools.
#Kraken-Tools is free software; you can redistribute it and/or modify
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
#Updated: 03/19/2019
#
#This program reads in a Kraken report file and generates
#an mpa-format (MetaPhlAn) style report. 
#
#Parameters:
#   -h, --help...............show help message
#   -r X, --report-file X....input kraken report filename 
#   -o X, --output X.........output mpa-report filename 
#   --display-header.........include header in mpa-report file [default:false]
#   --intermediate-ranks.....include non-traditional taxonomic ranks 
#   --no-intermediate-ranks..do not include non-traditional ranks [default] 
#Input file format (tab-delimited)
#   - percentage of total reads
#   - number of reads (including reads within subtree)
#   - number of reads (only at this level)
#   - taxonomic classification of level (U, D, P, C, O, F, G, S,...etc) 
#   - NCBI taxonomic ID 
#   - name of level   
#Output file format (tab-delimited)
#   - Taxonomy tree levels |-delimited, with level type [k,p,c,o,f,g,s,x]
#   - Number of reads within subtree of the specified level
#
#Methods
#   - main
#   - process_kraken_report
####################################################################
import os, sys, argparse

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
    try:
        int(split_str[1])
    except ValueError:
        return []
    all_reads = int(split_str[1])
    level_type = split_str[3]
    #Get name and spaces 
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break
    #Determine level based on number of spaces
    level_num = spaces/2
    return [name, level_num, level_type, all_reads]

#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--report-file', '--report', required=True,
        dest='r_file', help='Input kraken report file for converting')
    parser.add_argument('-o', '--output', required=True,
        dest='o_file', help='Output mpa-report file name')
    parser.add_argument('--display-header', action='store_true', 
        dest='add_header', default=False, required=False,
        help='Include header in mpa-report file [default: false]') 
    parser.add_argument('--intermediate-ranks', action='store_true',
        dest='x_include', default=False, required=False,
        help='Include non-traditional taxonomic ranks in output')
    parser.add_argument('--no-intermediate-ranks', action='store_false',
        dest='x_include', default=False, required=False,
        help='Do not include non-traditional taxonomic ranks in output [default: no intermediate ranks]')
    args=parser.parse_args()

    #Process report file and output 
    curr_path = [] 
    prev_lvl_num = -1
    r_file = open(args.r_file, 'r')
    o_file = open(args.o_file, 'w')
    #Print header
    if args.add_header:
        o_file.write("#Clade\tReads\n")
    
    #Read through report file 
    main_lvls = ['D','P','C','O','F','G','S']
    for line in r_file:
        report_vals = process_kraken_report(line)
        #If header line, skip
        if len(report_vals) < 4: 
            continue
        #Get relevant information from the line 
        [name, level_num, level_type, all_reads] = report_vals
        if level_type == 'U':
            continue
        #Create level name 
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "D":
            level_type = "K"
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
            if (level_type == "x" and args.x_include) or level_type != "x":
                #Print all ancestors of current level followed by |
                for string in curr_path:
                    if (string[0] == "x" and args.x_include) or string[0] != "x":
                        if string[0] != "r": 
                            o_file.write(string + "|")
                #Print final level and then number of reads
                o_file.write(level_str + "\t" + str(all_reads) + "\n")
            #Update
            curr_path.append(level_str)
            prev_lvl_num = level_num
    o_file.close()
    r_file.close()

if __name__ == "__main__":
    main()
