#!/usr/bin/env python
################################################################
#bray_curtis.py takes multiple reports of classifications or
#abundances and calculates a matrix of bray-curtis dissimilarity
#metrics
#Copyright (C) 2019 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of KrakenTools
#KrakenTools is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

#################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Christopher Pockrandt, pockrandt@jhu.edu
#Updated: 08/12/2019
#
#This program takes multiple reports of classifications or
#abundances and calculates a matrix of bray-curtis dissimilarity
#metrics
#
#Parameters:
#   -h, --help................show help message.
#   -i X, --input-files X.....all (at least 2) input files (separated by spaces)
#Input File options (all input files must be of the same format)
#   --single [default]........all samples are within a single tab-delimited file
#                             where the first column is the category and the rest are counts
#   --simple .................input files must be tab-delimited with the first
#                             two columns representing categories and then counts
#   --bracken.................input files are Bracken outputs: col  2=taxonomy, 6=counts
#   --kraken..................input files are Kraken reports: col 5=taxonomy, 3=counts
#   --krona...................input files are Krona format: col 1=counts, 2=taxonomy
#Input options:
#   --level [S, G, etc].......user specifies which level to measure at
#                             (for kraken, krona, or bracken input files)
####################################################################
import os, sys, argparse
import operator
from time import gmtime
from time import strftime
import numpy as np
####################################################################
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input','--input-files',
        '--inputs', required=True,dest='in_files',nargs='+',
        help='Input files (one per community) for which to compare for \
        bray-curtis dissimiliarity metrics')
    parser.add_argument('--type', required=False, default='single', dest='filetype',
        choices=['single','simple','bracken','kreport','kreport2','krona'],
        help='Type of input file[s]: single, simple [tab-delimited, specify --cols], \
            bracken, kreport, kreport2, krona. See docs for details')
    parser.add_argument('--cols','--columns', dest='cols',required=False, default='1,2',
        help='Specify category/counts separated by single comma: cat,counts (1 = first col)')
    parser.add_argument('--level', '-l', dest='lvl', required=False, default='all',choices=['all', 'S', 'G', 'F', 'O'],
        help='For Kraken or Krona files, taxonomy level for which to compare samples. Default: all')
    args=parser.parse_args()

    #################################################
    #Test input files
    in2counts = {}
    if args.filetype == 'single' and len(args.in_files) > 1:
        sys.stderr.write("Please specify only one file for '--type simple'\n")
        exit(1)
    for f in args.in_files:
        if not os.path.isfile(f):
            sys.stderr.write("File %s not found\n" % f)
            exit(1)

    #################################################
    #Determine columns for extracting
    categ_col = -1
    count_col = -1
    if args.filetype in ['single','simple']:
        if ',' not in args.cols:
            sys.stderr.write("Please specify column as 'a,b' where a = column of category, \
            b = column of first count\n")
            exit(1)
        else:
            [categ_col, count_col] = args.cols.split(',')
            if not categ_col.isdigit():
                sys.stderr.write("%s is not an integer\n" % categ_col)
                exit(1)
            elif not count_col.isdigit():
                sys.stderr.write("%s is not an integer\n" % count_col)
                exit(1)
            categ_col = int(categ_col) - 1
            count_col = int(count_col) - 1
    elif args.filetype == "bracken":
        categ_col = 0 # TODO taxid (col 1) does not seem to be properly set
        count_col = 5
        taxlvl_col = 2
    elif args.filetype == "kreport" or args.filetype == "kreport2": # TODO: what about kuniq reports?
        categ_col = 4
        count_col = 2
        taxlvl_col = 3
    elif args.filetype == "krona":
        categ_col = 1
        count_col = 0
    #################################################
    #STEP 1: READ IN SAMPLES
    i2totals = {}
    i2counts = {}
    i2names = {}
    num_samples = 0
    num_categories = 0
    if args.filetype == "single":
        #ALL SAMPLE COUNTS WITHIN A SINGLE FILE
        #sys.stdout.write(">>STEP 1: READING INPUT FILE\n")
        #sys.stdout.write("\t....reading line 1 as header\n")
        header = True
        i_file = open(args.in_files[0],'r')
        for line in i_file:
            l_vals = line.strip().split("\t")
            #Read header
            if header:
                s_count = 0
                for i in range(count_col, len(l_vals)):
                    i2names[s_count] = l_vals[i]
                    i2counts[s_count] = {}
                    i2totals[s_count] = 0
                    s_count += 1
                num_samples = s_count
                header = False
                #sys.stdout.write("\t....reading %i samples\n" % s_count)
            else:
                #Otherwise, save_counts
                s_count = 0
                curr_categ = l_vals[categ_col]
                for i in range(count_col, len(l_vals)):
                    if int(l_vals[i]) > 0:
                        i2totals[s_count] += int(l_vals[i])
                        i2counts[s_count][curr_categ] = int(l_vals[i])
                    s_count += 1
                num_categories += 1
        i_file.close()
        #sys.stdout.write("\t....finished reading counts for %i samples in %i categories\n" % (num_samples,num_categories))
    else: # for braken, kraken, kraken2 and krona
        num_samples = 0
        i2names = {}
        i2totals = {}
        i2counts = {}
        genus = {}

        for f in args.in_files:
            i_file = open(f,'r')
            i2names[num_samples] = f
            i2totals[num_samples] = 0
            i2counts[num_samples] = {}
            genus[num_samples] = {}

            for line in i_file:
                l_vals = line.strip().split("\t")

                # empty line, header line or line is a comment
                if len(l_vals) == 0 or (not l_vals[count_col].isdigit()) or l_vals[0] == '#':
                    continue

                if int(l_vals[count_col]) > 0:
                    if args.filetype == "krona":
                        # we don't know which column will be the genus type (might change for every line)
                        # hence, we iterate over it
                        for i in range(count_col, len(l_vals)):
                            if l_vals[i].startswith(args.lvl.lower() + "__") or args.lvl == "all":
                                categ_col = i
                                tax_name = l_vals[categ_col]
                                i2totals[num_samples] += int(l_vals[count_col])
                                if not tax_name in i2counts[num_samples]:
                                    i2counts[num_samples][tax_name] = 0
                                i2counts[num_samples][tax_name] += int(l_vals[count_col])
                                # sys.stdout.write("%s\t%s\t%i\n" % (f, tax_name, int(l_vals[count_col])))
                    else:
                        if l_vals[taxlvl_col][0] == args.lvl or args.lvl == "all": # TODO: what if it is G1?
                            # TODO: cant do this because of broken bracken files: tax_id = int(l_vals[categ_col])
                            tax_id = l_vals[categ_col]
                            genus[num_samples][tax_id] = l_vals[0]
                            i2totals[num_samples] += int(l_vals[count_col])
                            if tax_id not in i2counts[num_samples]:
                                i2counts[num_samples][tax_id] = 0
                            i2counts[num_samples][tax_id] += int(l_vals[count_col])
                            # sys.stdout.write("%s\t%s\t%i\n" % (f, line, i2counts[num_samples][tax_id]))
            i_file.close()
            num_samples += 1
    #################################################
    #STEP 2: CALCULATE BRAY-CURTIS DISSIMILARITIES
    #sys.stdout.write(">>STEP 2: COMPARING SAMPLES TO CALCULATE DISSIMILARITIES\n")

    bc = np.zeros((num_samples,num_samples))
    for i in range(0,num_samples):
        i_tot = i2totals[i]
        for j in range(i+1, num_samples):
            j_tot = i2totals[j]
            C_ij = 0.0
            for cat in i2counts[i]:
                if cat in i2counts[j]:
                    C_ij += min(i2counts[i][cat], i2counts[j][cat])
            #Calculate bray-curtis dissimilarity
            bc_ij = 1.0 - ((2.0*C_ij)/float(i_tot+j_tot))
            #sys.stdout.write("%s\t%s\t%i\t%i\t%i\n" % (i2names[i], i2names[j], C_ij, i_tot, j_tot))
            bc[i][j] = bc_ij
            bc[j][i] = bc_ij

    # for i in range(0,num_samples):
    #     sys.stdout.write("Totals: %s\t%i\n" % (i2names[i], i2totals[i]))
    #     for cat in i2counts[i]:
    #         sys.stdout.write("%i\t%s\t%i\n" % (cat, genus[i][cat], i2counts[i][cat]))

    #################################################
    #sys.stdout.write(">>STEP 3: PRINTING MATRIX OF BRAY_CURTIS DISSIMILARITIES\n")
    #sys.stdout.flush()
    #Print samples
    for i in i2names:
        sys.stdout.write("#%i\t%s (%i reads)\n" % (i,i2names[i],i2totals[i]))
    #Print headers
    sys.stdout.write("x")
    for i in range(num_samples):
        sys.stdout.write("\t%i" % i)
    sys.stdout.write("\n")
    #Print matrix
    for i in range(num_samples):
        sys.stdout.write("%i" % i)
        for j in range(num_samples):
            if i <= j:
                sys.stdout.write("\t%0.3f" % bc[i][j])
            else:
                sys.stdout.write("\tx.xxx")
        sys.stdout.write("\n")

####################################################################
if __name__ == "__main__":
    main()
