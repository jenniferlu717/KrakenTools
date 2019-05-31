# Kraken Tools
For news and updates, refer to the github page: https://github.com/jennifer.lu717/Kraken-Tools/

Kraken-Tools is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using Kraken-Tools with any of the listed programs. 

# Running Scripts:
No installation required. 
All scripts are scripts that are run on the command line as follows

## combine\_kreports.py 
Usage: python complete\_kreports.py 
    -r REPORT1.KREPORT REPORT2.KREPORT 
    -o COMBINED.KREPORT 

This script combines multiple Kraken reports into a combined report file.

Percentage is only reported for the summed read counts, not for each individual sample. 

The output file therefore contains the following tab-delimited columns:
    perc...........percentage of total reads rooted at this clade 
    tot_all .......total reads rooted at this clade (including reads at more specific clades) 
    tot_lvl........total reads at this clade  (not including reads at more specific clades)
    1_all..........reads from Sample 1 rooted at this clade 
    1_lvl..........reads from Sample 1 at this clade 
    2_all..........""
    2_lvl..........""
    etc..
    lvl_type.......Clade level type (R, D, P, C, O, F, G, S....) 
    taxid..........taxonomy ID of this clade
    name...........name of this clade 
# kreport2mpa.py  


# kreport2krona.py 

