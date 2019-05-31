# Kraken Tools
For news and updates, refer to the github page: https://github.com/jennifer.lu717/Kraken-Tools/

Kraken-Tools is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using Kraken-Tools with any of the listed programs. 

# Running Scripts:
No installation required. 
All scripts are scripts that are run on the command line as follows

## combine\_kreports.py 
### combine\_kreports.py USAGE
Usage: python complete\_kreports.py 
    -r REPORT1.KREPORT REPORT2.KREPORT 
    -o COMBINED.KREPORT 

This script combines multiple Kraken reports into a combined report file.

### combine\_kreports.py OUTPUT 
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
## kreport2mpa.py  


## kreport2krona.py 
### kreport2krona.py USAGE
Usage: python kreport2krona.py
    -r/--report/................Kraken report file 
    -o/--output ................Output Krona text file
Optional:
    --no-intermediate-ranks.....only output standard levels [D,P,C,O,F,G,S] 
    --intermediate-ranks........[default] include non-standard levels

This program takes a Kraken report file and prints out a krona-compatible TEXT file
### kreport2krona.py EXAMPLE 
Example usage:
    kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT 
        --paired SAMPLE\_1.FASTA SAMPLE\_2.FASTA > MYSAMPLE.KRAKEN2
    python kreport2krona.py -r MYSAMPLE.KREPORT -o MYSAMPLE.krona 
    ktImportText MYSAMPLE.krona -o MYSAMPLE.krona.html
    
Krona information: see https://github.com/marbl/Krona. 

