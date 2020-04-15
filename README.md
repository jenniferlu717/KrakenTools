# Kraken Tools
For news and updates, refer to the github page: https://github.com/jenniferlu717/KrakenTools/

KrakenTools is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using KrakenTools with any of the listed programs. 

## Scripts included in KrakenTools 
1. [extract\_kraken\_reads.py](#extract\_kraken\_readspy)
2. [combine\_kreports.py](#combine\_kreportspy)
3. [kreport2krona.py](#kreport2kronapy)
4. [filter\_bracken\_out.py](#filter\_bracken\_outy)
4. [fix\_unmapped.py](#fix\_unmappedpy)

# Running Scripts:
No installation required. 
All scripts are run on the command line as described.

Users can make scripts executable by running

    chmod +x myscript.py
    ./myscript.py -h 

# extract\_kraken\_reads.py

This program extract reads classified at any user-specified taxonomy IDs. User
must specify the Kraken output file, the sequence file(s), and at least one
taxonomy ID. Additional options are specified below.

## USAGE/OPTIONS
    
    `python extract_kraken_reads.py`
    *   `-k, --kraken MYFILE.KRAKEN.............`Kraken output file
    *   `-s, -s1, -1, -U SEQUENCE.FILE..........`FASTA/FASTQ sequence file (may be gzipped)
    *   `-s2, -2 SEQUENCE2.FILE.................`FASTA/FASTQ sequence file (for paired reads, may be gzipped)
    *   `-o, --output2 OUTPUT.FASTA.............`output FASTA/Q file with extracted seqs
    *   `-t, --taxid TID TID2 etc...............`list of taxonomy IDs to extract (separated by spaces)        

    Optional:
    *   `-o2, --output2 OUTPUT.FASTA.............`second output FASTA/Q file with extracted seqs (for paired reads)
    *   `--fastq-output..........................`Instead of producing FASTA files, print FASTQ files (requires FASTQ input)
    *   `--exclude...............................`Instead of finding reads matching specified taxids, finds reads NOT matching specified taxids.
    *   `-r, --report MYFILE.KREPORT.............`Kraken report file (required if specifying --include-children or --include-parents)
    *   `--include-children......................`include reads classified at more specific levels than specified taxonomy ID levels. 
    *   `--include-parents.......................`include reads classified at all taxonomy levels between root and the specified taxonomy ID levels.
    *   `--max #.................................`maximum number of reads to save.
    *   `--append................................`if output file exists, appends reads
    *   `--noappend..............................`[default] rewrites existing output file
    
## INPUT FILES: SEQUENCE FILES

    Input sequence files must be either FASTQ or FASTA files. Input files
    can be gzipped or not. The program will automatically detect whether
    the file is gzipped and whether it is FASTQ or FASTA formatted based on
    the first character in the file (">" for FASTA, "@" for FASTQ)

## PAIRED INPUT/OUTPUT
    
    Users that ran Kraken using paired reads should input both read files into
    extract_kraken_reads.py as follows:
        
    `extract_kraken_reads.py -k myfile.kraken -s1 read1.fq -s2 reads2.fq`
    
    Given paired reads, the script requires users to provide two output file 
    names to contain extracted reads:

    `extract_kraken_reads.py -k myfile.kraken -s1 read1.fq -s2 reads2.fq -o extracted1.fq -o2 extracted2.fq`
    
    The delimiter (`--delimiter` or `-d`) option has been removed.
    
    `extract_kraken_reads.py -k myfile.kraken ... -o reads_S1.fa -o2 reads_s2.fa`

## `--exclude` flag
    
    By default, reads classified at specified taxonomy IDs will be extracted (and any taxids selected using `--include-parents`/`--include-children`. However, specifying `--exclude` will cause the reads NOT classified at any specified taxonomy IDs.

    For example: 
    1. `extract_kraken_reads.py -k myfile.kraken ... --taxid 9606 --exclude` ==> extract all reads NOT classified as Human (taxid 9606). 
    2. `extract_kraken_reads.py -k myfile.kraken ... --taxid 2 --exclude --include-children` ==> extract all reads NOT classified as Bacteria (taxid 2) or any classification in the Bacteria subtree. 
    3. `extract_kraken_reads.py -k myfile.kraken ... --taxid 9606 --exclude --include-parents` ==> extract all reads NOT classified as Human or any classification in the direct ancestry of Human (e.g. will exclude reads classified at the Primate, Chordata, or Eukaryota levels). 


## `--include-parents`/`--include-children` flags
    
    By default, only reads classified exactly at the specified taxonomy IDs
    will be extracted. Options --include-children and --include parents can be
    used to extract reads classified within the same lineage as a specified
    taxonomy ID. For example, given a Kraken report containing the following:

        [%]     [reads] [lreads][lvl]   [tid]       [name]
        100     1000    0       R       1           root
        100     1000    0       R1      131567        cellular organisms
        100     1000    50      D       2               Bacteria
        0.95    950     0       P       1224              Proteobacteria
        0.95    950     0       C       1236                Gammaproteobacteria
        0.95    950     0       O       91347                 Enterobacterales
        0.95    950     0       F       543                     Enterobacteriaceae
        0.95    950     0       G       561                       Escherichia
        0.95    950     850     S       562                         Escherichia coli
        0.05    50      50      S1      498388                        Escherichia coli C
        0.05    50      50      S1      316401                        Escherichia coli ETEC

    
    1.  `extract_kraken_reads.py  [options] -t 562` ==> 850 reads classified as _E. coli_ will be extracted
    2.  `extract_kraken_reads.py  [options] -t 562 --include-parents` ==> 900 reads classified as _E. coli_ or Bacteria will be extracted
    3.  `extract_kraken_reads.py  [options] -t 562 --include-children` ==> 950 reads classified as _E. coli_, _E. coli C_, or _E. coli ETEC_ will be extracted
    4.  `extract_kraken_reads.py  [options] -t 498388` ==> 50 reads classified as _E. coli C_ will be extracted
    5.  `extract_kraken_reads.py  [options] -t 498388 --include-parents` ==> 950 reads classified as _E. coli C_, _E. coli_, or Bacteria will be extracted
    6.  `extract_kraken_reads.py  [options] -t 1 --include-children` ==> All classified reads will be extracted 

# combine\_kreports.py 

This script combines multiple Kraken reports into a combined report file.

## USAGE/OPTIONS
    
    `python complete_kreports.py`
    *    `-r 1.KREPORT 2.KREPORT........................`Kraken-style reports to combine 
    *    `-o COMBINED.KREPORT...........................`Output file 

    Optional:
    *   `--display-headers..............................`include headers describing the samples and columns [all headers start with #]
    *   `--no-headers...................................`do not include headers in output
    *   `--sample-names.................................`give abbreviated names for each sample [default: S1, S2, ... etc]
    *   `--only-combined................................`output uses exact same columns as a single Kraken-style report file. Only total numbers for read counts and percentages will be used. Reads from individual reports will not be included.

## OUTPUT 
    Percentage is only reported for the summed read counts, not for each individual sample. 

    The output file therefore contains the following tab-delimited columns:
    *    `perc............`percentage of total reads rooted at this clade 
    *    `tot_all ........`total reads rooted at this clade (including reads at more specific clades) 
    *    `tot_lvl.........`total reads at this clade  (not including reads at more specific clades)
    *    `1_all...........`reads from Sample 1 rooted at this clade 
    *    `1_lvl...........`reads from Sample 1 at this clade 
    *    `2_all...........`""
    *    `2_lvl...........`""
    *    etc..
    *    `lvl_type........`Clade level type (R, D, P, C, O, F, G, S....) 
    *    `taxid...........`taxonomy ID of this clade
    *    `name............`name of this clade 

# kreport2krona.py 

This program takes a Kraken report file and prints out a krona-compatible TEXT file

## USAGE/OPTIONS
    
    `python kreport2krona.py`
    *    `-r/--report MYFILE.KREPORT........`Kraken report file 
    *    `-o/--output MYFILE.KRONA..........`Output Krona text file
    
    Optional:
    *    `--no-intermediate-ranks...........`only output standard levels [D,P,C,O,F,G,S] 
    *    `--intermediate-ranks..............`[default] include non-standard levels

## EXAMPLE USAGE 
    
        kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT \
            --paired SAMPLE_1.FASTA SAMPLE_2.FASTA > MYSAMPLE.KRAKEN2
        python kreport2krona.py -r MYSAMPLE.KREPORT -o MYSAMPLE.krona 
        ktImportText MYSAMPLE.krona -o MYSAMPLE.krona.html
    
    Krona information: see https://github.com/marbl/Krona. 

# filter\_bracken\_out.py

This program takes the output file of a Bracken report and filters the desired taxonomy IDs. 

## USAGE/OPTIONS
    
    `python filter_bracken_out.py`
    *   `-i/--input MYFILE.BRACKEN..........`Bracken output file
    *   `-o/--output MYFILE.BRACKEN_NEW.....`Bracken-style output file with filtered taxids
    *   `--include TID TID2.................`taxonomy IDs to include in output file [space-delimited]
    *   `--exclude TID TID2.................`taxonomy IDs to exclude in output file [space-delimited]

    User should specify either taxonomy IDs with `--include` or `--exclude`. If
    both are specified, taxonomy IDs should not be in both lists and only
    taxonomies to include will be evaluated. 
    
    When specifying the --include flag, only lines for the included taxonomy
    IDs will be extracted to the filtered output file. The percentages in the
    filtered file will be re-calculated so the total percentage in the output 
    file will sum to 100%. 

    When specifying the --exclude flag alone, all lines in the Bracken file
    will be preserved EXCEPT for the lines matching taxonomy IDs provided. 
    
## EXAMPLE USAGE

    This program can be useful for isolating a subset of species to 
    better understand the distribution of those particular species in the sample. 
    
    For example:

    * `python filter_bracken_out.py [options] --include 1764 1769 1773 1781
      39689` will allow users to get the relative percentages of 
      _Mycobacterium avium, marinum, tuberculosis, leprae, and gallinarum_ in their samples.

    In other cases, users may want to focus on the distribution of all species
    that are NOT the host species in a given sample. This program can then
    recalculate percentage distributions for species when excluding reads for
    the host.
    
    For example, given this output:
        
        name                     tax_id      tax_lvl     kraken....  added...   new.... fraction...
        Homo sapiens             9606        S           ...         ....       999000  0.999000
        Streptococcus pyogenes   1314        S           ...         ....       10      0.000001
        Streptococcus agalactiae 1311        S           ...         ....       5       0.000000
        Streptococcus pneumoniae 1313        S           ...         ....       3       0.000000
        Bordetella pertussis     520         S           ...         ....       20      0.000002
        ...
    
    Users may not be interested in the 999,000 reads that are host DNA, but
    would rather know the percentage of non-host reads for each of the non-host
    species.  Using `python filter_bracken_out.py [options] --exclude 9606`
    allows better resolution of the non-host species, allowing each of the
    fraction of reads to be recalculated out of 1,000 instead of 1,000,000
    reads in the above example. The output would then be:

        name                     tax_id      tax_lvl     kraken....  added...   new.... fraction...
        Streptococcus pyogenes   1314        S           ...         ....       10      0.01000
        Streptococcus agalactiae 1311        S           ...         ....       5       0.05000
        Streptococcus pneumoniae 1313        S           ...         ....       3       0.03000
        Bordetella pertussis     520         S           ...         ....       200     0.20000
        ...
    
# fix\_unmapped.py 
    When building a Kraken database, an "unmapped.txt" file may be generated if a taxonomy
    for a given sequence is not found. This script can search through any accession2taxid files 
    provided and the unmapped.txt file and generate a seqid2taxid.map file to be appended to
    the one already generated. 

## USAGE/OPTIONS
    
    `python fix_unmapped.py`
    *    `-i/--input unmapped.txt...........`Any file containing accession IDs to map
    *    `--accession2taxid REF_FILES.......`Any tab-delimited file with 4 columns, 
          accessions = column 1, taxonomy IDs = column 3
    *    `-o/--output OUT_FILE..............`Output tab-delimited file with 2 columns: accessions and taxids
    
    Optional:
    *    `-r/--remaining....................`file containing any unmapped accession IDs after search [default: `still_unmapped.txt`]

## EXAMPLE USAGE 
    
        rm *.k2d 
        mv seqid2taxid.map seqid2taxid_1.map
        python fix_unmapped -i unmapped.txt --accession2taxid taxonomy/*accession2taxid -o seqid2taxid_temp.map 
        cat seqid2taxid_1.map seqid2taxid_temp.map
        kraken2-build --build --db . --threads 4
    
