# Diversity Kraken Tools
For news and updates, refer to the github page: https://github.com/jenniferlu717/KrakenTools/

KrakenToolsDiversity is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using KrakenTools with any of the listed programs. 

## Links to Kraken github pages
1. [Kraken 1](https://github.com/DerrickWood/kraken)
2. [Kraken 2](https://github.com/DerrickWood/kraken2)
3. [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
4. [Bracken](https://github.com/jenniferlu717/Bracken) 

For issues with any of the above programs, 
please open a github issue on their respective github pages. 
This github repository is dedicated to only the scripts provided here. 

---------------------------------------------------------
## Scripts included in Diversity KrakenTools
1. [alpha\_diversity.py](#alpha\_diversitypy)
 
# Running Scripts:
No installation required. 
All scripts are run on the command line as described.

Users can make scripts executable by running

    chmod +x myscript.py
    ./myscript.py -h 

---------------------------------------------------------
# alpha\_diversity.py

This program calculates alpha diversity, from the Bracken abundance estimation file. 
User must specify Bracken output file, and type of 
alpha diversity to be calculated. Specific options are specified below. 

## 1. alpha\_diversity.py usage/options

`python alpha_diversity.py`
*   `-f, --filename MYFILE.BRACKEN..........`Bracken output file
*   `-a, --alpha TYPE.......................`Single letter alpha diversity type (S, BP, F)

## 2. alpha\_diversity.py input file

Input Bracken file must be in standard Bracken output file format and must be run
specifically for Abundance Estimation. Example format:
 
	name                     tax_id      tax_lvl     kraken....  added...   new.... fraction...                                                                                                         
        Homo sapiens             9606        S           ...         ....       999000  0.999000                                                                                                            
        Streptococcus pyogenes   1314        S           ...         ....       10      0.000001                                                                                                            
        Streptococcus agalactiae 1311        S           ...         ....       5       0.000000                                                                                                            
        Streptococcus pneumoniae 1313        S           ...         ....       3       0.000000                                                                                                            
        Bordetella pertussis     520         S           ...         ....       20      0.000002 

Link to Brack github for reference: [Bracken](https://github.com/jenniferlu717/Bracken)

## 4. alpha\_diversity.py alpha type input
    
By default, the program will calculate Shannon's alpha:

	python alpha_diversity.py -f myfile.bracken
 
Users can specify which type of alpha diversity from this set:

*   S.......Shannon's alpha diversity
*   BP.....Berger-Parker's alpha
*   F.......Fischer's index

To calculate berger-parker's alpha:

	python alpha_diversity.py -f myfile.bracken -a BP

---------------------------------------------------------
