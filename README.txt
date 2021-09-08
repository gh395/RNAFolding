to get a sequence to run
1) download RNAStrAlign database https://rna.urmc.rochester.edu/publications.html. (Citation 112)
2) run scraper.py in terminal, specifying a ct file from RNAStrAlign database
e.g. >python scraper.py -f \file\path\to\file.ct

run nussinov.py, zuker.py, and ILM.py in terminal
e.g. >python nussinov.py -s GGCGGUGUAGCUCAGCUGGCUAGAGCGUCCGGUUCAUACCCGGGAGGUCGAGGGUUCGAUCCCCCCCGCCGCUA


to visualize RNA folds (must not be a pseudoknot), visit http://rna.tbi.univie.ac.at/forna/, input sequence and dot bracket output from one of the above.


table 1 in write up was computed using compare_pairs.py (The first output of compare_pairs function is the accuracy. I didn't end up using the second output but it is the ratio of bases found in experimental that match ground truth/ bases found in ground truth)

