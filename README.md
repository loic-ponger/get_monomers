# get_monomers
Program used to search and extract monomer sequences from tandem repeats 

- `get_monomers.py`: the main program
- `get_monomers_lib.py` the library of functions (should be in the same directory)
- `reference_DF0000029_DF0000014_DF0000015_AM237211.fst`: a fasta file including 4 alpha satellite sequences which are used as references for the search

Syntaxe:
```
./get_monomers.py -s my_file.fasta -a search_and_extract -p 1 -r ./reference_DF0000029_DF0000014_DF0000015_AM237211.fst -x my_prefix -v 9
```
