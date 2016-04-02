./mkliglist.py ../data/astex/ligands/ > ../data/astex/ligands/list.txt
head -n 2 ../data/astex/ligands/list.txt > ../data/astex/ligands/ligs.txt

./mkprtlist.py ../data/astex/proteins/ > ../data/astex/proteins/list.txt
head -n 1 ../data/astex/proteins/list.txt > ../data/astex/proteins/prts.txt

