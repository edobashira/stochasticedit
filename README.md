Implementation of Ristad & Yianilos' estimation procedure for
memoryless stochastic transducers using OpenFst.

Got some ideas from Martin Jansche's implementiaion (ry98lsed.py)

Format
--------

Uses the *NEWS* format. One line per pair. Source and Target are seperated
by tabs and each token is seperated by space.


Usage
-------

    ./edit --num_iter=4 train.txt out.txt

Will write the best alignment to out.txt
