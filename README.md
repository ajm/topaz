# TOPAZ: asymmetric suffix array neighbourhood search for massive protein databases

TOPAZ is a high-performance homology search method based on asymmetric suffix array neighbourhood search (SANS), scored seeds and optimal substitution ordering. A manuscript describing its operation is currently under review.

The TOPAZ source code is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Contact

If you have an questions about installing or running TOPAZ, please contact <a href="mailto:amedlar AT gmail DOT com">Alan Medlar</a>

## Installation

Download TOPAZ source code and compile:

    git clone https://github.com/ajm/topaz.git
    cd topaz/src
    make

## Indexing

TOPAZ has an index command to generate databases. Here we create a database containing all the sequences in database.fasta starting with the prefix "DB":

    topaz index -f database.fasta -p DB

## Searching

Once the database has finished indexing, we can search the database for homologous proteins to the sequences in queries.fasta. (If you did not just do the indexing, then you should run the following command to ensure the database is in the OS buffer cache):

    cat DB.ind DB.phr DB.pin DB.psq DB.SA > /dev/null

By default, TOPAZ searches for 100 hits per query sequence with an E-values less than 1e-3:

    topaz search -f queries.fasta -p DB -T 16 > results.txt

In tests TOPAZ was found to have the same sensitivity as DIAMOND. But we also offer a fast mode (--fast) that is almost as sensitivie:

    topaz search -f queries.fasta -p DB -T 16 --fast > results.txt

By default, TOPAZ only outputs the ids of the query and subject, the E-value and the bitscore. The --blasttab option outputs BLAST format 6:

    time topaz search -f queries.fasta -p DB -T 16 --blasttab > results.txt

Or change the number of hits (-H) or the E-value threshold (-E):

    topaz search -f queries.fasta -p DB -T 16 -H 1000 -E 1e-9 > results.txt

Additional command line options can be found in the help information:

    topaz -h

