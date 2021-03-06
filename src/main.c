#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>

#include "generic.h"
#include "options.h"
#include "substitution.h"
#include "index.h"
#include "engine.h"


void usage() {
    fprintf(stderr, "usage: %s <index|search> [-h] [-f FASTA] [-p STR]\n"
"\nGeneric:\n"
"\t-h, --help\n"
"\t-f FASTA, --file FASTA\n"
"\t-p STR, --prefix STR\n"
"\t-T INT, --threads INT (number of threads, default = %d)\n"
"\nOutput:\n"
"\t-H INT, --hits INT (default = %d)\n"
"\t-E FLOAT, --evalue FLOAT (default = %.1e)\n"
"\nAsymmetric SANS:\n"
"\t-S INT, --seeds INT (default = %d)\n"
"\t-A INT, --alignments INT (default = %d)\n"
"\t-l INT, --minlength INT (minimum length common substring, default = %d)\n"
"\t-s INT, --minscore INT (minimum raw score of common substring, default = %d)\n"
"\t-F, --fast (parameter defaults for 'fast' mode: seeds = %d, alignments = %d)"
"\nAlignment:\n"
"\t-o INT, --gapopen INT (default = %d)\n"
"\t-e INT, --gapext INT (default = %d)\n"
"\t-m STR, --matrix STR (substitution matrix, PAM*, BLOSUM*, default = %s)\n"
"\t-a FLOAT, --lambda FLOAT (E-value parameter, replace with lookup table like BLAST)\n"
"\t-k FLOAT, --kappa FLOAT (E-value parameter, replace with lookup table like BLAST)\n"
"\t-g FLOAT, --entropy FLOAT (E-value parameter, replace with lookup table like BLAST)\n"
"\t-B, --blasttab (output BLAST outfmt 6 tabular format at the expense of speed)\n"
"\t-M, --noseg (do not mask sequences with SEG during seeding)\n"
"\n", PROGRAM, DEFAULT_THREADS,
      DEFAULT_HITS, DEFAULT_EVALUE,
      DEFAULT_SEEDS, DEFAULT_ALIGNMENTS, DEFAULT_MINLENGTH, DEFAULT_MINSCORE,
      FAST_SEEDS, FAST_ALIGNMENTS,
      DEFAULT_GAPOPEN, DEFAULT_GAPEXTEND, DEFAULT_SUBSTITUTION);
}

void str2sm(char *name, options_t *opt) {
    opt->sub_name = name;
    opt->substitution = NULL;

    if(strncmp("PAM", name, 3) == 0) {
        switch(strtol(name + 3, NULL, 10)) {
            case 30 :
                opt->substitution = &PAM30;
                break;
            case 70 :
                opt->substitution = &PAM70;
                break;
            case 250 :
                opt->substitution = &PAM250;
                break;
            default :
                break;
        }
        return;
    }

    if(strncmp("BLOSUM", name, 6) == 0) {
        switch(strtol(name + 6, NULL, 10)) {
            case 45 :
                opt->substitution = &BLOSUM45;
                break;
            case 50 :
                opt->substitution = &BLOSUM50;
                break;
            case 62 :
                opt->substitution = &BLOSUM62;
                break;
            case 80 :
                opt->substitution = &BLOSUM80;
                break;
            case 90 :
                opt->substitution = &BLOSUM90;
                break;
            default :
                break;
        }
    }
}

void options_defaults(options_t *opt) {
    opt->num_threads = DEFAULT_THREADS;
    opt->num_hits = DEFAULT_HITS;
    opt->evalue = DEFAULT_EVALUE;
    opt->gap_open = DEFAULT_GAPOPEN;
    opt->gap_extend = DEFAULT_GAPEXTEND;
    opt->lambda = DEFAULT_LAMBDA;
    opt->k = DEFAULT_K;
    opt->H = DEFAULT_H;
    opt->db_size = DEFAULT_DBSIZE;
    opt->db_proteins = DEFAULT_DBPROTEINS;
    str2sm(DEFAULT_SUBSTITUTION, opt);
    opt->min_hit_length = DEFAULT_MINLENGTH;
    opt->min_hit_score = DEFAULT_MINSCORE;
    opt->number_of_seeds = DEFAULT_SEEDS;
    opt->number_of_alignments = DEFAULT_ALIGNMENTS;
    opt->superfast = 1;
    opt->seg_enabled = 1;
} 

int main(int argc, char** argv) {
    static struct option longopts[] = {
        { "help",           no_argument,        NULL, 'h' },
/* generic */
        { "file",           required_argument,  NULL, 'f' },
        { "prefix",         required_argument,  NULL, 'p' },
        { "threads",        required_argument,  NULL, 'T' },
/* output */
        { "hits",           required_argument,  NULL, 'H' },
        { "evalue",         required_argument,  NULL, 'E' },
/* alignment */
        { "gapopen",        required_argument,  NULL, 'o' },
        { "gapextend",      required_argument,  NULL, 'e' },
        { "matrix",         required_argument,  NULL, 'm' },
/* temporary */
        { "lambda",         required_argument,  NULL, 'a' },
        { "kappa",          required_argument,  NULL, 'k' },
        { "entropy",        required_argument,  NULL, 'g' },
/* seeding */
        { "minlength",      required_argument,  NULL, 'l' },
        { "minscore",       required_argument,  NULL, 's' },
        { "blasttab",       no_argument,        NULL, 'B' },
        { "alignments",     required_argument,  NULL, 'A' },
        { "seeds",          required_argument,  NULL, 'S' },
        { "noseg",          no_argument,        NULL, 'M' },
        { "fast",           no_argument,        NULL, 'F' },
        { 0, 0, 0, 0 }
    };
    int c;
    int command = -1;
    char* fname = NULL;
    char* prefix = NULL;
    char* tmp = NULL;
    options_t opt;
    search_t* search;


    if(argc < 2) {
        usage();
        exit(EXIT_FAILURE);
    }
    
    if(strcmp(argv[1], "index") == 0) {
        command = INDEX_COMMAND;
    }
    else if(strcmp(argv[1], "search") == 0) {
        command = SEARCH_COMMAND;
    }
    else {
        if((strcmp(argv[1], "-h") != 0) && (strcmp(argv[1], "--help")) != 0) {
            fprintf(stderr, "Unrecognised command: '%s' (valid commands are 'index' and 'search')\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        else {
            usage();
            exit(EXIT_SUCCESS);
        }
    }

    options_defaults(&opt);

    while((c = getopt_long(argc - 1, argv + 1, "hf:p:T:H:E:o:e:m:a:k:g:l:s:BA:S:MF", longopts, 0)) != -1) {
        switch(c) {
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'f':
                fname = optarg;
                break;
            case 'p':
                prefix = optarg;
                break;
            case 'H':
                opt.num_hits = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --hits, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                opt.gap_open = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --gapopen, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'e':
                opt.gap_extend = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --gapextend, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                str2sm(optarg, &opt);
                if(opt.substitution == NULL) {
                    fprintf(stderr, "Error: invalid argument for --substitution, (\"%s\" is not the name of a substitution matrix)\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'E':
                opt.evalue = strtod(optarg, &tmp);
                if(tmp == optarg || *tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --evalue, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'a':
                opt.lambda = strtod(optarg, &tmp);
                if(tmp == optarg || *tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --lambda, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'k':
                opt.k = strtod(optarg, &tmp);
                if(tmp == optarg || *tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --kappa, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'g':
                opt.H = strtod(optarg, &tmp);
                if(tmp == optarg || *tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --entropy, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'l':
                opt.min_hit_length = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --minlength, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 's':
                opt.min_hit_score = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --minscore, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'T':
                opt.num_threads = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --threads, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'B':
                opt.superfast = 0;
                break;
            case 'A':
                opt.number_of_alignments = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --alignments, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'S':
                opt.number_of_seeds = strtol(optarg, &tmp, 10);
                if(*tmp != '\0') {
                    fprintf(stderr, "Error: invalid argument for --seeds, %s (%s)\n", strerror(errno), optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'M':
                opt.seg_enabled = 0;
                fprintf(stderr, "WARNING: SEG disabled\n");
                break;
            case 'F' :
                opt.number_of_seeds = FAST_SEEDS;
                opt.number_of_alignments = FAST_ALIGNMENTS;
                fprintf(stderr, "fast mode enabled\n");
                break;
            case '?':
            default :
                exit(EXIT_FAILURE);
        }
    }

    argc -= optind;
    argv += optind;

    if(!fname) {
        fprintf(stderr, "ERROR: you must specify a FASTA file as input\n");
        exit(EXIT_FAILURE);
    }

    if(!prefix) {
        fprintf(stderr, "ERROR: you must specify a database prefix\n");
        exit(EXIT_FAILURE);
    }

    if(opt.min_hit_length < 3) {
        fprintf(stderr, "ERROR: minimum seed length must be at least 3\n");
        exit(EXIT_FAILURE);
    }

    if(opt.number_of_alignments < 0) {
        fprintf(stderr, "ERROR: number of alignments must be positive\n");
        exit(EXIT_FAILURE);
    }

    if(opt.number_of_seeds < 0) {
        fprintf(stderr, "ERROR: number of seeds must be positive\n");
        exit(EXIT_FAILURE);
    }

#ifdef CHATTY
    fprintf(stderr, "\nWARNING: this code uses a different internal alphabet, \n         internal2aa() in lexicographical.h needs to be called to output sequences\n\n");
#endif

    switch(command) {
        case INDEX_COMMAND :
            serif_index(fname, prefix);
            break;
        case SEARCH_COMMAND :
            search = search_alloc(fname, prefix, &opt);
            search_pthread(search);
            search_free(search);
            break;
        default :
            break;
    }


    return EXIT_SUCCESS;
}

