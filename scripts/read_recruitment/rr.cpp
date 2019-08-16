#include <cassert>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <zlib.h>
#include "edlib.h"
#include "kseq/kseq.h"
KSEQ_INIT(gzFile, gzread)


char complement(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
    return ' ';
}

// https://stackoverflow.com/a/1129028
void reverse(char word[])
{
    int len=strlen(word);
    for (int i=0;i<len/2;i++)
    {
        word[i]^=word[len-i-1];
        word[len-i-1]^=word[i];
        word[i]^=word[len-i-1];
    }
}


int main(int argc, char ** argv) {
    if ((argc != 5) or (argv[1] == "-h")) {
        printf("Usage: ./rr unit.fasta reads.fasta.gz output.fasta edit_distance_threshold");
        exit(0);
    }

    auto unit_fn = argv[1];
    auto read_fn = argv[2];
    auto output_fn = argv[3];
    auto threshold = std::atoi(argv[4]);

    // read unit
    FILE * unit_file = fopen(unit_fn, "r");
    gzFile unit_fp = gzdopen(fileno(unit_file), "r");
    kseq_t *unit = kseq_init(unit_fp);
    kseq_read(unit);
    char * unit_reverse = new char[strlen(unit->seq.s)];
    std::transform(unit->seq.s,
                   unit->seq.s + unit->seq.l,
                   unit_reverse,
                   complement);
    reverse(unit_reverse);

    // open output file
    FILE * output_file = fopen(output_fn, "w");

    // open read file
    FILE * read_file = fopen(read_fn, "r");
    gzFile read_fp = gzdopen(fileno(read_file), "r");

    int n = 0;
    kseq_t *read = kseq_init(read_fp);
    while (kseq_read(read) >= 0) {
        EdlibAlignResult result = edlibAlign(unit->seq.s, unit->seq.l,
                                             read->seq.s, read->seq.l,
                                             edlibNewAlignConfig(threshold, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        EdlibAlignResult rc_result = edlibAlign(unit_reverse, unit->seq.l,
                                                read->seq.s, read->seq.l,
                                                edlibNewAlignConfig(threshold, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        // if (n % 1000 == 0) {
        //     printf("%d\n", n);
        // }
        if (result.status == EDLIB_STATUS_OK or rc_result.status == EDLIB_STATUS_OK) {
            if (result.editDistance != -1 or rc_result.editDistance != -1) {
                fprintf(output_file, ">%s\n%s\n", read->name.s, read->seq.s);
            }
        }
        edlibFreeAlignResult(result);
        ++n;
    }

    kseq_destroy(unit);
    kseq_destroy(read);
    fclose(read_file);
    gzclose(read_fp);
    gzclose(unit_fp);
}
