
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>

#include "geauxdock.h"
#include "size.h"
#include "toggle.h"
#include "util_print.h"
#include "record.h"

#include <yeah/cpp/timer.hpp>


void
MonteCarlo_d (Complex * complex, Record * rec, const int s1, const int s2max);



void
Dock (Complex *complex, Record *record)
{
    yeah::Timer e[2];

    // data for analysis
    std::map < int, std::vector < LigRecordSingleStep > > multi_reps_records;

    e[0].Start ();

    // sizes
    const int steps_total = complex->mcpara.steps_total;
    const int steps_per_dump = complex->mcpara.steps_per_dump;

    printf ("Start kernels on msg %02d\n", complex->signal);


    MonteCarlo_d (complex, record, 0, 1); // calculate initial energy

    e[1].Start ();
    for (int s1 = 0; s1 < steps_total; s1 += steps_per_dump) {
        printf ("\t%d / %d \n", s1, steps_total);
        MonteCarlo_d (complex, record, s1, steps_per_dump);
    }
    e[1].Stop ();



    /*
       for (int s = 0; s < ligrecord[rep].next_ptr; ++s) {
       LigRecordSingleStep my_step = ligrecord[rep].step[s];
       multi_reps_records[rep].push_back(my_step);
       }
     */

    e[0].Stop ();



    // print results
    const int iter_begin = 0;
    const int iter_end = 0;
    const int arg = 2;
    //for (int r = 0; r < ph.complexsize->n_rep; ++r)
    for (int r = 0; r <= 0; ++r)
        PrintRecord (record, steps_per_dump, r, iter_begin, iter_end, arg);
    printf ("  0 0 0.6434 -0.0367 -0.2079 -0.1845 0.8521 -0.8880 0.0523 0.1744 -1.0000 0.7735 Ref Result\n"); 


    // print performance
    const float mcpersec = complex->mcpara.steps_total * complex->size.n_rep / (e[1].Span () / 1000);
    //printf ("mc kernel time\t\t\t%8.3f s\n", e[1].Span () / 1000);
    //printf ("time per MC sweep per replica\t%8.3f * 1e-6 s\n", 1e6 / mcpersec);
    //printf ("MC sweeps per second\t\t%8.3f\n", mcpersec);
    printf ("kernel wall time\t\t%8.3f ms\n", e[0].Span ());
    printf ("kernel monte carlo time\t\t%8.3f ms\n", e[1].Span ());
    printf ("speedup over prototype code\t%8.3f X\n", mcpersec / 1023.5);
    printf ("\n");
    printf ("\n");

    //PrintSummary (complex);
}



