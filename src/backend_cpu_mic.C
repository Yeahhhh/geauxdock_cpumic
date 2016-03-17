
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <dock.h>
#include <size.h>
#include <toggle.h>
#include <util_print.h>
#include <record.h>

#include <yeah/cpp/timer.hpp>

#include "backend_cpu_mic_mckernel.C"


void
Dock (Complex *complex, Record *record)
{
  yeah::Timer e[16];

  // data for analysis
  std::map < int, std::vector < LigRecordSingleStep > > multi_reps_records;

  e[10].Start ();

  // sizes
  const int steps_total = complex->mcpara.steps_total;
  const int steps_per_dump = complex->mcpara.steps_per_dump;

  printf ("Start kernels on msg %02d\n", complex->signal);
  printf ("steps_per_dump = %d\n", steps_per_dump);
  printf ("steps_total = %d\n", steps_total);

  e[3].Start ();
  // calculate initial energy
  MonteCarlo_d (complex, record, 0, 1);
  e[3].Stop ();

  e[4].Start ();
  for (int s1 = 0; s1 < steps_total; s1 += steps_per_dump) {
    printf ("\t%d / %d \n", s1, steps_total);
    // fflush (stdout);
    MonteCarlo_d (complex, record, s1, steps_per_dump);
    #include <kernel_dump.C>
  }
  e[4].Stop ();



/*
	for (int s = 0; s < ligrecord[rep].next_ptr; ++s) {
		LigRecordSingleStep my_step = ligrecord[rep].step[s];
		multi_reps_records[rep].push_back(my_step);
	}
*/



  e[10].Stop ();


  Complex *ch = complex;
#include "kernel_print.C"
#include "kernel_print_timer.C"
  //PrintSummary (complex);
}

