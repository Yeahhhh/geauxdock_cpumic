#include <iostream>
#include <yeah/c/timing.h>
#include <yeah/cpp/timer.hpp>


int main ()
{
    yeah::Timer t0;
    yeah::TimerAuto t_auto ("Main");

    for (int i = 0; i < 10; ++i) {
        t0.Start ();
        std::cout << "aaa" << std::endl;
        t0.Stop ();
    }

    t0.Print ();
    std::cout << HostTimeNow () << std::endl;


    return 0;
}

