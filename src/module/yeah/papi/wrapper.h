
#ifndef _YEAH_PAPI_WRAPPER_H_
#define _YEAH_PAPI_WRAPPER_H_

#include <memory.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <papi.h>


#define PAPI_ERR(err) __GetPapiError1 (err, __FILE__, __LINE__)





#ifdef __cplusplus
extern "C" {
#endif

inline void
__GetPapiError1 (int err, const char * const file, int line)
{
    if (err < PAPI_OK) {
        printf ("%s\tFAILED\nLine # %d\n", file, line);
        if (err == PAPI_ESYS) {
            char buf[128];
            memset (buf, '\0', sizeof (buf));
            sprintf (buf, "System error:");
            perror (buf);
        }
        else if (err > 0) {
            printf ("Error calculating: \n");
        }
        else {
            printf ("Error: %s\n", PAPI_strerror (err));
        }

        exit (EXIT_FAILURE);
    }
}

#ifdef __cplusplus
}
#endif


#endif
