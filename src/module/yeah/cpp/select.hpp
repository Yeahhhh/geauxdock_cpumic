#ifndef _YEAH_CPP_SELECT_HPP_
#define _YEAH_CPP_SELECT_HPP_

namespace yeah {


    template <class T>
    T Min_2 (T a, T b)
    {
        return a < b ? a : b;
    }

    template <class T>
    T Max_2 (T a, T b)
    {
        return a > b ? a : b;
    }



// reference:
// std::max_element ()
// std::min_element ()
// std::maxmin_element ()



    template <class T>
    int Max_array_idx (T *first, T *last)
    {
        T *a = first;
        for (T *i = a; i <= last; ++i)
            if (*i > *a)
                a = i;
        return a - first;
    }


    template <class T>
    int Min_array_idx (T *first, T *last)
    {
        T *a = first;
        for (T *i = a; i <= last; ++i)
            if (*i < *a)
                a = i;
        return a - first;
    }



    template <class T>
    T Max_array (T *first, T *last)
    {
        const int idx = Max_array_idx <T> (first, last);
        return first[idx];
    }


    template <class T>
    T Min_array (T *first, T *last)
    {
        const int idx = Min_array_idx <T> (first, last);
        return first[idx];
    }





    template <class T>
    double Avg_array (T *first, T *last)
    {
        double a = 0;
        int counter = 0;
        for (T *i = first; i <= last; ++i) {
            a += *i;
            counter++;
        }
        a /= (double) counter;
        return a;
    }





//if first == last, ignore "exclude_idx"
    template <class T>
    int Min_array_idx_exclude (T *first, T *last, const int exclude_idx)
    {
        T *a = first;
        if (first != last && exclude_idx == 0)
            a = first + 1;
        for (T *i = first; i <= last; ++i)
            if (i != first + exclude_idx && *i < *a)
                a = i;
        return a - first;
    }


    template <class T>
    int Max_array_idx_exclude (T *first, T *last, const int exclude_idx)
    {
        T *a = first;
        if (first != last && exclude_idx == 0)
            a = first + 1;
        for (T *i = first; i <= last; ++i)
            if (i != first + exclude_idx && *i > *a)
                a = i;
        return a - first;
    }


    template <class T>
    T Min_array_exclude (T *first, T *last, const int exclude_idx)
    {
        const int idx = Min_array_idx_exclude <T> (first, last, exclude_idx);
        return first[idx];
    }

    template <class T>
    T Max_array_exclude (T *first, T *last, const int exclude_idx)
    {
        const int idx = Max_array_idx_exclude <T> (first, last, exclude_idx);
        return first[idx];
    }



}



#endif
