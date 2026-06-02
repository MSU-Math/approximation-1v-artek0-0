#include "bessel.h"

namespace bessel {

int method_init(int n, const double *x, const double *f, double *a)
{
    int i;
    (void)x;
    if (n < 2)
        return -1;
    for (i = 0; i < n; i++)
        a[i] = f[i];
    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i, j_start, num, j, k;
    double y[4], dd[4], result;
    (void)seg_a;
    (void)seg_b;

    if (n < 2)
        return 0.0;

    lo = 0;
    hi = n - 1;
    if (xval <= x_array[0]) {
        lo = 0;
    } else if (xval >= x_array[n - 1]) {
        lo = n - 2;
    } else {
        while (hi - lo > 1) {
            mid = lo + (hi - lo) / 2;
            if (x_array[mid] <= xval)
                lo = mid;
            else
                hi = mid;
        }
    }
    i = lo;

    num = (n < 4) ? n : 4;

    if (n >= 4) {
        j_start = i - 1;
        if (j_start < 0)
            j_start = 0;
        if (j_start > n - 4)
            j_start = n - 4;
    } else {
        j_start = 0;
    }

    for (j = 0; j < num; j++) {
        y[j] = x_array[j_start + j];
        dd[j] = a_array[j_start + j];
    }

    for (k = 1; k < num; k++)
        for (j = num - 1; j >= k; j--)
            dd[j] = (dd[j] - dd[j - 1]) / (y[j] - y[j - k]);

    result = dd[num - 1];
    for (j = num - 2; j >= 0; j--)
        result = result * (xval - y[j]) + dd[j];

    return result;
}

}
