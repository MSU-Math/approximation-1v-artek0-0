#include "parab.h"

namespace parab {

int method_init(int n, const double *x, const double *f, double *a, const double *d2)
{
    int i;
    double hi, hi1;

    if (n < 2)
        return -1;

    for (i = 0; i < n; i++)
        a[i] = f[i];

    a[n] = d2[0] * 0.5;

    for (i = 0; i < n - 2; i++) {
        hi = x[i + 1] - x[i];
        hi1 = x[i + 2] - x[i + 1];
        a[n + i + 1] = ((f[i + 2] - f[i + 1]) / hi1 - (f[i + 1] - f[i]) / hi - a[n + i] * hi)
                       / hi1;
    }

    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i;
    double h, bi, dx;
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

    h = x_array[i + 1] - x_array[i];
    bi = (a_array[i + 1] - a_array[i]) / h - a_array[n + i] * h;
    dx = xval - x_array[i];
    return a_array[i] + bi * dx + a_array[n + i] * dx * dx;
}

}
