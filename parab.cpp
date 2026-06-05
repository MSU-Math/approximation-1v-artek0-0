#include "parab.h"

namespace parab {

static double xihalf(int j, int n, const double *x)
{
    if (j <= 0)
        return x[0] - 0.5 * (x[1] - x[0]);
    if (j >= n)
        return x[n - 1] + 0.5 * (x[n - 1] - x[n - 2]);
    return 0.5 * (x[j - 1] + x[j]);
}

static void thomas(int m, double *sub, double *diag, double *sup, double *rhs)
{
    int i;
    double w;
    for (i = 1; i < m; i++) {
        w = sub[i] / diag[i - 1];
        diag[i] -= w * sup[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    rhs[m - 1] /= diag[m - 1];
    for (i = m - 2; i >= 0; i--)
        rhs[i] = (rhs[i] - sup[i] * rhs[i + 1]) / diag[i];
}

int method_init(int n, const double *x, const double *f, double *a, const double *d2)
{
    int i, m;
    double xhm, xhc, xhp;
    double al, be, HL, C, D, HR, H;
    double *sub, *diag, *sup, *rhs;

    if (n < 2)
        return -1;

    for (i = 0; i < n; i++)
        a[i] = f[i];

    m = n + 1;
    sub  = new double[m]();
    diag = new double[m]();
    sup  = new double[m]();
    rhs  = new double[m]();

    {
        double xh0 = xihalf(0, n, x);
        double xh1 = xihalf(1, n, x);
        al   = x[0] - xh0;
        be   = xh1 - x[0];
        H    = xh1 - xh0;
        diag[0] = 1.0 / al;
        sup[0]  = 1.0 / be;
        rhs[0]  = f[0] * (1.0 / al + 1.0 / be) + d2[0] * H * 0.5;
    }

    for (i = 1; i < n; i++) {
        xhm = xihalf(i - 1, n, x);
        xhc = xihalf(i,     n, x);
        xhp = xihalf(i + 1, n, x);
        al = x[i - 1] - xhm;
        be = xhc - x[i - 1];
        HL = xhc - xhm;
        C  = x[i] - xhc;
        D  = xhp - x[i];
        HR = xhp - xhc;
        sub[i]  = 1.0 / al - 1.0 / HL;
        diag[i] = 1.0 / be + 1.0 / HL + 1.0 / C + 1.0 / HR;
        sup[i]  = 1.0 / D  - 1.0 / HR;
        rhs[i]  = (1.0 / al + 1.0 / be) * f[i - 1]
                + (1.0 / C  + 1.0 / D)  * f[i];
    }

    {
        double xhm2 = xihalf(n - 1, n, x);
        double xhn  = xihalf(n,     n, x);
        al = x[n - 1] - xhm2;
        be = xhn - x[n - 1];
        H  = xhn - xhm2;
        sub[n]  = 1.0 / al;
        diag[n] = 1.0 / be;
        rhs[n]  = f[n - 1] * (1.0 / al + 1.0 / be) + d2[1] * H * 0.5;
    }

    thomas(m, sub, diag, sup, rhs);

    for (i = 0; i < m; i++)
        a[n + i] = rhs[i];

    delete[] sub;
    delete[] diag;
    delete[] sup;
    delete[] rhs;

    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n,
                      const double *x_array, const double *a_array)
{
    int lo, hi, mid, j;
    double xhj, xhj1, al, be, H;
    double c1, c2, c3, vj, vj1, fj, dx;
    (void)seg_a;
    (void)seg_b;

    if (n < 2)
        return 0.0;

    if (xval <= 0.5 * (x_array[0] + x_array[1])) {
        j = 0;
    } else if (xval >= 0.5 * (x_array[n - 2] + x_array[n - 1])) {
        j = n - 1;
    } else {
        lo = 1;
        hi = n - 2;
        while (lo < hi) {
            mid = lo + (hi - lo + 1) / 2;
            if (0.5 * (x_array[mid - 1] + x_array[mid]) <= xval)
                lo = mid;
            else
                hi = mid - 1;
        }
        j = lo;
    }

    xhj  = (j == 0)     ? x_array[0] - 0.5 * (x_array[1] - x_array[0])
                        : 0.5 * (x_array[j - 1] + x_array[j]);
    xhj1 = (j == n - 1) ? x_array[n - 1] + 0.5 * (x_array[n - 1] - x_array[n - 2])
                        : 0.5 * (x_array[j] + x_array[j + 1]);

    fj  = a_array[j];
    vj  = a_array[n + j];
    vj1 = a_array[n + j + 1];

    al = x_array[j] - xhj;
    be = xhj1 - x_array[j];
    H  = xhj1 - xhj;

    c3 = (1.0 / H) * ((vj1 - fj) / be - (fj - vj) / al);
    c2 = (fj - vj) / al - al * c3;
    c1 = vj;

    dx = xval - xhj;
    return c1 + c2 * dx + c3 * dx * dx;
}

}
