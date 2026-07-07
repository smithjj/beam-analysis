/* msquared_mex.cpp
 * MEX drop-in replacement for Msquared.calculate() and
 * Msquared.calculate_pulse_flat().
 *
 * Usage:
 *   results = msquared_mex(field, xgrid, ygrid, wavelength)
 *   results = msquared_mex(field, xgrid, ygrid, wavelength, 'pulse_flat')     // fluence-based pulse analysis
 *   results = msquared_mex(field, xgrid, ygrid, wavelength, 'second_moment')  // simple second-moment widths only
 * (The string 'pulse' is accepted as a no-op alias for the default per-slice mode.)
 *   Note: field units are arbitrary; xgrid, ygrid, and wavelength must be in meters.
 *   
 *
 * Copyright 2026 Jesse Smith (jesse.smith@as-photonics.com) and Arlee Smith (arlee.smith@as-photonics.com)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
 * and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of 
 * the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
 * THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 *
 *Install MATLAB Support for MinGW-w64 C/C++/Fortran Compiler via
 * the add-ons explorer and compile with: mex msquared_mex.cpp
 */
// changelog: July 7, 2026: negate Rx, Ry to make convention positive roc is focusing
#define _USE_MATH_DEFINES
#include "mex.h"
#include <cmath>
#include <cstring>
#include <complex>
#include <vector>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_FFTW
#include "fftw3.h"
#include <string>
#include <sys/stat.h>
#endif

#ifdef MEX_PROFILE
#define MEX_TIMER_START(name) auto __t_##name = std::chrono::high_resolution_clock::now()
#define MEX_TIMER_STOP(name, label) { \
    auto __dt_##name = std::chrono::high_resolution_clock::now() - __t_##name; \
    double __ms_##name = std::chrono::duration<double, std::milli>(__dt_##name).count(); \
    mexPrintf("[MEX_PROFILE] " label ": %.3f ms\n", __ms_##name); }
#else
#define MEX_TIMER_START(name)
#define MEX_TIMER_STOP(name, label)
#endif

using cdouble = std::complex<double>;

/* Column-major indexing: element at row iy, col ix -> iy + ix*Nx */
#define IDX(iy, ix, Nx) ((iy) + (ix) * (Nx))

/* Call MATLAB's fft2 from MEX */
static void call_matlab_fft2(const cdouble* in, cdouble* out, mwSize Nx, mwSize Ny) {
    mxArray *fin = mxCreateNumericMatrix((mwSize)Nx, (mwSize)Ny, mxDOUBLE_CLASS, mxCOMPLEX);
    double *pr = mxGetPr(fin);
    double *pi = mxGetPi(fin);
    for (mwSize i = 0; i < Nx * Ny; i++) {
        pr[i] = in[i].real();
        pi[i] = in[i].imag();
    }
    mxArray *fout = nullptr;
    mexCallMATLAB(1, &fout, 1, &fin, "fft2");
    const double *opr = mxGetPr(fout);
    const double *opi = mxGetPi(fout);  /* NULL if fft2 returns real (e.g. zero input) */
    if (opi) {
        for (mwSize i = 0; i < Nx * Ny; i++) {
            out[i] = cdouble(opr[i], opi[i]);
        }
    } else {
        for (mwSize i = 0; i < Nx * Ny; i++) {
            out[i] = cdouble(opr[i], 0.0);
        }
    }
    mxDestroyArray(fin);
    mxDestroyArray(fout);
}

/* In-place 2D fftshift (column-major) */
static void fftshift2(cdouble* arr, mwSize Nx, mwSize Ny) {
    std::vector<cdouble> tmp(Nx * Ny);
    mwSize hx = Nx / 2, hy = Ny / 2;
    for (mwSize ix = 0; ix < Ny; ix++) {
        for (mwSize iy = 0; iy < Nx; iy++) {
            mwSize si = ((iy + hx) % Nx) + ((ix + hy) % Ny) * Nx;
            tmp[iy + ix * Nx] = arr[si];
        }
    }
    std::memcpy(arr, tmp.data(), Nx * Ny * sizeof(cdouble));
}

/* In-place fftshift for real-valued data (column-major) */
static void fftshift2_real(double* arr, mwSize Nx, mwSize Ny) {
    std::vector<double> tmp(Nx * Ny);
    mwSize hx = Nx / 2, hy = Ny / 2;
    for (mwSize ix = 0; ix < Ny; ix++) {
        for (mwSize iy = 0; iy < Nx; iy++) {
            mwSize si = ((iy + hx) % Nx) + ((ix + hy) % Ny) * Nx;
            tmp[iy + ix * Nx] = arr[si];
        }
    }
    std::memcpy(arr, tmp.data(), Nx * Ny * sizeof(double));
}

#ifdef USE_FFTW
/* Return a stable path for FFTW wisdom caching.
 * On Linux/WSL:  $HOME/.cache/fftw/msquared_wisdom
 * Falls back to /tmp if HOME is unset. */
static const char* fftw_wisdom_path() {
    static std::string path;
    if (path.empty()) {
        const char* home = getenv("HOME");
        if (home) {
            path = std::string(home) + "/.cache/fftw/msquared_wisdom";
        } else {
            path = "/tmp/msquared_fftw_wisdom";
        }
    }
    return path.c_str();
}

/* In-memory FFTW wisdom cache.  Plans are accumulated in a static
 * std::string so there is zero disk I/O.  Wisdom survives MEX reloads
 * within the same MATLAB process because we import it on every entry. */
static std::string& fftw_wisdom_string() {
    static std::string s;
    return s;
}

static void fftw_import_wisdom_string() {
    const std::string& s = fftw_wisdom_string();
    if (!s.empty()) {
        fftw_import_wisdom_from_string(s.c_str());
    }
}

static void fftw_export_wisdom_string() {
    char* str = fftw_export_wisdom_to_string();
    if (str) {
        fftw_wisdom_string() = str;   // copy into std::string
        fftw_free(str);
    }
}

/* In-place batched 2-D complex FFT using MATLAB's libmwfftw3.
 * Data is stored in MATLAB column-major order with logical dimensions
 * Nx (1st dim) x Ny (2nd dim).  We pass {Ny, Nx} to FFTW so that the
 * contiguous dimension (Nx) is the last/fastest one.
 *
 * On first call for a given geometry FFTW creates a plan (ESTIMATE is
 * fast).  The plan metadata is exported to an in-memory string so that
 * later calls — even after clear mex — can reuse it without any disk
 * access. */
static void fftw_batch_2d(cdouble* data, mwSize Nx, mwSize Ny, mwSize Nt) {
    fftw_import_wisdom_string();
    int n[2] = { static_cast<int>(Ny), static_cast<int>(Nx) };
    int np = static_cast<int>(Nx * Ny);
    int howmany = static_cast<int>(Nt);
    fftw_plan plan = fftw_plan_many_dft(2, n, howmany,
                                        reinterpret_cast<fftw_complex*>(data), nullptr, 1, np,
                                        reinterpret_cast<fftw_complex*>(data), nullptr, 1, np,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan) {
        mexErrMsgIdAndTxt("msquared_mex:fftw", "fftw_plan_many_dft failed");
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_export_wisdom_string();
}
#endif

struct SliceResult {
    double xBar, yBar, wx, wy;
    double kxBar, kyBar, wkx, wky;
    double Rx, Ry;
    double wx0sq, wy0sq;
    double M2x, M2y;
};

static SliceResult compute_slice(const cdouble* Exy, const cdouble* Ek,
                                  const double* xmat, const double* ymat,
                                  const double* kxmat, const double* kymat,
                                  mwSize Nx, mwSize Ny,
                                  double dx, double dy, double k0,
                                  cdouble* Exy_flat_out) {
    /* Inner OpenMP helps when this function is called outside the outer
     * parallel slice loop (e.g. single-slice / few-slice runs).  When
     * already nested we stay serial to avoid oversubscription. */
#ifdef _OPENMP
    bool nested = (omp_get_level() > 0);
    const mwSize omp_threshold = 65536;  /* ~256 x 256 */
#else
    bool nested = false;
    const mwSize omp_threshold = 0;
#endif
    mwSize np = Nx * Ny;

    /* Pass 1: Ixy, centroids, widths (flat loop for OpenMP) */
    double U = 0;
    double xBar = 0, yBar = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:U,xBar,yBar) schedule(static) if(!nested && np >= omp_threshold)
#endif
    for (mwSize i = 0; i < np; i++) {
        double I = Exy[i].real() * Exy[i].real() + Exy[i].imag() * Exy[i].imag();
        U += I;
        xBar += xmat[i] * I;
        yBar += ymat[i] * I;
    }

    if (U == 0) {
        /* Zero-power field (e.g. unseeded wave). Return zeros/NaN for all metrics. */
        SliceResult r;
        r.xBar = 0; r.yBar = 0; r.wx = 0; r.wy = 0;
        r.kxBar = 0; r.kyBar = 0; r.wkx = 0; r.wky = 0;
        r.Rx = 0; r.Ry = 0; r.wx0sq = 0; r.wy0sq = 0;
        r.M2x = 0; r.M2y = 0;
        if (Exy_flat_out) std::memset(Exy_flat_out, 0, Nx * Ny * sizeof(cdouble));
        return r;
    }
    xBar /= U;
    yBar /= U;

    double wxsq = 0, wysq = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:wxsq,wysq) schedule(static) if(!nested && np >= omp_threshold)
#endif
    for (mwSize i = 0; i < np; i++) {
        double I = Exy[i].real() * Exy[i].real() + Exy[i].imag() * Exy[i].imag();
        double ddx = xmat[i] - xBar;
        double ddy = ymat[i] - yBar;
        wxsq += ddx * ddx * I;
        wysq += ddy * ddy * I;
    }
    wxsq = 4.0 * wxsq / U;
    wysq = 4.0 * wysq / U;

    double Uk = 0, kxBar = 0, kyBar = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:Uk,kxBar,kyBar) schedule(static) if(!nested && np >= omp_threshold)
#endif
    for (mwSize i = 0; i < np; i++) {
        double Ik = Ek[i].real() * Ek[i].real() + Ek[i].imag() * Ek[i].imag();
        Uk += Ik;
        kxBar += kxmat[i] * Ik;
        kyBar += kymat[i] * Ik;
    }
    kxBar /= Uk;
    kyBar /= Uk;

    double wkxsq = 0, wkysq = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:wkxsq,wkysq) schedule(static) if(!nested && np >= omp_threshold)
#endif
    for (mwSize i = 0; i < np; i++) {
        double Ik = Ek[i].real() * Ek[i].real() + Ek[i].imag() * Ek[i].imag();
        double ddx = kxmat[i] - kxBar;
        double ddy = kymat[i] - kyBar;
        wkxsq += ddx * ddx * Ik;
        wkysq += ddy * ddy * Ik;
    }
    wkxsq = 4.0 * wkxsq / Uk;
    wkysq = 4.0 * wkysq / Uk;

    /* Iterative curvature removal (3 passes) */
    double Rx_vals[3], Ry_vals[3];
    std::vector<cdouble> Ework(Exy, Exy + Nx * Ny);

    for (int P = 0; P < 3; P++) {
        /* Derivative in x (row direction, 1st dim).
         * MATLAB: dE_dx = (circshift(Exy,[-1,0]) - circshift(Exy,[1,0]))/(2*dx)
         * circshift along 1st dim = shift rows.
         * dE_dx([1,end],:) = 0 zeros first and last rows.
         *
         * MATLAB: Ax = 1i/k0 * sum((xmat-xBar).*(EdE_star - conj(EdE_star)), 'all') / U
         * EdE_star - conj(EdE_star) = 2i*imag(EdE_star)
         * 1i * 2i = -2
         * So: Ax = -2/k0 * sum((xmat-xBar)*imag(EdE_star)) / U
         */
        double Ax_sum = 0;
        for (mwSize ix = 0; ix < Ny; ix++) {
            for (mwSize iy = 1; iy < Nx - 1; iy++) {
                mwSize idx = IDX(iy, ix, Nx);
                double Pr = Ework[idx].real(), Qi = Ework[idx].imag();
                double ap = Ework[IDX(iy + 1, ix, Nx)].real(), aq = Ework[IDX(iy + 1, ix, Nx)].imag();
                double bp = Ework[IDX(iy - 1, ix, Nx)].real(), bq = Ework[IDX(iy - 1, ix, Nx)].imag();
                double imag_part = (Qi * (ap - bp) - Pr * (aq - bq)) / dx;
                Ax_sum += (xmat[idx] - xBar) * imag_part;
            }
        }
        double Ax = -1.0 * Ax_sum / k0 / U;
        double Rx = wxsq / (2.0 * Ax);
        Rx_vals[P] = Rx;

        /* Derivative in y (column direction, 2nd dim).
         * MATLAB: dE_dy = (circshift(Exy,[0,-1]) - circshift(Exy,[0,1]))/(2*dy)
         * Simplified: imag_part = (Q*(r-t) - P*(s-u)) / dy
         */
        double Ay_sum = 0;
        for (mwSize iy = 0; iy < Nx; iy++) {
            for (mwSize ix = 1; ix < Ny - 1; ix++) {
                mwSize idx = IDX(iy, ix, Nx);
                double Pr = Ework[idx].real(), Qi = Ework[idx].imag();
                double ap = Ework[IDX(iy, ix + 1, Nx)].real(), aq = Ework[IDX(iy, ix + 1, Nx)].imag();
                double bp = Ework[IDX(iy, ix - 1, Nx)].real(), bq = Ework[IDX(iy, ix - 1, Nx)].imag();
                double imag_part = (Qi * (ap - bp) - Pr * (aq - bq)) / dy;
                Ay_sum += (ymat[idx] - yBar) * imag_part;
            }
        }
        double Ay = -1.0 * Ay_sum / k0 / U;
        double Ry = wysq / (2.0 * Ay);
        Ry_vals[P] = Ry;

        /* Remove curvature (skip on last pass) */
        if (P < 2) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) if(!nested && np >= omp_threshold)
#endif
            for (mwSize i = 0; i < np; i++) {
                double xc = xmat[i] - xBar;
                double yc = ymat[i] - yBar;
                double phase = -k0 * (xc * xc / (2.0 * Rx) + yc * yc / (2.0 * Ry));
                double sp = std::sin(phase);
                double cp = std::cos(phase);
                double er = Ework[i].real(), ei = Ework[i].imag();
                Ework[i] = cdouble(er * cp - ei * sp, er * sp + ei * cp);
            }
        }
    }

    /* Harmonic sum of R */
    double Rx_inv = 0, Ry_inv = 0;
    for (int P = 0; P < 3; P++) {
        Rx_inv += 1.0 / Rx_vals[P];
        Ry_inv += 1.0 / Ry_vals[P];
    }
    double Rx_c = 1.0 / Rx_inv;
    double Ry_c = 1.0 / Ry_inv;

    /* Beam waist and M-squared */
    double wx0sq = wxsq - wxsq * wxsq * k0 * k0 / (wkxsq * Rx_c * Rx_c);
    double wy0sq = wysq - wysq * wysq * k0 * k0 / (wkysq * Ry_c * Ry_c);
    if (wx0sq < 0) wx0sq = 0;
    if (wy0sq < 0) wy0sq = 0;
    double M2x = std::sqrt(wx0sq) * std::sqrt(wkxsq) / 2.0;
    double M2y = std::sqrt(wy0sq) * std::sqrt(wkysq) / 2.0;

    if (Exy_flat_out) {
        std::memcpy(Exy_flat_out, Ework.data(), Nx * Ny * sizeof(cdouble));
    }

    SliceResult r;
    r.xBar = xBar; r.yBar = yBar;
    r.wx = std::sqrt(wxsq); r.wy = std::sqrt(wysq);
    r.kxBar = kxBar; r.kyBar = kyBar;
    r.wkx = std::sqrt(wkxsq); r.wky = std::sqrt(wkysq);
    r.Rx = -Rx_c; r.Ry = -Ry_c;
    r.wx0sq = wx0sq; r.wy0sq = wy0sq;
    r.M2x = M2x; r.M2y = M2y;
    return r;
}

/* Fluence-based pulse analysis (Msquared.calculate_pulse_flat).
 * Computes M^2, RoC, waist, etc. from the time-integrated fluence.
 * All outputs are scalars (single value for the whole pulse).
 */
static void compute_pulse_flat(
    const cdouble* field, mwSize Nx, mwSize Ny, mwSize Nt,
    const double* xmat, const double* ymat,
    const double* kxmat, const double* kymat,
    double dx, double dy, double k0, double wavelength,
    double& M2x, double& M2y,
    double& Rx_out, double& Ry_out,
    double& wx0_out, double& wy0_out,
    double& xBar, double& yBar,
    double& wx_out, double& wy_out,
    double& kxBar, double& kyBar,
    double& wkx_out, double& wky_out,
    double& z0x, double& z0y,
    cdouble* flattened_out) {

    mwSize np = Nx * Ny;

    /* ---- Fluence in xy space: Fxy = sum_t |E|^2 ---- */
    MEX_TIMER_START(fxy);
    std::vector<double> Fxy(np, 0.0);
    for (mwSize k = 0; k < Nt; k++) {
        const cdouble* src = field + k * np;
        for (mwSize i = 0; i < np; i++) {
            Fxy[i] += std::norm(src[i]);
        }
    }
    MEX_TIMER_STOP(fxy, "pulse_flat Fxy accumulation");

    /* ---- Batch FFT of all time slices ----
     * Use FFTW directly when available; otherwise fall back to MATLAB's fft2.
     */
    MEX_TIMER_START(fft);
#ifdef USE_FFTW
    std::vector<cdouble> fft_field(field, field + np * Nt);
    fftw_batch_2d(fft_field.data(), Nx, Ny, Nt);
#else
    mwSize fdims[3] = {Nx, Ny, Nt};
    mxArray *field_in = mxCreateNumericArray(3, fdims, mxDOUBLE_CLASS, mxCOMPLEX);
    double *fpr = mxGetPr(field_in);
    double *fpi = mxGetPi(field_in);
    for (mwSize i = 0; i < np * Nt; i++) {
        fpr[i] = field[i].real();
        fpi[i] = field[i].imag();
    }
    mxArray *fft_out = nullptr;
    mexCallMATLAB(1, &fft_out, 1, &field_in, "fft2");
    mxDestroyArray(field_in);
    const double *fft_pr = mxGetPr(fft_out);
    const double *fft_pi = mxGetPi(fft_out);
#endif
    MEX_TIMER_STOP(fft, "pulse_flat batch FFT");

    /* ---- Fluence in k-space: Fkxky = sum_t |fftshift(fft2(E))|^2 ---- */
    MEX_TIMER_START(fkxky);
    std::vector<double> Fkxky(np, 0.0);
    for (mwSize k = 0; k < Nt; k++) {
#ifdef USE_FFTW
        cdouble *Ek = fft_field.data() + k * np;
        fftshift2(Ek, Nx, Ny);
        for (mwSize i = 0; i < np; i++) {
            Fkxky[i] += std::norm(Ek[i]);
        }
#else
        std::vector<cdouble> Ek_buf(np);
        for (mwSize i = 0; i < np; i++) {
            double re = fft_pr[k * np + i];
            double im = fft_pi ? fft_pi[k * np + i] : 0.0;
            Ek_buf[i] = cdouble(re, im);
        }
        fftshift2(Ek_buf.data(), Nx, Ny);
        for (mwSize i = 0; i < np; i++) {
            Fkxky[i] += std::norm(Ek_buf[i]);
        }
#endif
    }
#ifndef USE_FFTW
    mxDestroyArray(fft_out);
#endif
    MEX_TIMER_STOP(fkxky, "pulse_flat Fkxky accumulation");

    /* ---- Total energy (proportional) ---- */
    double U = 0, Uk = 0;
    for (mwSize i = 0; i < np; i++) {
        U  += Fxy[i];
        Uk += Fkxky[i];
    }

    if (U == 0 || Uk == 0) {
        M2x = M2y = Rx_out = Ry_out = 0;
        wx0_out = wy0_out = wx_out = wy_out = 0;
        xBar = yBar = kxBar = kyBar = 0;
        wkx_out = wky_out = z0x = z0y = 0;
        if (flattened_out) std::memset(flattened_out, 0, np * Nt * sizeof(cdouble));
        return;
    }

    /* ---- Centroids from fluence (Eq 7.112, 7.113) ---- */
    xBar = 0; yBar = 0;
    kxBar = 0; kyBar = 0;
    for (mwSize i = 0; i < np; i++) {
        xBar  += xmat[i]  * Fxy[i];
        yBar  += ymat[i]  * Fxy[i];
        kxBar += kxmat[i] * Fkxky[i];
        kyBar += kymat[i] * Fkxky[i];
    }
    xBar  /= U;
    yBar  /= U;
    kxBar /= Uk;
    kyBar /= Uk;

    /* ---- Second-moment widths from fluence (Eq 7.115, 7.116) ---- */
    double wxsq = 0, wysq = 0;
    double wkxsq = 0, wkysq = 0;
    for (mwSize i = 0; i < np; i++) {
        double dx_ = xmat[i] - xBar;
        double dy_ = ymat[i] - yBar;
        wxsq += dx_ * dx_ * Fxy[i];
        wysq += dy_ * dy_ * Fxy[i];
        double dkx = kxmat[i] - kxBar;
        double dky = kymat[i] - kyBar;
        wkxsq += dkx * dkx * Fkxky[i];
        wkysq += dky * dky * Fkxky[i];
    }
    wxsq  = 4.0 * wxsq  / U;
    wysq  = 4.0 * wysq  / U;
    wkxsq = 4.0 * wkxsq / Uk;
    wkysq = 4.0 * wkysq / Uk;

    wx_out  = std::sqrt(wxsq);
    wy_out  = std::sqrt(wysq);
    wkx_out = std::sqrt(wkxsq);
    wky_out = std::sqrt(wkysq);

    /* ---- Iterative curvature removal (3 passes on full 3D field) ---- */
    MEX_TIMER_START(curv);
    std::vector<cdouble> Ework(field, field + np * Nt);
    double Rx_vals[3], Ry_vals[3];

    for (int P = 0; P < 3; P++) {
        /* Accumulator: sum over time of 2*imag(E*conj(dE/dx)) per spatial point */
        std::vector<double> xImSum(np, 0.0);
        std::vector<double> yImSum(np, 0.0);

#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<double>> xImSum_thr(nthreads, std::vector<double>(np, 0.0));
        std::vector<std::vector<double>> yImSum_thr(nthreads, std::vector<double>(np, 0.0));

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            #pragma omp for
            for (mwSize k = 0; k < Nt; k++) {
                const cdouble* Ek = Ework.data() + k * np;

                /* x-derivative (along rows, 1st dim) */
                for (mwSize ix = 0; ix < Ny; ix++) {
                    for (mwSize iy = 1; iy < Nx - 1; iy++) {
                        mwSize idx = IDX(iy, ix, Nx);
                        double Pr = Ek[idx].real(), Qi = Ek[idx].imag();
                        double ap = Ek[IDX(iy + 1, ix, Nx)].real();
                        double aq = Ek[IDX(iy + 1, ix, Nx)].imag();
                        double bp = Ek[IDX(iy - 1, ix, Nx)].real();
                        double bq = Ek[IDX(iy - 1, ix, Nx)].imag();
                        xImSum_thr[tid][idx] += (Qi * (ap - bp) - Pr * (aq - bq)) / dx;
                    }
                }

                /* y-derivative (along columns, 2nd dim) */
                for (mwSize iy = 0; iy < Nx; iy++) {
                    for (mwSize ix = 1; ix < Ny - 1; ix++) {
                        mwSize idx = IDX(iy, ix, Nx);
                        double Pr = Ek[idx].real(), Qi = Ek[idx].imag();
                        double ap = Ek[IDX(iy, ix + 1, Nx)].real();
                        double aq = Ek[IDX(iy, ix + 1, Nx)].imag();
                        double bp = Ek[IDX(iy, ix - 1, Nx)].real();
                        double bq = Ek[IDX(iy, ix - 1, Nx)].imag();
                        yImSum_thr[tid][idx] += (Qi * (ap - bp) - Pr * (aq - bq)) / dy;
                    }
                }
            }
        }

        /* Reduce per-thread accumulators */
        for (int t = 0; t < nthreads; t++) {
            for (mwSize i = 0; i < np; i++) {
                xImSum[i] += xImSum_thr[t][i];
                yImSum[i] += yImSum_thr[t][i];
            }
        }
#else
        for (mwSize k = 0; k < Nt; k++) {
            const cdouble* Ek = Ework.data() + k * np;

            /* x-derivative (along rows, 1st dim) */
            for (mwSize ix = 0; ix < Ny; ix++) {
                for (mwSize iy = 1; iy < Nx - 1; iy++) {
                    mwSize idx = IDX(iy, ix, Nx);
                    double Pr = Ek[idx].real(), Qi = Ek[idx].imag();
                    double ap = Ek[IDX(iy + 1, ix, Nx)].real();
                    double aq = Ek[IDX(iy + 1, ix, Nx)].imag();
                    double bp = Ek[IDX(iy - 1, ix, Nx)].real();
                    double bq = Ek[IDX(iy - 1, ix, Nx)].imag();
                    xImSum[idx] += (Qi * (ap - bp) - Pr * (aq - bq)) / dx;
                }
            }

            /* y-derivative (along columns, 2nd dim) */
            for (mwSize iy = 0; iy < Nx; iy++) {
                for (mwSize ix = 1; ix < Ny - 1; ix++) {
                    mwSize idx = IDX(iy, ix, Nx);
                    double Pr = Ek[idx].real(), Qi = Ek[idx].imag();
                    double ap = Ek[IDX(iy, ix + 1, Nx)].real();
                    double aq = Ek[IDX(iy, ix + 1, Nx)].imag();
                    double bp = Ek[IDX(iy, ix - 1, Nx)].real();
                    double bq = Ek[IDX(iy, ix - 1, Nx)].imag();
                    yImSum[idx] += (Qi * (ap - bp) - Pr * (aq - bq)) / dy;
                }
            }
        }
#endif

        /* Spatial sum for Ax, Ay */
        double Ax_sum = 0, Ay_sum = 0;
        for (mwSize i = 0; i < np; i++) {
            Ax_sum += (xmat[i] - xBar) * xImSum[i];
            Ay_sum += (ymat[i] - yBar) * yImSum[i];
        }

        double Ax = -1.0 * Ax_sum / k0 / U;
        double Ay = -1.0 * Ay_sum / k0 / U;
        double Rx = wxsq / (2.0 * Ax);
        double Ry = wysq / (2.0 * Ay);
        Rx_vals[P] = Rx;
        Ry_vals[P] = Ry;

        /* Remove curvature from FULL 3D field (skip on last pass) */
        if (P < 2) {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for (mwSize k = 0; k < Nt; k++) {
                cdouble* Ek = Ework.data() + k * np;
                for (mwSize ix = 0; ix < Ny; ix++) {
                    for (mwSize iy = 0; iy < Nx; iy++) {
                        mwSize idx = IDX(iy, ix, Nx);
                        double xc = xmat[idx] - xBar;
                        double yc = ymat[idx] - yBar;
                        double phase = -k0 * (xc * xc / (2.0 * Rx) + yc * yc / (2.0 * Ry));
                        double sp = std::sin(phase);
                        double cp = std::cos(phase);
                        double er = Ek[idx].real(), ei = Ek[idx].imag();
                        Ek[idx] = cdouble(er * cp - ei * sp, er * sp + ei * cp);
                    }
                }
            }
        }
    }

    MEX_TIMER_STOP(curv, "pulse_flat curvature removal (3 passes)");

    /* Harmonic sum of R */
    double Rx_inv = 0, Ry_inv = 0;
    for (int P = 0; P < 3; P++) {
        Rx_inv += 1.0 / Rx_vals[P];
        Ry_inv += 1.0 / Ry_vals[P];
    }
    double Rx_c = 1.0 / Rx_inv;
    double Ry_c = 1.0 / Ry_inv;
    Rx_out = -Rx_c;
    Ry_out = -Ry_c;

    /* Beam waist and M-squared (Eq 7.119, 7.120) */
    double wx0sq = wxsq - wxsq * wxsq * k0 * k0 / (wkxsq * Rx_c * Rx_c);
    double wy0sq = wysq - wysq * wysq * k0 * k0 / (wkysq * Ry_c * Ry_c);
    if (wx0sq < 0) wx0sq = 0;
    if (wy0sq < 0) wy0sq = 0;
    wx0_out = std::sqrt(wx0sq);
    wy0_out = std::sqrt(wy0sq);
    M2x = std::sqrt(wx0sq) * std::sqrt(wkxsq) / 2.0;
    M2y = std::sqrt(wy0sq) * std::sqrt(wkysq) / 2.0;

    /* Rayleigh range — matches MATLAB class calculate_pulse_flat */
    z0x = k0 * k0 * wxsq / (Rx_vals[0] * wkxsq);
    z0y = k0 * k0 * wysq / (Ry_vals[0] * wkysq);

    /* Flattened field */
    if (flattened_out) {
        std::memcpy(flattened_out, Ework.data(), np * Nt * sizeof(cdouble));
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 4) {
        mexErrMsgIdAndTxt("msquared_mex:nargin",
            "Need at least 4 inputs: field, xgrid, ygrid, wavelength");
    }

    const mxArray *field_mx = prhs[0];
    const double *xgrid = mxGetPr(prhs[1]);
    const double *ygrid = mxGetPr(prhs[2]);
    double wavelength = mxGetScalar(prhs[3]);

    mwSize Nx = mxGetM(field_mx);
    const mwSize *dims = mxGetDimensions(field_mx);
    mwSize Ny = dims[1];
    mwSize numel = mxGetNumberOfElements(field_mx);
    mwSize N3 = (mxGetNumberOfDimensions(field_mx) >= 3) ? dims[2] : 1;

    const double *fpr = mxGetPr(field_mx);
    const double *fpi = mxGetPi(field_mx);
    std::vector<cdouble> field(numel);
    if (fpi) {
        for (mwSize i = 0; i < numel; i++) {
            field[i] = cdouble(fpr[i], fpi[i]);
        }
    } else {
        for (mwSize i = 0; i < numel; i++) {
            field[i] = cdouble(fpr[i], 0.0);
        }
    }

    /* Determine mode */
    bool pulse_flat_mode = false;
    bool second_moment_mode = false;
    for (int i = 4; i < nrhs; i++) {
        if (mxIsChar(prhs[i])) {
            char buf[64];
            mxGetString(prhs[i], buf, sizeof(buf));
            if (strcmp(buf, "pulse_flat") == 0 || strcmp(buf, "calculate_pulse_flat") == 0)
                pulse_flat_mode = true;
            else if (strcmp(buf, "second_moment") == 0 || strcmp(buf, "second") == 0)
                second_moment_mode = true;
            /* 'pulse' is accepted but is a no-op: default per-slice mode already
             * matches the class's calculate() method. */
        }
    }

    /* Build meshgrid. Msquared expects:
     *   xgrid varies in 1st dim (constant in 2nd)
     *   ygrid varies in 2nd dim (constant in 1st)
     * Standard MATLAB meshgrid(xvec, yvec) gives the opposite.
     * We detect orientation and swap if needed. */
    std::vector<double> xmat(Nx * Ny), ymat(Nx * Ny);
    mwSize xgM = mxGetM(prhs[1]), xgN = mxGetN(prhs[1]);
    mwSize ygM = mxGetM(prhs[2]), ygN = mxGetN(prhs[2]);

    if (xgM == Nx && xgN == Ny && ygM == Nx && ygN == Ny) {
        /* Check orientation by seeing which dimension is constant.
         * xmat should be constant in 2nd dim (col): xmat(iy,ix) == xmat(iy,ix+1)
         * ymat should be constant in 1st dim (row): ymat(iy,ix) == ymat(iy+1,ix) */
        bool x_const_in_col = true;
        bool y_const_in_row = true;
        for (mwSize ix = 0; ix < Ny - 1 && x_const_in_col; ix++) {
            if (std::abs(xgrid[IDX(0, ix, Nx)] - xgrid[IDX(0, ix + 1, Nx)]) > 1e-15)
                x_const_in_col = false;
        }
        for (mwSize iy = 0; iy < Nx - 1 && y_const_in_row; iy++) {
            if (std::abs(ygrid[IDX(iy, 0, Nx)] - ygrid[IDX(iy + 1, 0, Nx)]) > 1e-15)
                y_const_in_row = false;
        }

        if (x_const_in_col && y_const_in_row) {
            /* Correct orientation */
            std::memcpy(xmat.data(), xgrid, Nx * Ny * sizeof(double));
            std::memcpy(ymat.data(), ygrid, Nx * Ny * sizeof(double));
        } else {
            /* Swapped: swap xgrid and ygrid */
            for (mwSize i = 0; i < Nx * Ny; i++) {
                xmat[i] = ygrid[i];
                ymat[i] = xgrid[i];
            }
        }
    } else {
        /* Vectors -> build meshgrid with correct orientation.
         * xmat(iy,ix) = xvec(iy), ymat(iy,ix) = yvec(ix) */
        std::vector<double> xv(Nx), yv(Ny);
        if (xgM == Nx) {
            for (mwSize i = 0; i < Nx; i++) xv[i] = xgrid[i];
        } else {
            for (mwSize i = 0; i < Nx; i++) xv[i] = xgrid[i];
        }
        if (ygN == Ny) {
            for (mwSize i = 0; i < Ny; i++) yv[i] = ygrid[i];
        } else {
            for (mwSize i = 0; i < Ny; i++) yv[i] = ygrid[i];
        }
        for (mwSize ix = 0; ix < Ny; ix++) {
            for (mwSize iy = 0; iy < Nx; iy++) {
                xmat[IDX(iy, ix, Nx)] = xv[iy];
                ymat[IDX(iy, ix, Nx)] = yv[ix];
            }
        }
    }

    /* ---- SECOND MOMENT mode: simple spatial widths only (no FFT / M2) ---- */
    if (second_moment_mode) {
        mwSize np_sm = Nx * Ny;
        mxArray* wx_arr = mxCreateDoubleMatrix(N3, 1, mxREAL);
        mxArray* wy_arr = mxCreateDoubleMatrix(N3, 1, mxREAL);
        double* wx = mxGetPr(wx_arr);
        double* wy = mxGetPr(wy_arr);

#ifdef _OPENMP
        #pragma omp parallel for if(N3 > 1)
#endif
        for (mwSize K = 0; K < N3; K++) {
            const cdouble* src = field.data() + K * np_sm;
            double I_sum = 0.0;
            double x_sum = 0.0, y_sum = 0.0;
            double x2_sum = 0.0, y2_sum = 0.0;
            for (mwSize idx = 0; idx < np_sm; idx++) {
                double I = std::norm(src[idx]);
                I_sum += I;
                x_sum += xmat[idx] * I;
                y_sum += ymat[idx] * I;
                x2_sum += xmat[idx] * xmat[idx] * I;
                y2_sum += ymat[idx] * ymat[idx] * I;
            }
            if (I_sum > 0.0) {
                double xBar = x_sum / I_sum;
                double yBar = y_sum / I_sum;
                wx[K] = std::sqrt(4.0 * (x2_sum / I_sum - xBar * xBar));
                wy[K] = std::sqrt(4.0 * (y2_sum / I_sum - yBar * yBar));
            } else {
                wx[K] = 0.0;
                wy[K] = 0.0;
            }
        }

        const char* sm_names[] = {"wx", "wy"};
        plhs[0] = mxCreateStructMatrix(1, 1, 2, sm_names);
        mxSetField(plhs[0], 0, "wx", wx_arr);
        mxSetField(plhs[0], 0, "wy", wy_arr);
        return;
    }

    double dx = xmat[IDX(1, 0, Nx)] - xmat[IDX(0, 0, Nx)];
    double dy = ymat[IDX(0, 1, Nx)] - ymat[IDX(0, 0, Nx)];

    /* k-space grids */
    double dkx = 2.0 * M_PI / dx / (double)Nx;
    double dky = 2.0 * M_PI / dy / (double)Ny;
    std::vector<double> kxvec(Nx), kyvec(Ny);
    for (mwSize i = 0; i < Nx; i++)
        kxvec[i] = ((double)i - ((Nx % 2) ? ((double)Nx - 1.0) * 0.5 : (double)Nx * 0.5)) * dkx;
    for (mwSize i = 0; i < Ny; i++)
        kyvec[i] = ((double)i - ((Ny % 2) ? ((double)Ny - 1.0) * 0.5 : (double)Ny * 0.5)) * dky;

    std::vector<double> kxmat(Nx * Ny), kymat(Nx * Ny);
    for (mwSize ix = 0; ix < Ny; ix++) {
        for (mwSize iy = 0; iy < Nx; iy++) {
            kxmat[IDX(iy, ix, Nx)] = kxvec[iy];
            kymat[IDX(iy, ix, Nx)] = kyvec[ix];
        }
    }

    double k0 = 2.0 * M_PI / wavelength;
    mwSize np = Nx * Ny;

    /* ---- PULSE FLAT mode: fluence-based analysis ---- */
    if (pulse_flat_mode && N3 > 1) {
        double M2x_s, M2y_s, Rx_s, Ry_s, wx0_s, wy0_s;
        double xBar_s, yBar_s, wx_s, wy_s;
        double kxBar_s, kyBar_s, wkx_s, wky_s, z0x_s, z0y_s;
        std::vector<cdouble> flattened(np * N3);

        compute_pulse_flat(field.data(), Nx, Ny, N3,
            xmat.data(), ymat.data(),
            kxmat.data(), kymat.data(),
            dx, dy, k0, wavelength,
            M2x_s, M2y_s, Rx_s, Ry_s,
            wx0_s, wy0_s, xBar_s, yBar_s,
            wx_s, wy_s, kxBar_s, kyBar_s,
            wkx_s, wky_s, z0x_s, z0y_s,
            flattened.data());

        /* Build output struct with scalar fields (using pulse_ names) */
        const char *pfnames[] = {
            "pulse_M2_x", "pulse_M2_y", "pulse_Rx", "pulse_Ry",
            "pulse_wx0", "pulse_wy0", "pulse_xBar", "pulse_yBar",
            "pulse_wx", "pulse_wy", "pulse_kxBar", "pulse_kyBar",
            "pulse_wkx", "pulse_wky", "pulse_z0x", "pulse_z0y",
            "pulse_flattened_Exyt"
        };
        int nfields = 17;
        plhs[0] = mxCreateStructMatrix(1, 1, nfields, pfnames);

        auto set_scalar = [&](int idx, double val) {
            mxArray *arr = mxCreateDoubleScalar(val);
            mxSetField(plhs[0], 0, pfnames[idx], arr);
        };

        set_scalar(0,  M2x_s);   
        set_scalar(1,  M2y_s);
        set_scalar(2,  Rx_s);    
        set_scalar(3,  Ry_s);
        set_scalar(4,  wx0_s);   
        set_scalar(5,  wy0_s);
        set_scalar(6,  xBar_s);  
        set_scalar(7,  yBar_s);
        set_scalar(8,  wx_s);    
        set_scalar(9,  wy_s);
        set_scalar(10, kxBar_s); 
        set_scalar(11, kyBar_s);
        set_scalar(12, wkx_s);   
        set_scalar(13, wky_s);
        set_scalar(14, z0x_s);   
        set_scalar(15, z0y_s);

        /* Flattened field (3D complex) */
        mwSize fdims3[3] = {(mwSize)Nx, (mwSize)Ny, (mwSize)N3};
        mxArray *farr = mxCreateNumericArray(3, fdims3, mxDOUBLE_CLASS, mxCOMPLEX);
        double *fpr_out = mxGetPr(farr);
        double *fpi_out = mxGetPi(farr);
        for (mwSize i = 0; i < np * N3; i++) {
            fpr_out[i] = flattened[i].real();
            fpi_out[i] = flattened[i].imag();
        }
        mxSetField(plhs[0], 0, pfnames[16], farr);
        return;
    }

    /* ---- Standard (per-slice) mode ---- */
    /* Batch FFT all slices at once (direct FFTW, or MATLAB fft2 fallback) */
    MEX_TIMER_START(std_fft);
#ifdef USE_FFTW
    std::vector<cdouble> fft_field(field);
    fftw_batch_2d(fft_field.data(), Nx, Ny, N3);
#else
    mxArray *fft_field_mx = nullptr;
    mxArray *field_in = const_cast<mxArray*>(field_mx);
    mexCallMATLAB(1, &fft_field_mx, 1, &field_in, "fft2");
    const double *fft_pr = mxGetPr(fft_field_mx);
    const double *fft_pi = mxGetPi(fft_field_mx);
#endif
    MEX_TIMER_STOP(std_fft, "standard batch FFT");

    MEX_TIMER_START(std_sliceloop);

    /* Process slices */
    std::vector<double> M2x(N3), M2y(N3), Rx(N3), Ry(N3);
    std::vector<double> wx0(N3), wy0(N3), xBar(N3), yBar(N3);
    std::vector<double> wx(N3), wy(N3), kxBar(N3), kyBar(N3);
    std::vector<double> wkx(N3), wky(N3), z0x(N3), z0y(N3);
    std::vector<cdouble> flattened(Nx * Ny * N3);

    mwSize slice_size = Nx * Ny;
#ifdef _OPENMP
    #pragma omp parallel if(N3 > 1)
#endif
    {
        std::vector<cdouble> slice(slice_size);
        std::vector<cdouble> Ek_buf(slice_size);

#ifdef _OPENMP
        #pragma omp for
#endif
        for (mwSize K = 0; K < N3; K++) {
            const cdouble *src = field.data() + K * slice_size;
            std::memcpy(slice.data(), src, slice_size * sizeof(cdouble));

#ifdef USE_FFTW
            std::memcpy(Ek_buf.data(), fft_field.data() + K * slice_size, slice_size * sizeof(cdouble));
#else
            if (fft_pi) {
                for (mwSize i = 0; i < slice_size; i++) {
                    Ek_buf[i] = cdouble(fft_pr[K * slice_size + i], fft_pi[K * slice_size + i]);
                }
            } else {
                for (mwSize i = 0; i < slice_size; i++) {
                    Ek_buf[i] = cdouble(fft_pr[K * slice_size + i], 0.0);
                }
            }
#endif
            fftshift2(Ek_buf.data(), Nx, Ny);

            SliceResult r = compute_slice(slice.data(), Ek_buf.data(), xmat.data(), ymat.data(),
                kxmat.data(), kymat.data(), Nx, Ny, dx, dy, k0,
                flattened.data() + K * slice_size);

            M2x[K] = r.M2x; M2y[K] = r.M2y;
            Rx[K] = r.Rx;   Ry[K] = r.Ry;
            wx0[K] = std::sqrt(r.wx0sq); wy0[K] = std::sqrt(r.wy0sq);
            xBar[K] = r.xBar; yBar[K] = r.yBar;
            wx[K] = r.wx;     wy[K] = r.wy;
            kxBar[K] = r.kxBar; kyBar[K] = r.kyBar;
            wkx[K] = r.wkx;   wky[K] = r.wky;

            double wx_diff = wx[K] * wx[K] - wx0[K] * wx0[K];
            double wy_diff = wy[K] * wy[K] - wy0[K] * wy0[K];
            z0x[K] = (wx_diff > 0) ? M_PI * wx0[K] / (M2x[K] * wavelength) * std::sqrt(wx_diff) : 0;
            z0y[K] = (wy_diff > 0) ? M_PI * wy0[K] / (M2y[K] * wavelength) * std::sqrt(wy_diff) : 0;
        }
    }

    MEX_TIMER_STOP(std_sliceloop, "standard slice loop");

    MEX_TIMER_START(std_output);
    /* Build output struct */
    const char *fnames[] = {
        "M2_x", "M2_y", "Rx", "Ry", "wx0", "wy0",
        "xBar", "yBar", "wx", "wy", "kxBar", "kyBar",
        "wkx", "wky", "z0x", "z0y", "flattened_Exyz"
    };
    const char *pfnames[] = {
        "pulse_M2_x", "pulse_M2_y", "pulse_Rx", "pulse_Ry",
        "pulse_wx0", "pulse_wy0", "pulse_xBar", "pulse_yBar",
        "pulse_wx", "pulse_wy", "pulse_kxBar", "pulse_kyBar",
        "pulse_wkx", "pulse_wky", "pulse_z0x", "pulse_z0y",
        "pulse_flattened_Exyt"
    };
    const char **names = fnames;
    int nfields = 17;

    plhs[0] = mxCreateStructMatrix(1, 1, nfields, names);

    auto set_field = [&](int idx, const std::vector<double> &data) {
        mxArray *arr = mxCreateDoubleMatrix((mwSize)N3, 1, mxREAL);
        std::memcpy(mxGetPr(arr), data.data(), N3 * sizeof(double));
        mxSetField(plhs[0], 0, names[idx], arr);
    };

    set_field(0, M2x); set_field(1, M2y);
    set_field(2, Rx);  set_field(3, Ry);
    set_field(4, wx0); set_field(5, wy0);
    set_field(6, xBar); set_field(7, yBar);
    set_field(8, wx);  set_field(9, wy);
    set_field(10, kxBar); set_field(11, kyBar);
    set_field(12, wkx); set_field(13, wky);
    set_field(14, z0x); set_field(15, z0y);

    /* Flattened field (3D complex) */
    mwSize fdims3[3] = {(mwSize)Nx, (mwSize)Ny, (mwSize)N3};
    mxArray *farr = mxCreateNumericArray(3, fdims3, mxDOUBLE_CLASS, mxCOMPLEX);
    double *fpr_out = mxGetPr(farr);
    double *fpi_out = mxGetPi(farr);
    for (mwSize i = 0; i < Nx * Ny * N3; i++) {
        fpr_out[i] = flattened[i].real();
        fpi_out[i] = flattened[i].imag();
    }
    mxSetField(plhs[0], 0, names[16], farr);
    MEX_TIMER_STOP(std_output, "standard output build");
}
