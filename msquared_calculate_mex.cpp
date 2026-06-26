/*
 * msquared_calculate_mex.cpp
 *
 * MEX replacement for the Msquared_class MATLAB class.  Implements the three
 * public methods:
 *   'calculate'            -> results struct for M2 vs z/t
 *   'calculate_pulse_flat' -> results struct for fluence-based pulse-averaged M2
 *   'calculate_second_moment' -> simple x-y second-moment widths
 *
 * The mathematics is duplicated from Msquared_class.m (Arlee Smith,
 * Crystal Nonlinear Optics, chapter 7.7).
 *
 * FFTs are delegated back to MATLAB via mexCallMATLAB, so the MEX file only
 * needs a C++ compiler and does not require linking external FFTW libraries.
 *
 * Compile from MATLAB with:
 *   mex -R2018a -O msquared_calculate_mex.cpp -output msquared_calculate_mex
 */

#include <mex.h>
#include <complex>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

inline double sq(double x) { return x*x; }

/* Return a pointer to complex interleaved data for a double array.
 * If the input is already complex, no copy is made.  If it is real, a
 * temporary complex copy is created with zero imaginary part.  The caller
 * must call mxDestroyArray on tmpHolder if it is non-NULL. */
mxComplexDouble* get_complex_data_robust(const mxArray* arr, mxArray*& tmpHolder) {
    if (mxIsComplex(arr)) {
        tmpHolder = NULL;
        return mxGetComplexDoubles(arr);
    }
    tmpHolder = mxDuplicateArray(arr);
    mxMakeArrayComplex(tmpHolder);
    return mxGetComplexDoubles(tmpHolder);
}

/* Helper: fft2 followed by fftshift along dims 1 and 2.
 * The output has the same dimensions as the input. */
void fft2_shifted(mxArray* in, mxArray** out) {
    mxArray* fft_result = NULL;
    mexCallMATLAB(1, &fft_result, 1, &in, "fft2");

    mxArray* dim1 = mxCreateDoubleScalar(1.0);
    mxArray* dim2 = mxCreateDoubleScalar(2.0);

    mxArray* args1[2] = { fft_result, dim1 };
    mxArray* shifted1 = NULL;
    mexCallMATLAB(1, &shifted1, 2, args1, "fftshift");

    mxArray* args2[2] = { shifted1, dim2 };
    mexCallMATLAB(1, out, 2, args2, "fftshift");

    mxDestroyArray(fft_result);
    mxDestroyArray(shifted1);
    mxDestroyArray(dim1);
    mxDestroyArray(dim2);
}

/* Extract monotonic x and y vectors from the supplied xgrid/ygrid arrays.
 * Accepts Nx x Ny meshgrid matrices or Nx/Ny vectors. */
void extract_grid_vectors(const mxArray* mx_xgrid, const mxArray* mx_ygrid,
                          int Nx, int Ny,
                          std::vector<double>& xvec,
                          std::vector<double>& yvec) {
    if (!mxIsDouble(mx_xgrid) || !mxIsDouble(mx_ygrid)) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:invalidGrid",
            "xgrid and ygrid must be double.");
    }

    size_t sx_x = mxGetM(mx_xgrid);
    size_t sy_x = mxGetN(mx_xgrid);
    size_t sx_y = mxGetM(mx_ygrid);
    size_t sy_y = mxGetN(mx_ygrid);

    double* xg = mxGetPr(mx_xgrid);
    double* yg = mxGetPr(mx_ygrid);

    xvec.assign(Nx, 0.0);
    yvec.assign(Ny, 0.0);

    if (mxGetNumberOfElements(mx_xgrid) == (size_t)Nx * Ny && sx_x == (size_t)Nx && sy_x == (size_t)Ny) {
        for (int i = 0; i < Nx; ++i) xvec[i] = xg[i]; /* first column */
    } else if (mxGetNumberOfElements(mx_xgrid) == (size_t)Nx) {
        for (int i = 0; i < Nx; ++i) xvec[i] = xg[i];
    } else {
        mexErrMsgIdAndTxt("msquared_calculate_mex:invalidXGrid",
            "xgrid size does not match field dimensions.");
    }

    if (mxGetNumberOfElements(mx_ygrid) == (size_t)Nx * Ny && sx_y == (size_t)Nx && sy_y == (size_t)Ny) {
        for (int j = 0; j < Ny; ++j) yvec[j] = yg[j * Nx]; /* first row */
    } else if (mxGetNumberOfElements(mx_ygrid) == (size_t)Ny) {
        for (int j = 0; j < Ny; ++j) yvec[j] = yg[j];
    } else {
        mexErrMsgIdAndTxt("msquared_calculate_mex:invalidYGrid",
            "ygrid size does not match field dimensions.");
    }

    if (Nx > 1 && xvec[1] < xvec[0]) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:nonMonotonicX",
            "xgrid must be monotonically increasing.");
    }
    if (Ny > 1 && yvec[1] < yvec[0]) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:nonMonotonicY",
            "ygrid must be monotonically increasing.");
    }
}

/* Central-difference derivative along x (rows).  First and last rows set to 0. */
void derivative_x(const std::complex<double>* E, int Nx, int Ny, double dx,
                  std::complex<double>* dEdx) {
    for (int j = 0; j < Ny; ++j) {
        int idx0 = j * Nx;
        dEdx[idx0] = 0.0;
        dEdx[idx0 + Nx - 1] = 0.0;
        for (int i = 1; i < Nx - 1; ++i) {
            dEdx[idx0 + i] = (E[idx0 + i + 1] - E[idx0 + i - 1]) / (2.0 * dx);
        }
    }
}

/* Central-difference derivative along y (columns).  First and last columns set to 0. */
void derivative_y(const std::complex<double>* E, int Nx, int Ny, double dy,
                  std::complex<double>* dEdy) {
    for (int i = 0; i < Nx; ++i) {
        dEdy[i] = 0.0;
        dEdy[i + (Ny - 1) * Nx] = 0.0;
    }
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 0; i < Nx; ++i) {
            dEdy[i + j * Nx] = (E[i + (j + 1) * Nx] - E[i + (j - 1) * Nx]) / (2.0 * dy);
        }
    }
}

/* Create an output struct field containing a real column vector of length N. */
mxArray* create_real_vector(size_t N) {
    mwSize dims[2] = { (mwSize)N, 1 };
    return mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
}

/* Same but complex 3-D array Nx x Ny x N3. */
mxArray* create_complex_3d(int Nx, int Ny, int N3) {
    mwSize dims[3] = { (mwSize)Nx, (mwSize)Ny, (mwSize)N3 };
    return mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
}

/* Copy std::complex<double> buffer into an mxCOMPLEX array of the same total size. */
void copy_complex_to_mxarray(const std::complex<double>* src, mxArray* arr) {
    mxComplexDouble* dst = mxGetComplexDoubles(arr);
    size_t N = mxGetNumberOfElements(arr);
    for (size_t i = 0; i < N; ++i) {
        dst[i].real = src[i].real();
        dst[i].imag = src[i].imag();
    }
}

/* Compute k-space second moments from already fftshift-ed (dims 1 and 2) data.
 * kdata points to one slice offset. */
void kspace_moments_shifted(const mxComplexDouble* kdata,
                            int Nx, int Ny,
                            const std::vector<double>& kxvec,
                            const std::vector<double>& kyvec,
                            double& kxBar, double& kyBar,
                            double& wkx2, double& wky2) {
    double ksum = 0.0;
    double kx_sum = 0.0, ky_sum = 0.0;
    double kx2_sum = 0.0, ky2_sum = 0.0;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            size_t idx = (size_t)i + (size_t)j * Nx;
            double re = kdata[idx].real;
            double im = kdata[idx].imag;
            double I = re*re + im*im;
            ksum += I;
            kx_sum += kxvec[i] * I;
            ky_sum += kyvec[j] * I;
            kx2_sum += sq(kxvec[i]) * I;
            ky2_sum += sq(kyvec[j]) * I;
        }
    }

    if (ksum > 0.0) {
        kxBar = kx_sum / ksum;
        kyBar = ky_sum / ksum;
        wkx2 = 4.0 * (kx2_sum / ksum - sq(kxBar));
        wky2 = 4.0 * (ky2_sum / ksum - sq(kyBar));
        if (wkx2 < 0.0) wkx2 = 0.0;
        if (wky2 < 0.0) wky2 = 0.0;
    } else {
        kxBar = kyBar = 0.0;
        wkx2 = wky2 = 0.0;
    }
}

/* Accumulate spatial intensity moments of E. */
void spatial_moments(const std::complex<double>* E,
                     int Nx, int Ny,
                     const std::vector<double>& xmat,
                     const std::vector<double>& ymat,
                     double& xBar, double& yBar,
                     double& wx2, double& wy2) {
    double I_sum = 0.0;
    double x_sum = 0.0, y_sum = 0.0;
    double x2_sum = 0.0, y2_sum = 0.0;
    size_t Ntot = (size_t)Nx * Ny;

    for (size_t idx = 0; idx < Ntot; ++idx) {
        double I = std::norm(E[idx]);
        I_sum += I;
        x_sum += xmat[idx] * I;
        y_sum += ymat[idx] * I;
        x2_sum += sq(xmat[idx]) * I;
        y2_sum += sq(ymat[idx]) * I;
    }

    if (I_sum > 0.0) {
        xBar = x_sum / I_sum;
        yBar = y_sum / I_sum;
        wx2 = 4.0 * (x2_sum / I_sum - sq(xBar));
        wy2 = 4.0 * (y2_sum / I_sum - sq(yBar));
        if (wx2 < 0.0) wx2 = 0.0;
        if (wy2 < 0.0) wy2 = 0.0;
    } else {
        xBar = yBar = 0.0;
        wx2 = wy2 = 0.0;
    }
}

/* Compute Ax and Rx for the current complex field slice E. */
void compute_curvature_x(const std::complex<double>* E,
                         const std::complex<double>* dEdx,
                         int Nx, int Ny,
                         const std::vector<double>& xmat,
                         double xBar, double I_sum,
                         double wx2, double k0,
                         double& Ax, double& Rx) {
    std::complex<double> accum(0.0, 0.0);
    size_t Ntot = (size_t)Nx * Ny;
    for (size_t idx = 0; idx < Ntot; ++idx) {
        std::complex<double> z = E[idx] * std::conj(dEdx[idx]);
        std::complex<double> diff = z - std::conj(z);
        accum += (xmat[idx] - xBar) * diff;
    }
    std::complex<double> val = std::complex<double>(0.0, 1.0) / k0 * accum / I_sum;
    Ax = val.real();
    if (std::abs(Ax) < 1e-40) Ax = 1e-40;
    Rx = wx2 / (2.0 * Ax);
}

void compute_curvature_y(const std::complex<double>* E,
                         const std::complex<double>* dEdy,
                         int Nx, int Ny,
                         const std::vector<double>& ymat,
                         double yBar, double I_sum,
                         double wy2, double k0,
                         double& Ay, double& Ry) {
    std::complex<double> accum(0.0, 0.0);
    size_t Ntot = (size_t)Nx * Ny;
    for (size_t idx = 0; idx < Ntot; ++idx) {
        std::complex<double> z = E[idx] * std::conj(dEdy[idx]);
        std::complex<double> diff = z - std::conj(z);
        accum += (ymat[idx] - yBar) * diff;
    }
    std::complex<double> val = std::complex<double>(0.0, 1.0) / k0 * accum / I_sum;
    Ay = val.real();
    if (std::abs(Ay) < 1e-40) Ay = 1e-40;
    Ry = wy2 / (2.0 * Ay);
}

/* Remove spherical curvature from a 3D array E (Nx x Ny x Nz) in place. */
void remove_curvature(std::complex<double>* E,
                      int Nx, int Ny, int Nz,
                      const std::vector<double>& xmat,
                      const std::vector<double>& ymat,
                      double xBar, double yBar,
                      double Rx, double Ry, double k0) {
    size_t plane = (size_t)Nx * Ny;
    for (int z = 0; z < Nz; ++z) {
        size_t off = (size_t)z * plane;
        for (size_t idx = 0; idx < plane; ++idx) {
            double phi = -k0 * ( sq(xmat[idx] - xBar) / (2.0 * Rx) +
                                  sq(ymat[idx] - yBar) / (2.0 * Ry) );
            std::complex<double> phase = std::polar(1.0, phi);
            E[off + idx] *= phase;
        }
    }
}

/* Fill xmat/ymat meshgrid vectors from extracted vectors. */
void fill_meshgrid(int Nx, int Ny,
                   const std::vector<double>& xvec,
                   const std::vector<double>& yvec,
                   std::vector<double>& xmat,
                   std::vector<double>& ymat) {
    xmat.resize((size_t)Nx * Ny);
    ymat.resize((size_t)Nx * Ny);
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            size_t idx = (size_t)i + (size_t)j * Nx;
            xmat[idx] = xvec[i];
            ymat[idx] = yvec[j];
        }
    }
}

void build_kvec(int N, double d, std::vector<double>& kvec) {
    kvec.resize(N);
    if (N % 2 == 0) {
        double half = -N / 2.0 * d;
        for (int i = 0; i < N; ++i) kvec[i] = half + i * d;
    } else {
        double half = -(N - 1) / 2.0 * d;
        for (int i = 0; i < N; ++i) kvec[i] = half + i * d;
    }
}

/* ------------------------------------------------------------------------- */
/* calculate (M2 vs z or t)                                                  */
/* ------------------------------------------------------------------------- */
void run_calculate(int nlhs, mxArray* plhs[],
                   const mxArray* prhs[]) {
    mxArray* field_tmp = NULL;
    mxComplexDouble* mx_field = get_complex_data_robust(prhs[1], field_tmp);
    const mwSize* dims = mxGetDimensions(prhs[1]);
    size_t ndims = mxGetNumberOfDimensions(prhs[1]);
    int Nx = (int)dims[0];
    int Ny = (int)dims[1];
    int N3 = (ndims > 2) ? (int)dims[2] : 1;

    std::vector<double> xvec, yvec;
    extract_grid_vectors(prhs[2], prhs[3], Nx, Ny, xvec, yvec);

    double wavelength = mxGetScalar(prhs[4]);
    double power_threshold = mxGetScalar(prhs[5]);
    double k0 = 2.0 * M_PI / wavelength;

    double dx = (Nx > 1) ? (xvec[1] - xvec[0]) : 1.0;
    double dy = (Ny > 1) ? (yvec[1] - yvec[0]) : 1.0;

    std::vector<double> xmat, ymat;
    fill_meshgrid(Nx, Ny, xvec, yvec, xmat, ymat);

    std::vector<double> kxvec, kyvec;
    build_kvec(Nx, 2.0 * M_PI / dx / Nx, kxvec);
    build_kvec(Ny, 2.0 * M_PI / dy / Ny, kyvec);

    /* Determine slice range if power_threshold is active. */
    std::vector<double> slice_power(N3, 0.0);
    double max_power = 0.0;
    if (power_threshold > 0.0) {
        size_t slice_stride = (size_t)Nx * Ny;
        for (int k = 0; k < N3; ++k) {
            double p = 0.0;
            for (size_t i = 0; i < slice_stride; ++i) {
                p += mx_field[k * slice_stride + i].real * mx_field[k * slice_stride + i].real +
                     mx_field[k * slice_stride + i].imag * mx_field[k * slice_stride + i].imag;
            }
            slice_power[k] = p;
            if (p > max_power) max_power = p;
        }
    }

    int k_start = 0, k_end = N3 - 1;
    bool truncated = false;
    if (power_threshold > 0.0) {
        double thr = power_threshold * max_power;
        bool found_start = false;
        for (int k = 0; k < N3; ++k) {
            if (slice_power[k] > thr) {
                if (!found_start) { k_start = k; found_start = true; }
                k_end = k;
            }
        }
        if (!found_start) {
            k_start = 0; k_end = N3 - 1;
        } else {
            truncated = (k_start != 0) || (k_end != N3 - 1);
        }
    }

    /* Compute the 2D FFT of the whole 3D field, shifted along x and y.
     * For calculate() we only need the k-space representation once per slice,
     * so one fft2 call on the whole field suffices. */
    mxArray* mx_k = NULL;
    fft2_shifted(const_cast<mxArray*>(prhs[1]), &mx_k);
    mxComplexDouble* kdata_full = mxGetComplexDoubles(mx_k);

    /* Create output struct. */
    const char* field_names[] = {
        "z0x2", "z0y2", "z0x", "z0y", "M2_x", "M2_y", "Rx", "Ry",
        "wx0", "wy0", "xBar", "yBar", "wx", "wy", "kxBar", "kyBar",
        "wkx", "wky", "flattened_Exyz"
    };
    int nfields = sizeof(field_names) / sizeof(field_names[0]);
    mxArray* S = mxCreateStructMatrix(1, 1, nfields, field_names);
    plhs[0] = S;

    mxArray* mx_M2x = create_real_vector(N3);
    mxArray* mx_M2y = create_real_vector(N3);
    mxArray* mx_Rx  = create_real_vector(N3);
    mxArray* mx_Ry  = create_real_vector(N3);
    mxArray* mx_wx0 = create_real_vector(N3);
    mxArray* mx_wy0 = create_real_vector(N3);
    mxArray* mx_xBar = create_real_vector(N3);
    mxArray* mx_yBar = create_real_vector(N3);
    mxArray* mx_wx  = create_real_vector(N3);
    mxArray* mx_wy  = create_real_vector(N3);
    mxArray* mx_kxBar = create_real_vector(N3);
    mxArray* mx_kyBar = create_real_vector(N3);
    mxArray* mx_wkx = create_real_vector(N3);
    mxArray* mx_wky = create_real_vector(N3);
    mxArray* mx_z0x2 = create_real_vector(N3);
    mxArray* mx_z0y2 = create_real_vector(N3);
    mxArray* mx_z0x = create_real_vector(N3);
    mxArray* mx_z0y = create_real_vector(N3);
    mxArray* mx_flat = create_complex_3d(Nx, Ny, N3);

    double* pM2x = mxGetPr(mx_M2x); double* pM2y = mxGetPr(mx_M2y);
    double* pRx = mxGetPr(mx_Rx);   double* pRy = mxGetPr(mx_Ry);
    double* pwx0 = mxGetPr(mx_wx0); double* pwy0 = mxGetPr(mx_wy0);
    double* pxBar = mxGetPr(mx_xBar); double* pyBar = mxGetPr(mx_yBar);
    double* pwx = mxGetPr(mx_wx);   double* pwy = mxGetPr(mx_wy);
    double* pkxBar = mxGetPr(mx_kxBar); double* pkyBar = mxGetPr(mx_kyBar);
    double* pwkx = mxGetPr(mx_wkx); double* pwky = mxGetPr(mx_wky);
    double* pz0x2 = mxGetPr(mx_z0x2); double* pz0y2 = mxGetPr(mx_z0y2);
    double* pz0x = mxGetPr(mx_z0x); double* pz0y = mxGetPr(mx_z0y);
    mxComplexDouble* pflat = mxGetComplexDoubles(mx_flat);

    /* The time-slice loop is embarrassingly parallel; each slice writes to
     * its own distinct region of the output arrays.  Per-thread scratch
     * vectors are grown once outside the loop. */
    std::vector<std::complex<double>> E_cur, dEdx, dEdy;
#ifdef _OPENMP
    #pragma omp parallel private(E_cur, dEdx, dEdy)
#endif
    {
        E_cur.resize((size_t)Nx * Ny);
        dEdx.resize((size_t)Nx * Ny);
        dEdy.resize((size_t)Nx * Ny);
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int k = k_start; k <= k_end; ++k) {
            size_t off = (size_t)k * Nx * Ny;
            for (size_t i = 0; i < E_cur.size(); ++i) {
                pflat[off + i].real = mx_field[off + i].real;
                pflat[off + i].imag = mx_field[off + i].imag;
                E_cur[i] = std::complex<double>(mx_field[off + i].real, mx_field[off + i].imag);
            }

            double xBar, yBar, wx2, wy2;
            spatial_moments(E_cur.data(), Nx, Ny, xmat, ymat, xBar, yBar, wx2, wy2);

            double kxBar, kyBar, wkx2, wky2;
            kspace_moments_shifted(&kdata_full[off], Nx, Ny, kxvec, kyvec,
                                   kxBar, kyBar, wkx2, wky2);

            double Rxmat[3], Rymat[3];
            for (int P = 0; P < 3; ++P) {
                derivative_x(E_cur.data(), Nx, Ny, dx, dEdx.data());
                derivative_y(E_cur.data(), Nx, Ny, dy, dEdy.data());

                double I_sum = 0.0;
                for (size_t i = 0; i < E_cur.size(); ++i) I_sum += std::norm(E_cur[i]);
                if (I_sum == 0.0) I_sum = 1.0;

                double Ax, Rx, Ay, Ry;
                compute_curvature_x(E_cur.data(), dEdx.data(), Nx, Ny, xmat, xBar, I_sum, wx2, k0, Ax, Rx);
                compute_curvature_y(E_cur.data(), dEdy.data(), Nx, Ny, ymat, yBar, I_sum, wy2, k0, Ay, Ry);
                Rxmat[P] = Rx;
                Rymat[P] = Ry;

                if (P < 2) {
                    remove_curvature(E_cur.data(), Nx, Ny, 1, xmat, ymat, xBar, yBar, Rx, Ry, k0);
                }
            }

            double Rx_corr = 1.0 / (1.0 / Rxmat[0] + 1.0 / Rxmat[1] + 1.0 / Rxmat[2]);
            double Ry_corr = 1.0 / (1.0 / Rymat[0] + 1.0 / Rymat[1] + 1.0 / Rymat[2]);

            double wx0sq = wx2 - sq(wx2) * k0 * k0 / (wkx2 * sq(Rx_corr));
            double wy0sq = wy2 - sq(wy2) * k0 * k0 / (wky2 * sq(Ry_corr));
            if (wx0sq < 0.0) wx0sq = 0.0;
            if (wy0sq < 0.0) wy0sq = 0.0;

            double wx_val = std::sqrt(wx2);
            double wy_val = std::sqrt(wy2);
            double wkx_val = std::sqrt(wkx2);
            double wky_val = std::sqrt(wky2);
            double wx0_val = std::sqrt(wx0sq);
            double wy0_val = std::sqrt(wy0sq);

            double M2x = wx0_val * wkx_val / 2.0;
            double M2y = wy0_val * wky_val / 2.0;
            if (M2x <= 0.0) M2x = 1e-40;
            if (M2y <= 0.0) M2y = 1e-40;

            double z0x2_val = -k0 * k0 * wx2 / (Rx_corr * wkx2);
            double z0y2_val = -k0 * k0 * wy2 / (Ry_corr * wky2);
            double z0x_val  = M_PI * wx0_val / (M2x * wavelength) * std::sqrt(std::max(0.0, wx2 - wx0sq));
            double z0y_val  = M_PI * wy0_val / (M2y * wavelength) * std::sqrt(std::max(0.0, wy2 - wy0sq));

            pM2x[k] = M2x; pM2y[k] = M2y;
            pRx[k] = Rx_corr; pRy[k] = Ry_corr;
            pwx0[k] = wx0_val; pwy0[k] = wy0_val;
            pxBar[k] = xBar; pyBar[k] = yBar;
            pwx[k] = wx_val; pwy[k] = wy_val;
            pkxBar[k] = kxBar; pkyBar[k] = kyBar;
            pwkx[k] = wkx_val; pwky[k] = wky_val;
            pz0x2[k] = z0x2_val; pz0y2[k] = z0y2_val;
            pz0x[k] = z0x_val; pz0y[k] = z0y_val;

            size_t flat_off = (size_t)k * Nx * Ny;
            for (size_t i = 0; i < E_cur.size(); ++i) {
                pflat[flat_off + i].real = E_cur[i].real();
                pflat[flat_off + i].imag = E_cur[i].imag();
            }
        }
    }

    mxDestroyArray(mx_k);
    if (field_tmp) mxDestroyArray(field_tmp);

    if (truncated) {
        for (int k = 0; k < k_start; ++k) {
            pM2x[k] = pM2x[k_start]; pM2y[k] = pM2y[k_start];
            pRx[k] = pRx[k_start]; pRy[k] = pRy[k_start];
            pwx0[k] = pwx0[k_start]; pwy0[k] = pwy0[k_start];
            pxBar[k] = pxBar[k_start]; pyBar[k] = pyBar[k_start];
            pwx[k] = pwx[k_start]; pwy[k] = pwy[k_start];
            pkxBar[k] = pkxBar[k_start]; pkyBar[k] = pkyBar[k_start];
            pwkx[k] = pwkx[k_start]; pwky[k] = pwky[k_start];
            pz0x2[k] = pz0x2[k_start]; pz0y2[k] = pz0y2[k_start];
            pz0x[k] = pz0x[k_start]; pz0y[k] = pz0y[k_start];
            size_t dst = (size_t)k * Nx * Ny;
            size_t src = (size_t)k_start * Nx * Ny;
            for (size_t i = 0; i < (size_t)Nx * Ny; ++i) {
                pflat[dst + i] = pflat[src + i];
            }
        }
        for (int k = k_end + 1; k < N3; ++k) {
            pM2x[k] = pM2x[k_end]; pM2y[k] = pM2y[k_end];
            pRx[k] = pRx[k_end]; pRy[k] = pRy[k_end];
            pwx0[k] = pwx0[k_end]; pwy0[k] = pwy0[k_end];
            pxBar[k] = pxBar[k_end]; pyBar[k] = pyBar[k_end];
            pwx[k] = pwx[k_end]; pwy[k] = pwy[k_end];
            pkxBar[k] = pkxBar[k_end]; pkyBar[k] = pkyBar[k_end];
            pwkx[k] = pwkx[k_end]; pwky[k] = pwky[k_end];
            pz0x2[k] = pz0x2[k_end]; pz0y2[k] = pz0y2[k_end];
            pz0x[k] = pz0x[k_end]; pz0y[k] = pz0y[k_end];
            size_t dst = (size_t)k * Nx * Ny;
            size_t src = (size_t)k_end * Nx * Ny;
            for (size_t i = 0; i < (size_t)Nx * Ny; ++i) {
                pflat[dst + i] = pflat[src + i];
            }
        }
    }

    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "M2_x"), mx_M2x);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "M2_y"), mx_M2y);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "Rx"), mx_Rx);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "Ry"), mx_Ry);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wx0"), mx_wx0);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wy0"), mx_wy0);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "xBar"), mx_xBar);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "yBar"), mx_yBar);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wx"), mx_wx);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wy"), mx_wy);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "kxBar"), mx_kxBar);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "kyBar"), mx_kyBar);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wkx"), mx_wkx);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "wky"), mx_wky);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "z0x2"), mx_z0x2);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "z0y2"), mx_z0y2);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "z0x"), mx_z0x);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "z0y"), mx_z0y);
    mxSetFieldByNumber(S, 0, mxGetFieldNumber(S, "flattened_Exyz"), mx_flat);
}

/* ------------------------------------------------------------------------- */
/* calculate_pulse_flat                                                      */
/* ------------------------------------------------------------------------- */
void run_calculate_pulse_flat(int nlhs, mxArray* plhs[],
                              const mxArray* prhs[]) {
    mxArray* field_tmp = NULL;
    mxComplexDouble* mx_field = get_complex_data_robust(prhs[1], field_tmp);
    const mwSize* dims = mxGetDimensions(prhs[1]);
    size_t ndims = mxGetNumberOfDimensions(prhs[1]);
    if (ndims < 3 || dims[2] < 2) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:pulseFlatField",
            "calculate_pulse_flat requires a 3D field with at least 2 slices.");
    }
    int Nx = (int)dims[0];
    int Ny = (int)dims[1];
    int Nt = (int)dims[2];

    std::vector<double> xvec, yvec;
    extract_grid_vectors(prhs[2], prhs[3], Nx, Ny, xvec, yvec);

    double wavelength = mxGetScalar(prhs[4]);
    double k0 = 2.0 * M_PI / wavelength;

    double dx = (Nx > 1) ? (xvec[1] - xvec[0]) : 1.0;
    double dy = (Ny > 1) ? (yvec[1] - yvec[0]) : 1.0;

    std::vector<double> xmat, ymat;
    fill_meshgrid(Nx, Ny, xvec, yvec, xmat, ymat);

    std::vector<double> kxvec, kyvec;
    build_kvec(Nx, 2.0 * M_PI / dx / Nx, kxvec);
    build_kvec(Ny, 2.0 * M_PI / dy / Ny, kyvec);

    size_t plane_stride = (size_t)Nx * Ny;
    size_t total = plane_stride * Nt;

    std::vector<std::complex<double>> Exyt(total);
    for (size_t i = 0; i < total; ++i) {
        Exyt[i] = std::complex<double>(mx_field[i].real, mx_field[i].imag);
    }

    const int Nrem = 3;
    double xBar_fl[Nrem], yBar_fl[Nrem], kxBar_fl[Nrem], kyBar_fl[Nrem];
    double wxsq_fl[Nrem], wysq_fl[Nrem], wkxsq_fl[Nrem], wkysq_fl[Nrem];
    double Ax[Nrem], Ay[Nrem], Rx[Nrem], Ry[Nrem];

    std::vector<std::complex<double>> dEdx(total);
    std::vector<std::complex<double>> dEdy(total);

    for (int P = 0; P < Nrem; ++P) {
        /* Build mxArray from current Exyt and take the 2D FFT shifted in x and y. */
        mxArray* mx_exyt = create_complex_3d(Nx, Ny, Nt);
        copy_complex_to_mxarray(Exyt.data(), mx_exyt);
        mxArray* mx_k = NULL;
        fft2_shifted(mx_exyt, &mx_k);
        mxDestroyArray(mx_exyt);
        mxComplexDouble* kdata = mxGetComplexDoubles(mx_k);

        /* Spatial fluence and moments over full 3D volume. */
        double U = 0.0;
        double x_sum = 0.0, y_sum = 0.0, x2_sum = 0.0, y2_sum = 0.0;
        for (size_t idx = 0; idx < total; ++idx) {
            double I = std::norm(Exyt[idx]);
            U += I;
            x_sum += xmat[idx % plane_stride] * I;
            y_sum += ymat[idx % plane_stride] * I;
            x2_sum += sq(xmat[idx % plane_stride]) * I;
            y2_sum += sq(ymat[idx % plane_stride]) * I;
        }
        if (U <= 0.0) U = 1.0;
        xBar_fl[P] = x_sum / U;
        yBar_fl[P] = y_sum / U;
        wxsq_fl[P] = 4.0 * (x2_sum / U - sq(xBar_fl[P]));
        wysq_fl[P] = 4.0 * (y2_sum / U - sq(yBar_fl[P]));
        if (wxsq_fl[P] < 0.0) wxsq_fl[P] = 0.0;
        if (wysq_fl[P] < 0.0) wysq_fl[P] = 0.0;

        /* k-space fluence and moments. */
        double Uk = 0.0;
        double kx_sum = 0.0, ky_sum = 0.0, kx2_sum = 0.0, ky2_sum = 0.0;
        for (int t = 0; t < Nt; ++t) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    size_t idx = (size_t)t * plane_stride + (size_t)i + (size_t)j * Nx;
                    double re = kdata[idx].real;
                    double im = kdata[idx].imag;
                    double I = re*re + im*im;
                    Uk += I;
                    kx_sum += kxvec[i] * I;
                    ky_sum += kyvec[j] * I;
                    kx2_sum += sq(kxvec[i]) * I;
                    ky2_sum += sq(kyvec[j]) * I;
                }
            }
        }
        mxDestroyArray(mx_k);

        if (Uk <= 0.0) Uk = 1.0;
        kxBar_fl[P] = kx_sum / Uk;
        kyBar_fl[P] = ky_sum / Uk;
        wkxsq_fl[P] = 4.0 * (kx2_sum / Uk - sq(kxBar_fl[P]));
        wkysq_fl[P] = 4.0 * (ky2_sum / Uk - sq(kyBar_fl[P]));
        if (wkxsq_fl[P] < 0.0) wkxsq_fl[P] = 0.0;
        if (wkysq_fl[P] < 0.0) wkysq_fl[P] = 0.0;

        /* Derivatives of the 3D field. */
        for (int t = 0; t < Nt; ++t) {
            size_t off = (size_t)t * plane_stride;
            derivative_x(&Exyt[off], Nx, Ny, dx, &dEdx[off]);
            derivative_y(&Exyt[off], Nx, Ny, dy, &dEdy[off]);
        }

        /* Ax/Rx using the spatial second moments from the current stage. */
        std::complex<double> accum_x(0.0, 0.0), accum_y(0.0, 0.0);
        for (size_t idx = 0; idx < total; ++idx) {
            std::complex<double> zx = Exyt[idx] * std::conj(dEdx[idx]);
            std::complex<double> zy = Exyt[idx] * std::conj(dEdy[idx]);
            accum_x += (xmat[idx % plane_stride] - xBar_fl[P]) * (zx - std::conj(zx));
            accum_y += (ymat[idx % plane_stride] - yBar_fl[P]) * (zy - std::conj(zy));
        }
        Ax[P] = (std::complex<double>(0.0, 1.0) / k0 * accum_x / U).real();
        Ay[P] = (std::complex<double>(0.0, 1.0) / k0 * accum_y / U).real();
        if (std::abs(Ax[P]) < 1e-40) Ax[P] = 1e-40;
        if (std::abs(Ay[P]) < 1e-40) Ay[P] = 1e-40;
        Rx[P] = wxsq_fl[P] / (2.0 * Ax[P]);
        Ry[P] = wysq_fl[P] / (2.0 * Ay[P]);

        if (P < Nrem - 1) {
            remove_curvature(Exyt.data(), Nx, Ny, Nt, xmat, ymat,
                             xBar_fl[P], yBar_fl[P], Rx[P], Ry[P], k0);
        }
    }

    double Rx_corr = 1.0 / (1.0 / Rx[0] + 1.0 / Rx[1] + 1.0 / Rx[2]);
    double Ry_corr = 1.0 / (1.0 / Ry[0] + 1.0 / Ry[1] + 1.0 / Ry[2]);

    double wx0sq = wxsq_fl[0] - (1.0 / wkxsq_fl[0]) * sq(wxsq_fl[0] * k0 / Rx_corr);
    double wy0sq = wysq_fl[0] - (1.0 / wkysq_fl[0]) * sq(wysq_fl[0] * k0 / Ry_corr);
    if (wx0sq < 0.0) wx0sq = 0.0;
    if (wy0sq < 0.0) wy0sq = 0.0;

    double M2x = std::sqrt(wx0sq) * std::sqrt(wkxsq_fl[0]) / 2.0;
    double M2y = std::sqrt(wy0sq) * std::sqrt(wkysq_fl[0]) / 2.0;
    double wx0 = std::sqrt(wx0sq);
    double wy0 = std::sqrt(wy0sq);

    double z0x_vec[Nrem], z0y_vec[Nrem];
    for (int P = 0; P < Nrem; ++P) {
        z0x_vec[P] = k0 * k0 * wxsq_fl[P] / (Rx[P] * wkxsq_fl[P]);
        z0y_vec[P] = k0 * k0 * wysq_fl[P] / (Ry[P] * wkysq_fl[P]);
    }

    const char* field_names[] = {
        "pulse_M2_x", "pulse_M2_y", "pulse_Rx", "pulse_Ry",
        "pulse_wx0", "pulse_wy0", "pulse_xBar", "pulse_yBar",
        "pulse_wx", "pulse_wy", "pulse_kxBar", "pulse_kyBar",
        "pulse_wkx", "pulse_wky", "pulse_z0x", "pulse_z0y",
        "pulse_flattened_Exyt"
    };
    int nfields = sizeof(field_names) / sizeof(field_names[0]);
    mxArray* S = mxCreateStructMatrix(1, 1, nfields, field_names);
    plhs[0] = S;

    auto set_scalar = [&](const char* name, double val) {
        mxArray* a = mxCreateDoubleScalar(val);
        mxSetField(S, 0, name, a);
    };
    auto set_vec3 = [&](const char* name, const double* v) {
        mxArray* a = mxCreateDoubleMatrix(3, 1, mxREAL);
        double* p = mxGetPr(a);
        p[0] = v[0]; p[1] = v[1]; p[2] = v[2];
        mxSetField(S, 0, name, a);
    };

    set_scalar("pulse_M2_x", M2x);
    set_scalar("pulse_M2_y", M2y);
    set_scalar("pulse_Rx", Rx_corr);
    set_scalar("pulse_Ry", Ry_corr);
    set_scalar("pulse_wx0", wx0);
    set_scalar("pulse_wy0", wy0);
    set_scalar("pulse_xBar", xBar_fl[0]);
    set_scalar("pulse_yBar", yBar_fl[0]);
    set_scalar("pulse_wx", std::sqrt(wxsq_fl[0]));
    set_scalar("pulse_wy", std::sqrt(wysq_fl[0]));
    set_scalar("pulse_kxBar", kxBar_fl[0]);
    set_scalar("pulse_kyBar", kyBar_fl[0]);
    set_scalar("pulse_wkx", std::sqrt(wkxsq_fl[0]));
    set_scalar("pulse_wky", std::sqrt(wkysq_fl[0]));
    set_vec3("pulse_z0x", z0x_vec);
    set_vec3("pulse_z0y", z0y_vec);

    mxArray* mx_flat = create_complex_3d(Nx, Ny, Nt);
    copy_complex_to_mxarray(Exyt.data(), mx_flat);
    mxSetField(S, 0, "pulse_flattened_Exyt", mx_flat);
    if (field_tmp) mxDestroyArray(field_tmp);
}

/* ------------------------------------------------------------------------- */
/* calculate_second_moment                                                   */
/* ------------------------------------------------------------------------- */
void run_calculate_second_moment(int nlhs, mxArray* plhs[],
                                 const mxArray* prhs[]) {
    mxArray* field_tmp = NULL;
    mxComplexDouble* mx_field = get_complex_data_robust(prhs[1], field_tmp);
    const mwSize* dims = mxGetDimensions(prhs[1]);
    size_t ndims = mxGetNumberOfDimensions(prhs[1]);
    int Nx = (int)dims[0];
    int Ny = (int)dims[1];
    int Nt = (ndims > 2) ? (int)dims[2] : 1;

    std::vector<double> xvec, yvec;
    extract_grid_vectors(prhs[2], prhs[3], Nx, Ny, xvec, yvec);

    std::vector<double> xmat, ymat;
    fill_meshgrid(Nx, Ny, xvec, yvec, xmat, ymat);

    mxArray* wx_a = mxCreateDoubleMatrix(Nt, 1, mxREAL);
    mxArray* wy_a = mxCreateDoubleMatrix(Nt, 1, mxREAL);
    double* wx = mxGetPr(wx_a);
    double* wy = mxGetPr(wy_a);

    size_t plane = (size_t)Nx * Ny;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int t = 0; t < Nt; ++t) {
        double I_sum = 0.0;
        double x_sum = 0.0, y_sum = 0.0;
        double x2_sum = 0.0, y2_sum = 0.0;
        for (size_t idx = 0; idx < plane; ++idx) {
            size_t global = t * plane + idx;
            double re = mx_field[global].real;
            double im = mx_field[global].imag;
            double I = re*re + im*im;
            I_sum += I;
            x_sum += xmat[idx] * I;
            y_sum += ymat[idx] * I;
            x2_sum += sq(xmat[idx]) * I;
            y2_sum += sq(ymat[idx]) * I;
        }
        if (I_sum > 0.0) {
            double xBar = x_sum / I_sum;
            double yBar = y_sum / I_sum;
            wx[t] = std::sqrt(4.0 * (x2_sum / I_sum - sq(xBar)));
            wy[t] = std::sqrt(4.0 * (y2_sum / I_sum - sq(yBar)));
        } else {
            wx[t] = 0.0;
            wy[t] = 0.0;
        }
    }

    const char* names[] = {"wx", "wy"};
    mxArray* S = mxCreateStructMatrix(1, 1, 2, names);
    mxSetField(S, 0, "wx", wx_a);
    mxSetField(S, 0, "wy", wy_a);
    plhs[0] = S;
    if (field_tmp) mxDestroyArray(field_tmp);
}

} /* end anonymous namespace */

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    if (nrhs < 5) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:nrhs",
            "Usage: S = msquared_calculate_mex('action', field, xgrid, ygrid, wavelength [, power_threshold]);");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:action",
            "First argument must be a char array action string.");
    }
    char* action_c = mxArrayToString(prhs[0]);
    std::string action(action_c);
    mxFree(action_c);

    if (!mxIsDouble(prhs[1]) || mxIsEmpty(prhs[1])) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:field",
            "field must be a non-empty double array.");
    }
    if (!mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3])) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:grid",
            "xgrid and ygrid must be double.");
    }
    if (!mxIsScalar(prhs[4])) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:wavelength",
            "wavelength must be a scalar.");
    }
    if (nrhs > 5 && !mxIsScalar(prhs[5])) {
        mexErrMsgIdAndTxt("msquared_calculate_mex:powerThreshold",
            "power_threshold must be a scalar.");
    }

    if (action == "calculate") {
        if (nrhs < 6) {
            mexErrMsgIdAndTxt("msquared_calculate_mex:calculateArgs",
                "'calculate' requires power_threshold as sixth argument.");
        }
        run_calculate(nlhs, plhs, prhs);
    } else if (action == "pulse_flat") {
        run_calculate_pulse_flat(nlhs, plhs, prhs);
    } else if (action == "second_moment") {
        run_calculate_second_moment(nlhs, plhs, prhs);
    } else {
        mexErrMsgIdAndTxt("msquared_calculate_mex:unknownAction",
            "Unknown action '%s'.  Use 'calculate', 'pulse_flat', or 'second_moment'.", action.c_str());
    }
}
