// Taken from LCG ROOT MathLib
// License info:
// Authors: Andras Zsenei & Lorenzo Moneta   06/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Langaus authors:
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//  Markus Friedl (Markus.Friedl@cern.ch)
//
//  Adaption for Python by David-Leon Pohl, pohl@physik.uni-bonn.de
//  Adaption for pybind11 by Bong-Hwi Lim, bong-hwi.lim@cern.ch

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <iostream>

namespace py = pybind11;
// Constants
const double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)

// Function prototypes
double landauPDF(const double& x, const double& x0, const double& xi)  // same algorithm is used in GSL
{
	static double p1[5] =
	{ 0.4259894875, -0.1249762550, 0.03984243700, -0.006298287635, 0.001511162253 };
	static double q1[5] =
	{ 1.0, -0.3388260629, 0.09594393323, -0.01608042283, 0.003778942063 };

	static double p2[5] =
	{ 0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411, 0.0001283617211 };
	static double q2[5] =
	{ 1.0, 0.7428795082, 0.3153932961, 0.06694219548, 0.008790609714 };

	static double p3[5] =
	{ 0.1788544503, 0.09359161662, 0.006325387654, 0.00006611667319, -0.000002031049101 };
	static double q3[5] =
	{ 1.0, 0.6097809921, 0.2560616665, 0.04746722384, 0.006957301675 };

	static double p4[5] =
	{ 0.9874054407, 118.6723273, 849.2794360, -743.7792444, 427.0262186 };
	static double q4[5] =
	{ 1.0, 106.8615961, 337.6496214, 2016.712389, 1597.063511 };

	static double p5[5] =
	{ 1.003675074, 167.5702434, 4789.711289, 21217.86767, -22324.94910 };
	static double q5[5] =
	{ 1.0, 156.9424537, 3745.310488, 9834.698876, 66924.28357 };

	static double p6[5] =
	{ 1.000827619, 664.9143136, 62972.92665, 475554.6998, -5743609.109 };
	static double q6[5] =
	{ 1.0, 651.4101098, 56974.73333, 165917.4725, -2815759.939 };

	static double a1[3] =
	{ 0.04166666667, -0.01996527778, 0.02709538966 };

	static double a2[2] =
	{ -1.845568670, -4.284640743 };

	if (xi <= 0)
		return 0;
	double v = (x - x0) / xi;
	double u, ue, us, denlan;
	if (v < -5.5) {
		u = std::exp(v + 1.0);
		if (u < 1e-10)
			return 0.0;
		ue = std::exp(-1 / u);
		us = std::sqrt(u);
		denlan = 0.3989422803 * (ue / us) * (1 + (a1[0] + (a1[1] + a1[2] * u) * u) * u);
	}
	else if (v < -1) {
		u = std::exp(-v - 1);
		denlan = std::exp(-u) * std::sqrt(u) * (p1[0] + (p1[1] + (p1[2] + (p1[3] + p1[4] * v) * v) * v) * v) / (q1[0] + (q1[1] + (q1[2] + (q1[3] + q1[4] * v) * v) * v) * v);
	}
	else if (v < 1) {
		denlan = (p2[0] + (p2[1] + (p2[2] + (p2[3] + p2[4] * v) * v) * v) * v) / (q2[0] + (q2[1] + (q2[2] + (q2[3] + q2[4] * v) * v) * v) * v);
	}
	else if (v < 5) {
		denlan = (p3[0] + (p3[1] + (p3[2] + (p3[3] + p3[4] * v) * v) * v) * v) / (q3[0] + (q3[1] + (q3[2] + (q3[3] + q3[4] * v) * v) * v) * v);
	}
	else if (v < 12) {
		u = 1 / v;
		denlan = u * u * (p4[0] + (p4[1] + (p4[2] + (p4[3] + p4[4] * u) * u) * u) * u) / (q4[0] + (q4[1] + (q4[2] + (q4[3] + q4[4] * u) * u) * u) * u);
	}
	else if (v < 50) {
		u = 1 / v;
		denlan = u * u * (p5[0] + (p5[1] + (p5[2] + (p5[3] + p5[4] * u) * u) * u) * u) / (q5[0] + (q5[1] + (q5[2] + (q5[3] + q5[4] * u) * u) * u) * u);
	}
	else if (v < 300) {
		u = 1 / v;
		denlan = u * u * (p6[0] + (p6[1] + (p6[2] + (p6[3] + p6[4] * u) * u) * u) * u) / (q6[0] + (q6[1] + (q6[2] + (q6[3] + q6[4] * u) * u) * u) * u);
	}
	else {
		u = 1 / (v - v * std::log(v) / (v + 1));
		denlan = u * u * (1 + (a2[0] + a2[1] * u) * u);
	}
	return denlan / xi;
}

double gaussPDF(const double& x, const double& mu, const double& sigma)
{
	return invsq2pi / sigma * exp(- pow((x - mu), 2) / (2 * pow(sigma, 2)));
}

double landauGaussPDF(const double& x, const double& mu, const double& eta, const double& sigma)
{
	// Numeric constants
	double mpshift = 0.; // -0.22278298;     // Landau maximum location shift in original code is wrong, since the shift does not depend on mu only

	// Control constants
	unsigned int np = 100;      // number of convolution steps
	const unsigned int sc = 8;        // convolution extends to +-sc Gaussian sigmas

	// Convolution steps have to be increased if sigma > eta * 5 to get stable solution that does not oscillate, addresses #1
	if (sigma > 3 * eta)
		np *= int(sigma / eta / 3.);
	if (np > 100000)  // Do not use too many convolution steps to save time
		np = 100000;

	// Variables
	double xx;
	double mpc;
	double fland;
	double sum = 0.0;
	double xlow, xupp;
	double step;

	// MP shift correction
	mpc = mu - mpshift;// * eta;

	// Range of convolution integral
	xlow = x - sc * sigma;
	xupp = x + sc * sigma;

	step = (xupp - xlow) / (double) np;

	// Discrete linear convolution of Landau and Gaussian
	for (unsigned int i = 1; i <= np / 2; ++i) {
		xx = xlow + double (i - 0.5) * step;
		fland = landauPDF(xx, mpc, eta) / eta;
		sum += fland * gaussPDF(x, xx, sigma);

		xx = xupp - double (i - 0.5) * step;
		fland = landauPDF(xx, mpc, eta) / eta;
		sum += fland * gaussPDF(x, xx, sigma);
	}

	const double norm = 0.398902310115109;  // normalization of the integral from -inf..intf

	return (step * sum);
}

// Modified to take a simple pointer (double* data) rather than a reference to a pointer (double*& data)
// This simplifies compatibility with pybind11 by reading data directly without modifying the original pointer
double* getLandauPDFData(double* data, const unsigned int& size, const double& mu, const double& eta) {
	// 
    double* result = new double[size];
    for (unsigned int i = 0; i < size; ++i) {
        result[i] = landauPDF(data[i], mu, eta);
    }
    return result;
}

// Modified to take a simple pointer (double* data) rather than a reference to a pointer (double*& data)
// This simplifies compatibility with pybind11 by reading data directly without modifying the original pointer
double* getLangauPDFData(double* data, const unsigned int& size, const double& mu, const double& eta, const double& sigma) {
    double* result = new double[size];
    for (unsigned int i = 0; i < size; ++i) {
        result[i] = landauGaussPDF(data[i], mu, eta, sigma);
    }
    return result;
}

// Helper function to convert C++ array to NumPy array
py::array_t<double> data_to_numpy_array_double(double* ptr, size_t N) {
    py::array_t<double> arr(N, ptr);
    delete[] ptr; // Clean up allocated memory
    return arr;
}

// Helper function to check and adjust parameters
std::tuple<double, double, double, double> check_parameter(double mpv, double eta, double sigma, double A = 1.0) {
    if (eta < 1e-9) {
        std::cerr << "Warning: eta < 1e-9 is not supported. Setting eta to 1e-9." << std::endl;
        eta = 1e-9;
    }
    if (sigma < 0) sigma *= -1;
    if (sigma > 100 * eta) {
        std::cerr << "Warning: sigma > 100 * eta can lead to oscillations. Check result." << std::endl;
    }
    if (A < 0.0) throw std::invalid_argument("A must be >= 0");
    return std::make_tuple(mpv, eta, sigma, A);
}

// Placeholder _scale_to_mpv function for optimization
std::tuple<double, double, double, double> scale_to_mpv(double mu, double eta, double sigma = 0.0, double A = 1.0) {
    double mpv_scaled = mu;  // Placeholder: replace with actual optimization if necessary
    double A_scaled = A;
    return std::make_tuple(mpv_scaled, eta, sigma, A_scaled);
}

// Single-value PDF functions
double get_landau_pdf(double value, double mu = 0, double eta = 1) {
    return landauPDF(value, mu, eta);
}

double get_gauss_pdf(double value, double mu = 0, double sigma = 1) {
    return gaussPDF(value, mu, sigma);
}

double get_langau_pdf(double value, double mu = 0, double eta = 1, double sigma = 1) {
    return landauGaussPDF(value, mu, eta, sigma);
}

// Array-based Landau PDF computation
py::array_t<double> landau_pdf(py::array_t<double> input, double mu = 0, double eta = 1) {
    auto buf = input.request();
    double* ptr = static_cast<double*>(buf.ptr);
    size_t size = buf.size;
    double* result = getLandauPDFData(ptr, size, mu, eta);
    return data_to_numpy_array_double(result, size);
}

// Array-based Langau PDF computation
py::array_t<double> langau_pdf(py::array_t<double> input, double mu = 0, double eta = 1, double sigma = 1) {
    auto buf = input.request();
    double* ptr = static_cast<double*>(buf.ptr);
    size_t size = buf.size;
    double* result = getLangauPDFData(ptr, size, mu, eta, sigma);
    return data_to_numpy_array_double(result, size);
}

// The `landau` function that applies amplitude scaling
py::array_t<double> landau(py::array_t<double> array, double mpv = 0, double eta = 1, double A = 1) {
    if (A == 0) return py::array_t<double>(array.size()); // Return zeroed array if A == 0

    auto [mpv_corrected, eta_corrected, sigma_corrected, A_corrected] = check_parameter(mpv, eta, 0, A);
    auto [mpv_scaled, eta_scaled, sigma_scaled, A_scaled] = scale_to_mpv(mpv_corrected, eta_corrected, 0, A_corrected);

    py::array_t<double> result = landau_pdf(array, mpv_scaled, eta_scaled);
    auto result_mut = result.mutable_unchecked<1>(); // Specify 1 for 1D array

    for (ssize_t i = 0; i < result_mut.size(); ++i) {
        result_mut(i) *= A_scaled;
    }
    return result;
}

// The `langau` function that applies amplitude scaling
py::array_t<double> langau(py::array_t<double> array, double mpv = 0, double eta = 1, double sigma = 1, double A = 1, bool scale_langau = true) {
    if (A == 0) return py::array_t<double>(array.size()); // Return zeroed array if A == 0

    auto [mpv_corrected, eta_corrected, sigma_corrected, A_corrected] = check_parameter(mpv, eta, sigma, A);

    double mpv_scaled, A_scaled;
    if (scale_langau) {
        std::tie(mpv_scaled, std::ignore, std::ignore, A_scaled) = scale_to_mpv(mpv_corrected, eta_corrected, sigma_corrected, A_corrected);
    } else {
        std::tie(mpv_scaled, std::ignore, std::ignore, A_scaled) = scale_to_mpv(mpv_corrected, eta_corrected, 0, A_corrected);
    }

    py::array_t<double> result = langau_pdf(array, mpv_scaled, eta_corrected, sigma_corrected);
    auto result_mut = result.mutable_unchecked<1>(); // Specify 1 for 1D array

    for (ssize_t i = 0; i < result_mut.size(); ++i) {
        result_mut(i) *= A_scaled;
    }
    return result;
}

// Define the pybind11 module
PYBIND11_MODULE(pylandau, m) {
    m.def("get_landau_pdf", &get_landau_pdf, "Single value Landau PDF", py::arg("value"), py::arg("mu") = 0, py::arg("eta") = 1);
    m.def("get_gauss_pdf", &get_gauss_pdf, "Single value Gaussian PDF", py::arg("value"), py::arg("mu") = 0, py::arg("sigma") = 1);
    m.def("get_langau_pdf", &get_langau_pdf, "Single value Landau-Gauss PDF", py::arg("value"), py::arg("mu") = 0, py::arg("eta") = 1, py::arg("sigma") = 1);
    m.def("landau_pdf", &landau_pdf, "Array-based Landau PDF", py::arg("input"), py::arg("mu") = 0, py::arg("eta") = 1);
    m.def("langau_pdf", &langau_pdf, "Array-based Langau PDF", py::arg("input"), py::arg("mu") = 0, py::arg("eta") = 1, py::arg("sigma") = 1);
    m.def("check_parameter", &check_parameter, "Check and adjust parameters", py::arg("mpv"), py::arg("eta"), py::arg("sigma"), py::arg("A") = 1.0);
    m.def("scale_to_mpv", &scale_to_mpv, "Scale function to maximum probable value", py::arg("mu"), py::arg("eta"), py::arg("sigma") = 0.0, py::arg("A") = 1.0);
    m.def("landau", &landau, "Higher-level Landau function", py::arg("array"), py::arg("mpv") = 0, py::arg("eta") = 1, py::arg("A") = 1);
    m.def("langau", &langau, "Higher-level Langau function", py::arg("array"), py::arg("mpv") = 0, py::arg("eta") = 1, py::arg("sigma") = 1, py::arg("A") = 1, py::arg("scale_langau") = true);
}
