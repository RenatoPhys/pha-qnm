#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>

namespace pha_qnm {

template <std::size_t Order>
struct Series {
  std::array<double, Order + 1> c{};

  Series() = default;
  Series(double value) { c[0] = value; }

  double& operator[](std::size_t i) { return c.at(i); }
  const double& operator[](std::size_t i) const { return c.at(i); }
};

template <std::size_t N>
Series<N> operator+(const Series<N>& a, const Series<N>& b) {
  Series<N> out;
  for (std::size_t i = 0; i <= N; ++i) out[i] = a[i] + b[i];
  return out;
}

template <std::size_t N>
Series<N> operator-(const Series<N>& a, const Series<N>& b) {
  Series<N> out;
  for (std::size_t i = 0; i <= N; ++i) out[i] = a[i] - b[i];
  return out;
}

template <std::size_t N>
Series<N> operator-(const Series<N>& a) {
  Series<N> out;
  for (std::size_t i = 0; i <= N; ++i) out[i] = -a[i];
  return out;
}

template <std::size_t N>
Series<N> operator*(const Series<N>& a, const Series<N>& b) {
  Series<N> out;
  for (std::size_t i = 0; i <= N; ++i)
    for (std::size_t j = 0; j + i <= N; ++j) out[i + j] += a[i] * b[j];
  return out;
}

template <std::size_t N>
Series<N> operator+(const Series<N>& a, double b) { return a + Series<N>(b); }
template <std::size_t N>
Series<N> operator+(double a, const Series<N>& b) { return Series<N>(a) + b; }
template <std::size_t N>
Series<N> operator-(const Series<N>& a, double b) { return a - Series<N>(b); }
template <std::size_t N>
Series<N> operator-(double a, const Series<N>& b) { return Series<N>(a) - b; }
template <std::size_t N>
Series<N> operator*(const Series<N>& a, double b) {
  Series<N> out;
  for (std::size_t i = 0; i <= N; ++i) out[i] = a[i] * b;
  return out;
}
template <std::size_t N>
Series<N> operator*(double a, const Series<N>& b) { return b * a; }

template <std::size_t N>
Series<N> inverse(const Series<N>& a) {
  if (a[0] == 0.0) throw std::domain_error("series inverse with zero constant");
  Series<N> out;
  out[0] = 1.0 / a[0];
  for (std::size_t n = 1; n <= N; ++n) {
    double sum = 0.0;
    for (std::size_t k = 1; k <= n; ++k) sum += a[k] * out[n - k];
    out[n] = -sum / a[0];
  }
  return out;
}

template <std::size_t N>
Series<N> operator/(const Series<N>& a, const Series<N>& b) { return a * inverse(b); }
template <std::size_t N>
Series<N> operator/(const Series<N>& a, double b) { return a * (1.0 / b); }
template <std::size_t N>
Series<N> operator/(double a, const Series<N>& b) { return Series<N>(a) / b; }

template <std::size_t N>
Series<N> derivative(const Series<N>& a) {
  Series<N> out;
  for (std::size_t i = 1; i <= N; ++i) out[i - 1] = static_cast<double>(i) * a[i];
  return out;
}

template <std::size_t N>
double evaluate(const Series<N>& a, double x) {
  double value = a[N];
  for (std::size_t i = N; i-- > 0;) value = value * x + a[i];
  return value;
}

namespace math {

inline double exp(double x) { return std::exp(x); }
inline double sinh(double x) { return std::sinh(x); }
inline double cosh(double x) { return std::cosh(x); }
inline double tanh(double x) { return std::tanh(x); }
inline double sech(double x) {
  const double e = std::exp(-std::abs(x));
  return 2.0 * e / (1.0 + e * e);
}

inline std::complex<double> exp(const std::complex<double>& x) { return std::exp(x); }
inline std::complex<double> sinh(const std::complex<double>& x) { return std::sinh(x); }
inline std::complex<double> cosh(const std::complex<double>& x) { return std::cosh(x); }
inline std::complex<double> tanh(const std::complex<double>& x) { return std::tanh(x); }
inline std::complex<double> sech(const std::complex<double>& x) {
  return 1.0 / std::cosh(x);
}

template <std::size_t N>
Series<N> exp(const Series<N>& x) {
  Series<N> y;
  y[0] = std::exp(x[0]);
  for (std::size_t n = 1; n <= N; ++n) {
    double sum = 0.0;
    for (std::size_t k = 1; k <= n; ++k)
      sum += static_cast<double>(k) * x[k] * y[n - k];
    y[n] = sum / static_cast<double>(n);
  }
  return y;
}

template <std::size_t N>
Series<N> sinh(const Series<N>& x) { return 0.5 * (exp(x) - exp(-x)); }

template <std::size_t N>
Series<N> cosh(const Series<N>& x) { return 0.5 * (exp(x) + exp(-x)); }

template <std::size_t N>
Series<N> tanh(const Series<N>& x) {
  if (x[0] > 20.0) {
    const auto e = exp(-2.0 * x);
    return (Series<N>(1.0) - e) / (Series<N>(1.0) + e);
  }
  if (x[0] < -20.0) return -tanh(-x);
  return sinh(x) / cosh(x);
}

template <std::size_t N>
Series<N> sech(const Series<N>& x) {
  if (x[0] > 40.0) {
    const auto e = exp(-x);
    return 2.0 * e / (Series<N>(1.0) + e * e);
  }
  if (x[0] < -40.0) return sech(-x);
  return inverse(cosh(x));
}

} // namespace math

} // namespace pha_qnm
