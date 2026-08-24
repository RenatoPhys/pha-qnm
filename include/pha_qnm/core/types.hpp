#pragma once

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace pha_qnm {

using Real = double;
using Complex = std::complex<double>;

template <class T>
class DenseMatrix {
public:
  DenseMatrix() = default;
  DenseMatrix(std::size_t rows, std::size_t cols, T value = {})
      : rows_(rows), cols_(cols), values_(rows * cols, value) {}

  [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
  [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
  [[nodiscard]] T* data() noexcept { return values_.data(); }
  [[nodiscard]] const T* data() const noexcept { return values_.data(); }

  T& operator()(std::size_t row, std::size_t col) {
    if (row >= rows_ || col >= cols_) throw std::out_of_range("DenseMatrix index");
    return values_[col * rows_ + row];
  }
  const T& operator()(std::size_t row, std::size_t col) const {
    if (row >= rows_ || col >= cols_) throw std::out_of_range("DenseMatrix index");
    return values_[col * rows_ + row];
  }

private:
  std::size_t rows_{};
  std::size_t cols_{};
  std::vector<T> values_;
};

using RealMatrix = DenseMatrix<Real>;
using ComplexMatrix = DenseMatrix<Complex>;

} // namespace pha_qnm

