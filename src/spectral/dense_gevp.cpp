#include "pha_qnm/spectral/dense_gevp.hpp"

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <vector>

extern "C" {
void zggev_(const char* jobvl, const char* jobvr, const int* n,
            std::complex<double>* a, const int* lda, std::complex<double>* b,
            const int* ldb, std::complex<double>* alpha,
            std::complex<double>* beta, std::complex<double>* vl, const int* ldvl,
            std::complex<double>* vr, const int* ldvr, std::complex<double>* work,
            const int* lwork, double* rwork, int* info);
}

namespace pha_qnm {

GeneralizedEigenResult generalized_eigenvalues(ComplexMatrix M0, ComplexMatrix M1,
                                               bool vectors) {
  if (M0.rows() != M0.cols() || M1.rows() != M0.rows() || M1.cols() != M0.cols())
    throw std::invalid_argument("generalized eigensolver requires equal square matrices");
  const int n = static_cast<int>(M0.rows());
  if (static_cast<std::size_t>(n) != M0.rows()) throw std::overflow_error("matrix too large for LAPACK int");
  GeneralizedEigenResult out{{}, {}, ComplexMatrix(vectors ? n : 1, vectors ? n : 1),
                              ComplexMatrix(vectors ? n : 1, vectors ? n : 1)};
  out.alpha.resize(static_cast<std::size_t>(n));
  out.beta.resize(static_cast<std::size_t>(n));
  const char job = vectors ? 'V' : 'N';
  const int ld = std::max(1, n);
  const int ldv = vectors ? ld : 1;
  std::vector<double> rwork(static_cast<std::size_t>(8 * std::max(1, n)));
  Complex query{};
  int lwork = -1;
  int info{};
  zggev_(&job, &job, &n, M0.data(), &ld, M1.data(), &ld, out.alpha.data(),
         out.beta.data(), out.left_eigenvectors.data(), &ldv,
         out.right_eigenvectors.data(), &ldv, &query, &lwork, rwork.data(), &info);
  if (info != 0) throw std::runtime_error("LAPACK ZGGEV workspace query failed");
  lwork = std::max(2 * n, static_cast<int>(query.real()));
  std::vector<Complex> work(static_cast<std::size_t>(std::max(1, lwork)));
  zggev_(&job, &job, &n, M0.data(), &ld, M1.data(), &ld, out.alpha.data(),
         out.beta.data(), out.left_eigenvectors.data(), &ldv,
         out.right_eigenvectors.data(), &ldv, work.data(), &lwork, rwork.data(), &info);
  if (info < 0) throw std::runtime_error("LAPACK ZGGEV received an invalid argument");
  if (info > 0) throw std::runtime_error("LAPACK QZ iteration failed to converge");
  return out;
}

} // namespace pha_qnm

