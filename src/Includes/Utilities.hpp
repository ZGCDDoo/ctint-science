#pragma once
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
//#define ARMA_NO_DEBUG

#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>
#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include <fstream>
#include <cassert>
// #include <complex>
//#include <complex.h>
#include <boost/random.hpp>
//#include <boost/filesystem.hpp>
#include <mm_malloc.h>
#include <vector>
#include <string>
#include <chrono>
#include <ccomplex>
//External Libraries
#include <armadillo>
#include "json.hpp"

using Json = nlohmann::json;

//Inspired by Patrick Sémon

typedef std::complex<double> cd_t;
typedef int Sign_t;
typedef size_t Site_t;
typedef double Tau_t;

enum AuxSpin_t
{
    Up,
    Down
};

const size_t DIMENSION = 2;
using SiteVectorCD_t = arma::cx_vec;
using SiteRowCD_t = arma::cx_rowvec;
using ClusterSitesCD_t = std::vector<arma::cx_vec>;
using ClusterMatrixCD_t = arma::cx_mat;
using ClusterCubeCD_t = arma::cx_cube;

using SiteVector_t = arma::vec;
using SiteRow_t = arma::rowvec;
using ClusterSites_t = std::vector<arma::vec>;
using ClusterMatrix_t = arma::mat;
using ClusterCube_t = arma::cube;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

namespace Utilities
{

extern "C"
{
    unsigned int dger_(unsigned int const *, unsigned int const *, double const *, double const *, unsigned int const *, double const *, unsigned int const *, double *, unsigned int const *);
    double ddot_(unsigned int const *, double const *, unsigned int const *, double const *, unsigned int const *);
    unsigned int dgemv_(char const *, unsigned int const *, unsigned int const *, double const *, double const *, unsigned int const *, double const *, unsigned int const *, double const *, double *, unsigned int const *);
    unsigned int dcopy_(unsigned int const *, double const *, unsigned int const *, double *, unsigned int const *);
}

typedef boost::mt19937 EngineType_t;
typedef boost::uniform_real<double> UniformDistribution_t;
typedef boost::variate_generator<EngineType_t &, UniformDistribution_t> UniformRng_t;

using arma::dot;
using arma::kron;

double DotVectors(const SiteVector_t &v1, const SiteVector_t &v2)
{
    unsigned int N = v1.n_elem;
    assert(v1.n_elem == v2.n_elem);
    unsigned int inc = 1;
    return ddot_(&N, v1.memptr(), &inc, v2.memptr(), &inc);
}

void MatrixVectorMult(const ClusterMatrix_t &A, const SiteVector_t &X, const double &alpha, SiteVector_t &Y)
{
    unsigned int N = A.n_cols;
    assert(N == X.n_elem);
    double beta = 0.0;
    char trans = 'n';
    unsigned int inc = 1;
    dgemv_(&trans, &N, &N, &alpha, A.memptr(), &N, X.memptr(), &inc, &beta, Y.memptr(), &inc);
}

void VectorMatrixMult(SiteVector_t const &X, ClusterMatrix_t const &A, const double &alpha, SiteVector_t &Y)
{
    unsigned int N = A.n_cols;
    assert(N == X.n_elem);
    double beta = 0.0;
    char trans = 't';
    unsigned int inc = 1;
    dgemv_(&trans, &N, &N, &alpha, A.memptr(), &N, X.memptr(), &inc, &beta, Y.memptr(), &inc);
}

double Dot(const SiteVector_t &v1, const ClusterMatrix_t &A, const SiteVector_t &v2)
{
    const size_t N = A.n_cols;
    SiteVector_t dummy(N);
    double one = 1.0;
    MatrixVectorMult(A, v2, one, dummy);
    return DotVectors(v1, dummy);
}

class Vertex
{

  public:
    Vertex(Tau_t tau, AuxSpin_t aux) : tau_(tau), aux_(aux){};
    // Getters
    Tau_t tau() const { return tau_; };
    AuxSpin_t aux() const { return aux_; };

    void FlipAux() { aux_ == Up ? aux_ = Down : aux_ = Up; };

  private:
    Tau_t tau_;
    AuxSpin_t aux_;
};

//Upgrade the matrix if the last element of the inverse is known (STilde)
void BlockRankOneUpgrade(ClusterMatrix_t &mk, const SiteVector_t &Q, const SiteVector_t &R, const double &STilde)
{
    const unsigned int k = mk.n_cols;
    const unsigned int kp1 = k + 1;
    const double one = 1.0;

    SiteVector_t mkQ(k);
    SiteVector_t Rmk(k);
    MatrixVectorMult(mk, Q, one, mkQ);
    VectorMatrixMult(R, mk, one, Rmk);

    SiteVector_t QTilde = -STilde * mkQ;
    SiteVector_t RTilde = -STilde * Rmk;

    const unsigned int inc = 1;
    dger_(&k, &k, &STilde, &(mkQ.memptr()[0]), &inc, &(Rmk.memptr()[0]), &inc, mk.memptr(), &k);

    //There is probably a faster way to recreate the matrix. m1.submat ? Test this. in fact use blas dcopy directly?

    mk.resize(kp1, kp1);
    dcopy_(&k, &(RTilde.memptr()[0]), &inc, &(mk.memptr()[k]), &kp1);
    dcopy_(&k, &(QTilde.memptr()[0]), &inc, &(mk.memptr()[kp1 * k]), &inc);

    mk(k, k) = STilde;
    return;
}

void BlockRankOneDowngrade(ClusterMatrix_t &m1, ClusterMatrix_t &dummy) //return dummy
{
    const unsigned int inc = 1;
    const unsigned int kk = m1.n_rows;
    const unsigned int kkm1 = kk - 1;
    dummy.set_size(kk - 1, kk - 1);

    if (kkm1 == 0)
    {
        m1.clear();
    }
    else
    {
        dummy = m1.submat(0, 0, kkm1 - 1, kkm1 - 1);
        double *m1_mem = m1.memptr();
        double alpha = -1.0 / m1(kkm1, kkm1);
        // this next line does S = S - A12 A22^(-1) A21;
        dger_(&kkm1, &kkm1, &alpha, &(m1_mem[kk * kkm1]), &inc, &(m1_mem[kkm1]), &kk, dummy.memptr(), &kkm1);
    }
}

} // namespace Utilities