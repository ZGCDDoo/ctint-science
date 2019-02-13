#pragma once
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
//#define ARMA_NO_DEBUG

#include <boost/random.hpp>
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

double DotVectors(const SiteVector_t &v1, const SiteVector_t &v2);

void MatrixVectorMult(const ClusterMatrix_t &A, const SiteVector_t &X, const double &alpha, SiteVector_t &Y);

void VectorMatrixMult(SiteVector_t const &X, ClusterMatrix_t const &A, const double &alpha, SiteVector_t &Y);

double Dot(const SiteVector_t &v1, const ClusterMatrix_t &A, const SiteVector_t &v2);

class Vertex
{

  public:
    Vertex(Tau_t tau, AuxSpin_t aux);
    // Getters
    Tau_t tau() const;
    AuxSpin_t aux() const;

    void FlipAux();

  private:
    Tau_t tau_;
    AuxSpin_t aux_;
};

//Upgrade the matrix if the last element of the inverse is known (STilde)
void BlockRankOneUpgrade(ClusterMatrix_t &mk, const SiteVector_t &Q, const SiteVector_t &R, const double &STilde);

void BlockRankOneDowngrade(ClusterMatrix_t &m1, ClusterMatrix_t &dummy);
} // namespace Utilities