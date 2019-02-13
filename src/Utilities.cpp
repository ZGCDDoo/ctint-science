
#include "ctint-science/Utilities.hpp"

namespace Utilities
{

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

Vertex::Vertex(Tau_t tau, AuxSpin_t aux) : tau_(tau), aux_(aux){};
// Getters
Tau_t Vertex::tau() const { return tau_; };
AuxSpin_t Vertex::aux() const { return aux_; };

void Vertex::FlipAux() { aux_ == Up ? aux_ = Down : aux_ = Up; };

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