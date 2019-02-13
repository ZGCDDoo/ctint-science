#pragma once

#include "Utilities.hpp"
#include "GreenMat.hpp"

using Vertex = Utilities::Vertex;

namespace Models
{

class Model_2D
{

  public:
    Model_2D(const Json &jj) : U_(jj["U"]),
                               delta_(jj["delta"]),
                               beta_(jj["beta"]),
                               mu_(jj["mu"])
    {
        std::cout << "Start Model constructor " << std::endl;

        FinishConstructor(jj);
        std::cout << "End Model constructor " << std::endl;
    };

    void FinishConstructor(const Json &jj)
    {
        const std::string hybNameUp = jj["HybFile"].get<std::string>() + std::string(".arma");
        ClusterMatrixCD_t hybtmpUp;
        assert(hybtmpUp.load(hybNameUp));
        assert(hybtmpUp.n_cols == 1);

        const size_t NNOld = hybtmpUp.n_rows;
        hybtmpUp.resize(jj["NMAT"].get<size_t>(), 1);
        for (size_t nn = NNOld; nn < hybtmpUp.n_rows; nn++)
        {
            const cd_t iwn(0.0, (2.0 * nn + 1) * M_PI / beta_);
            hybtmpUp(nn, 0) = 4.0 / iwn;
        }

        this->hybridizationMatUp_ = hybtmpUp;
        std::cout << "hybtmpUp.real() " << hybtmpUp(0, 0).real() << ",  " << hybtmpUp(0, 0).imag() << std::endl;

        //this is in fact greencluster tilde.
        this->greenCluster0MatUp_ = GreenMat::GreenCluster0Mat(this->hybridizationMatUp_, this->auxMu(), this->beta_);
    }

    ~Model_2D() = default;

    //Getters
    double mu() const { return mu_; };
    double U() const { return U_; };
    double delta() const { return delta_; };
    double beta() const { return beta_; };
    GreenMat::GreenCluster0Mat const greenCluster0MatUp() { return greenCluster0MatUp_; };
    ClusterMatrixCD_t const hybridizationMatUp() const { return hybridizationMatUp_; };

    double auxUp(const Vertex vertex) const { return ((vertex.aux() == Up) ? 1.0 + delta_ : -delta_); };
    double auxDown(const Vertex vertex) const { return ((vertex.aux() == Down) ? 1.0 + delta_ : -delta_); };
    double auxU() const { return U_ / 2.0; };
    double auxMu() const { return mu_ - U_ / 2.0; };
    double auxDO() const { return delta_ * (1.0 + delta_); };

  protected:
    const double U_;
    const double delta_;
    const double beta_;
    double mu_;

    ClusterMatrixCD_t hybridizationMatUp_;
    GreenMat::GreenCluster0Mat greenCluster0MatUp_;
};

} // namespace Models
