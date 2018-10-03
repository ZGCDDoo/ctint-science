#pragma once
#include "Utilities.hpp"

namespace GreenMat
{

class GreenCluster0Mat
{

  public:
    GreenCluster0Mat() : hyb_(),
                         mu_(),
                         beta_(), zm_(), fm_(), sm_(), tm_(), data_(){};

    GreenCluster0Mat(const GreenCluster0Mat &gf) : hyb_(gf.hyb_),
                                                   mu_(gf.mu_),
                                                   beta_(gf.beta_), zm_(gf.zm_), fm_(gf.fm_), sm_(gf.sm_), tm_(gf.tm_), data_(gf.data_){};

    GreenCluster0Mat(const ClusterMatrixCD_t &hyb, const double &mu, const double &beta) : hyb_(hyb),
                                                                                           mu_(mu),
                                                                                           beta_(beta),
                                                                                           data_()
    {
        std::cout << "Start GreenMat constructor " << std::endl;

        const size_t Nc = 1;
        const size_t NN = hyb.n_rows;

        data_.resize(NN, Nc);
        data_.zeros();

        zm_ = cd_t(0.0);
        fm_ = cd_t(1.0, 0.0);
        sm_ = cd_t(-mu_, 0.0);
        tm_ = cd_t(mu_ * mu_ + 4.0, 0.0);

        for (size_t n = 0; n < NN; n++)
        {
            const cd_t zz = cd_t(mu_, (2.0 * n + 1.0) * M_PI / beta_);
            const cd_t tmp = zz - hyb_(n, 0);
            data_(n, 0) = 1.0 / tmp;
        }

        Save();
        std::cout << "End GreenMat constructor " << std::endl;
    }

    ~GreenCluster0Mat() = default;

    const GreenCluster0Mat &operator=(const GreenCluster0Mat &gf)
    {
        if (this == &gf)
            return *this; //évite les boucles infinies
        data_ = gf.data_;
        zm_ = gf.zm_;
        fm_ = gf.fm_;
        sm_ = gf.sm_;
        tm_ = gf.tm_;
        mu_ = gf.mu_;
        beta_ = gf.beta_;
        hyb_ = gf.hyb_;
        return *this;
    }

    ClusterMatrixCD_t data() const { return data_; };
    cd_t zm() const { return zm_; };
    cd_t fm() const { return fm_; };
    cd_t sm() const { return sm_; };
    cd_t tm() const { return tm_; };
    double beta() const { return beta_; };
    size_t n_rows() const { return data_.n_rows; };
    size_t n_cols() const { return data_.n_cols; };

    void Save()
    {
        std::ofstream fout("g0iwn.dat");
        for (size_t nn = 0; nn < data_.n_rows; nn++)
        {
            const double wn = (2.0 * nn + 1.0) * M_PI / beta_;
            fout << wn << " " << data_(nn, 0).real() << " " << data_(nn, 0).imag() << std::endl;
        }
        fout.close();
    }

  private:
    ClusterMatrixCD_t hyb_;
    double mu_;
    double beta_;

    cd_t zm_;
    cd_t fm_;
    cd_t sm_;
    cd_t tm_;
    ClusterMatrixCD_t data_;
};
} // namespace GreenMat