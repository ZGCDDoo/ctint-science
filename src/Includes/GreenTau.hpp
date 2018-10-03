#pragma once

#include "GreenMat.hpp"
#include "Fourier.hpp"
#include "Utilities.hpp"

namespace GreenTau
{

using namespace GreenMat;

class GreenCluster0Tau
{

  public:
    GreenCluster0Tau() : gfMatCluster_(), NTau_(), beta_(){};
    GreenCluster0Tau(const GreenCluster0Mat &gfMatCluster, const size_t &NTau) : gfMatCluster_(gfMatCluster), NTau_(NTau), beta_(gfMatCluster.beta())
    {
        data_ = SiteVector_t(NTau).zeros();

        for (size_t tt = 0; tt < NTau_; tt++)
        {
            const Tau_t tau = gfMatCluster.beta() * (static_cast<double>(tt) + 1.0) / NTau_;
            double fm = gfMatCluster_.fm().real();
            double sm = gfMatCluster_.sm().real();
            double tm = gfMatCluster_.tm().real();
            data_(tt) = Fourier::MatToTauAnalytic(gfMatCluster_.data().col(0), tau, beta_, fm, sm, tm);
        }

        Save();
    };

    ~GreenCluster0Tau() = default;

    double operator()(Tau_t tau)
    {
        double aps = 1.0;
        if (tau < 0.0)
        {
            tau += beta_;
            aps = -1.0;
        }

        double nt = tau / beta_ * (static_cast<double>(NTau_) - 1.0);
        size_t n0 = static_cast<size_t>(nt);
        double greentau0 = aps * ((1.0 - (nt - n0)) * data_(n0) + (nt - n0) * data_(n0 + 1));
        return greentau0;
    }

    const GreenCluster0Tau &operator=(const GreenCluster0Tau &gf)
    {
        if (this == &gf)
            return *this; //évite les boucles infinies
        gfMatCluster_ = gf.gfMatCluster_;
        NTau_ = gf.NTau_;
        beta_ = gf.beta_;
        data_ = gf.data_;
        return *this;
    }

    void Save()
    {
        std::ofstream fout("g0tau.dat");
        for (size_t tt = 0; tt < NTau_; tt++)
        {
            fout << double(tt) / double(NTau_) * beta_ << " " << data_(tt) << std::endl;
        }
    }

  private:
    GreenCluster0Mat gfMatCluster_;
    size_t NTau_;
    double beta_;
    SiteVector_t data_;
};
} // namespace GreenTau
