#pragma once
#include "Utilities.hpp"
#include "GreenMat.hpp"

namespace Fourier
{

double MatToTau(const SiteVectorCD_t &greenMat, const double &tau, const double &beta) //Only for a "scalar green function, not a cluster green"
{
    double greenTau = 0.0;
    for (size_t n = 0; n < greenMat.size(); n++)
    {
        const double wn = (2.0 * n + 1.0) * M_PI / beta;
        greenTau += cos(wn * tau) * greenMat(n).real() + sin(wn * tau) * greenMat(n).imag();
    }

    return (2.0 * greenTau / beta);
}

double MatToTauAnalytic(SiteVectorCD_t greenMat, const double &tau, const double &beta, const double &fm, const double &sm, const double &tm)
{

    double result = 0.0;

    //result+= les moments en tau calculés analytiquement
    result += -0.5 * fm;
    result += (tau / 2.0 - beta / 4.0) * sm;
    result += -1.0 / 4.0 * (tau * (tau - beta)) * tm;

    //On transforme la greenMat moins ses moments
    for (size_t n = 0; n < greenMat.n_elem; n++)
    {
        cd_t wn(0.0, (2 * n + 1.0) * M_PI / beta);
        greenMat(n) -= fm / wn + sm / (wn * wn) + tm / (wn * wn * wn);
    }

    result += MatToTau(greenMat, tau, beta);

    return result;
}

} // namespace Fourier