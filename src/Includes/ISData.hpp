#pragma once

#include "Utilities.hpp"
#include "GreenTau.hpp"

namespace Markov
{

namespace Obs
{
using Vertex = Utilities::Vertex;

class ISDataCTINT
{
    using GreenTau_t = GreenTau::GreenCluster0Tau;

  public:
    ISDataCTINT(const double &beta, const GreenTau_t &gf0TauUp) : beta_(beta),
                                                                  green0CachedUp_(gf0TauUp),
                                                                  sign_(1),
                                                                  Mup_(new ClusterMatrix_t()),
                                                                  Mdown_(new ClusterMatrix_t()),
                                                                  vertices_()
    {
    }

    //private:
    const double beta_;
    GreenTau_t green0CachedUp_; //what the hell, this should be const. Better, Moved to Model Shit Fuck

    Sign_t sign_;
    std::shared_ptr<ClusterMatrix_t> Mup_; //The  inverse of the Matrix D defined in review Gull. See Patrick Semon thesis.
    std::shared_ptr<ClusterMatrix_t> Mdown_;
    std::vector<Vertex> vertices_;
};
} // namespace Obs
} // namespace Markov