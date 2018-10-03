#pragma once

#include "Utilities.hpp"
#include "Fourier.hpp"
#include "GreenTau.hpp"
#include "Observables.hpp"
#include "ISData.hpp"
#include "Model.hpp"

//#define DEBUG_TEST

namespace Markov
{

using Vertex = Utilities::Vertex;
using Utilities::Dot;
using Model_t = Models::Model_2D;

class MarkovChain
{

    using GreenTau_t = GreenTau::GreenCluster0Tau;

  public:
    const size_t Nc = 1;
    const double PROBINSERT = 0.5;
    const double PROBREMOVE = 1.0 - PROBINSERT;

    MarkovChain(const Json &jj, const size_t &seed) : model_(jj),
                                                      rng_(seed),
                                                      urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0)),
                                                      NMat_(model_.greenCluster0MatUp().data().n_rows),
                                                      NTau_(static_cast<size_t>(jj["NTAU"])),
                                                      totalNMeas_(static_cast<size_t>(jj["TOTALNMEAS"])),
                                                      datactint_(
                                                          new Obs::ISDataCTINT(
                                                              static_cast<double>(jj["beta"]),
                                                              GreenTau_t(model_.greenCluster0MatUp(), NTau_))),
                                                      obs_(model_, datactint_, totalNMeas_)
    {
    }

    ~MarkovChain(){};

    void DoStep()
    {
        urng_() < PROBINSERT ? InsertVertex() : RemoveVertex();
    }

    void InsertVertex()
    {
        Vertex vertex = Vertex(datactint_->beta_ * urng_(), urng_() < 0.5 ? Up : Down);

        double sUp = GetGreenTau0Up(vertex, vertex) - model_.auxUp(vertex);
        double sDown = GetGreenTau0Up(vertex, vertex) - model_.auxDown(vertex);

        if (datactint_->vertices_.size())
        {

            //The following should be in a seperate method.
            //TO DO.
            const size_t kkold = datactint_->vertices_.size();
            const size_t kknew = datactint_->vertices_.size() + 1;

            //by convention
            newLastColUp_.set_size(kkold);
            newLastRowUp_.set_size(kkold);
            newLastColDown_.set_size(kkold);
            newLastRowDown_.set_size(kkold);

            double sTildeUpI = sUp;
            double sTildeDownI = sDown;

            //Probably put this in a method
            for (size_t i = 0; i < kkold; i++)
            {
                newLastRowUp_(i) = GetGreenTau0Up(vertex, datactint_->vertices_.at(i));
                newLastColUp_(i) = GetGreenTau0Up(datactint_->vertices_.at(i), vertex);

                newLastRowDown_(i) = newLastRowUp_(i);
                newLastColDown_(i) = newLastColUp_(i);
            }

            sTildeUpI -= Dot(newLastRowUp_, (*(datactint_->Mup_)), newLastColUp_);
            sTildeDownI -= Dot(newLastRowDown_, (*(datactint_->Mdown_)), newLastColDown_);

            const double ratio = sTildeUpI * sTildeDownI;
            double probAcc = (-2.0 * Nc * datactint_->beta_ * model_.auxU()) / kknew * ratio;
            probAcc *= PROBREMOVE / PROBINSERT;

            if (urng_() < std::abs(probAcc))
            {
                if (probAcc < .0)
                {
                    datactint_->sign_ *= -1;
                }

                Utilities::BlockRankOneUpgrade((*(datactint_->Mup_)), newLastColUp_, newLastRowUp_, 1.0 / sTildeUpI);
                Utilities::BlockRankOneUpgrade((*(datactint_->Mdown_)), newLastColDown_, newLastRowDown_, 1.0 / sTildeDownI);
                datactint_->vertices_.push_back(vertex);
            }
        }
        else
        {
            double probAcc = -2.0 * Nc * datactint_->beta_ * model_.auxU() * sUp * sDown;
            probAcc *= PROBREMOVE / PROBINSERT;
            if (urng_() < std::abs(probAcc))
            {
                if (probAcc < .0)
                {
                    datactint_->sign_ *= -1;
                }
                (*(datactint_->Mup_)) = ClusterMatrix_t(1, 1);
                (*(datactint_->Mdown_)) = ClusterMatrix_t(1, 1);
                (*(datactint_->Mup_))(0, 0) = 1.0 / sUp;
                (*(datactint_->Mdown_))(0, 0) = 1.0 / sDown;
                datactint_->vertices_.push_back(vertex);
            }
        }
    }

    void RemoveVertex()
    {
        if (datactint_->vertices_.size())
        {
            const size_t pp = static_cast<int>(urng_() * datactint_->vertices_.size());

            double probAcc = static_cast<double>(datactint_->vertices_.size()) / (-2.0 * Nc * datactint_->beta_ * model_.auxU()) * ((*(datactint_->Mup_))(pp, pp) * (*(datactint_->Mdown_))(pp, pp));
            probAcc *= PROBINSERT / PROBREMOVE;
            if (urng_() < std::abs(probAcc))
            {
                if (probAcc < .0)
                {
                    datactint_->sign_ *= -1; //not to sure here, should it not just be sign = -1 ??
                }

                //The update matrices of size k-1 x k-1 with the pp row and col deleted and the last row and col now at index pp

                const size_t LL = datactint_->vertices_.size();
                (*(datactint_->Mup_)).swap_rows(pp, LL - 1);
                (*(datactint_->Mup_)).swap_cols(pp, LL - 1);
                (*(datactint_->Mdown_)).swap_rows(pp, LL - 1);
                (*(datactint_->Mdown_)).swap_cols(pp, LL - 1);
                Utilities::BlockRankOneDowngrade((*(datactint_->Mup_)), dummy_);
                (*(datactint_->Mup_)) = dummy_;
                Utilities::BlockRankOneDowngrade((*(datactint_->Mdown_)), dummy_);
                (*(datactint_->Mdown_)) = dummy_;
                std::iter_swap(datactint_->vertices_.begin() + pp, datactint_->vertices_.begin() + datactint_->vertices_.size() - 1); //swap the last vertex and the vertex pp in vertices.
                                                                                                                                      //to be consistent with the updated Mup and (*(datactint_->Mdown_))
                datactint_->vertices_.pop_back();
            }
        }
    }

    void CleanUpdate()
    {
        size_t kk = datactint_->vertices_.size();
        ClusterMatrix_t Dup(kk, kk);
        ClusterMatrix_t Ddown(kk, kk);

        for (size_t i = 0; i < kk; i++)
        {
            for (size_t j = 0; j < kk; j++)
            {
                Dup(i, j) = GetGreenTau0Up(datactint_->vertices_.at(i), datactint_->vertices_.at(j));
                Ddown(i, j) = GetGreenTau0Up(datactint_->vertices_.at(i), datactint_->vertices_.at(j));
                if (i == j)
                {
                    Dup(i, i) -= model_.auxUp(datactint_->vertices_.at(i));
                    Ddown(i, i) -= model_.auxDown(datactint_->vertices_.at(i));
                }
            }
        }

        (*(datactint_->Mup_)) = Dup.i();
        (*(datactint_->Mdown_)) = Ddown.i();

        Dup.clear();
        Ddown.clear();
    }

    double GetGreenTau0Up(const Vertex &vertexI, const Vertex &vertexJ) const
    {
        const double delta = 1e-12;
        Tau_t tauDiff = vertexI.tau() - (vertexJ.tau() + delta);
        return (datactint_->green0CachedUp_(tauDiff));
    }

    void Measure()
    {
        // std::cout << "Measuring, sign, k " << datactint_->sign_ << ",  " << datactint_->vertices_.size() << std::endl;
        obs_.Measure();
    }

    void Save()
    {
        obs_.Save();
    }

  private:
    //attributes
    Model_t model_;
    Utilities::EngineType_t rng_;
    Utilities::UniformRng_t urng_;

    const size_t NMat_;
    const size_t NTau_;
    const size_t totalNMeas_;

    std::shared_ptr<Obs::ISDataCTINT> datactint_;
    Obs::Observables obs_;

    ClusterMatrix_t dummy_;
    SiteVector_t newLastColUp_;
    SiteVector_t newLastRowUp_;
    SiteVector_t newLastColDown_;
    SiteVector_t newLastRowDown_;
};
} // namespace Markov
