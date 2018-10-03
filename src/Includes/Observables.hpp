#pragma once

#include "Utilities.hpp"
#include "Fourier.hpp"
#include "GreenTau.hpp"
#include "ISData.hpp"
#include "Model.hpp"

namespace Markov
{
namespace Obs
{
using namespace Utilities;
using Model_t = Models::Model_2D;

class Observables
{

      public:
        Observables(const Model_t &model, const std::shared_ptr<ISDataCTINT> &datactint,
                    const size_t &totalNMeas) : model_(model),
                                                datactint_(datactint),
                                                totalNMeas_(totalNMeas),
                                                NMat_(model_.greenCluster0MatUp().data().n_rows),
                                                signMeas_(0.0),
                                                expOrder_(0.0),
                                                greenMatsubaraUp_(NMat_, 1),
                                                NMeas_(0)

        {
                greenMatsubaraUp_.zeros();

                Maverage_.resize(NMat_, 0.0);
                const size_t Nc = 1;
        }

        //Getters
        double signMeas() const { return signMeas_; };
        double expOrder() const { return expOrder_; };
        ClusterMatrixCD_t greenMatsubaraUp() const { return greenMatsubaraUp_; };

        void Measure()
        {
                NMeas_++;
                Maveraged_ = 0.5 * ((*(datactint_->Mup_)) + (*(datactint_->Mdown_)));
                signMeas_ += static_cast<double>(datactint_->sign_) / static_cast<double>(totalNMeas_);
                double tmp = static_cast<double>(datactint_->vertices_.size()) * static_cast<double>(datactint_->sign_) / static_cast<double>(totalNMeas_);
                expOrder_ += tmp;

                MeasureGreen();
        }

        void MeasureGreen() //Part not needed to be sampled are added to the greens in Save()
        {
                const double epsilon = 1e-12;
                double fact = datactint_->sign_ / static_cast<double>(totalNMeas_);

                for (size_t nn = 0; nn < NMat_; nn++)
                {
                        const cd_t iwn = cd_t(0.0, (2.0 * nn + 1.0) * M_PI / model_.beta());
                        for (size_t kk = 0; kk < datactint_->vertices_.size(); kk++)
                        {
                                for (size_t mm = 0; mm < datactint_->vertices_.size(); mm++)
                                {
                                        Tau_t tau = datactint_->vertices_.at(kk).tau() + epsilon - datactint_->vertices_.at(mm).tau();

                                        const cd_t expkkmm = std::exp(iwn * tau);
                                        Maverage_.at(nn) += fact * expkkmm * Maveraged_(kk, mm);
                                }
                        }
                }
        }

        void Save()
        {
                const double expOrderResult = expOrder_ / signMeas_;

                Json jjout;

                //Je vais aussi ecrire les resultats dans un fichier pour des tests. Penser a mettre cela dans un
                //map<std::string, double> prochainement.
                jjout["sign"] = signMeas_;
                jjout["k"] = expOrderResult;

                const std::string fnameObs = std::string("Obs") + std::string(".json");
                std::ofstream fout(fnameObs);
                fout << std::setw(4) << jjout << std::endl;
                fout.close();

                //The greens
                FinalizeGreens();
                std::ofstream fup;
                std::string fnameGreenUp = "greenUp.dat";
                fup.open(fnameGreenUp, std::ios::out);

                for (size_t nn = 0; nn < greenMatsubaraUp_.n_rows; nn++)
                {
                        double iwn = (2.0 * nn + 1.0) * M_PI / datactint_->beta_;
                        fup << iwn << " ";
                        for (size_t ii = 0; ii < greenMatsubaraUp_.n_cols; ii++)
                        {
                                fup << greenMatsubaraUp_(nn, ii).real() << " " << greenMatsubaraUp_(nn, ii).imag() << " ";
                        }
                        fup << "\n";
                }
                fup.close();

                fnameGreenUp = "greenUp.arma";
                greenMatsubaraUp_.save(fnameGreenUp, arma::raw_ascii);
        }

        void FinalizeGreens()
        {

                //First use the fact that some sites are equivalent for getting  smaller MsigmaAverages matrices. Then do the dot product.

                for (size_t nn = 0; nn < NMat_; nn++)
                {
                        const cd_t tmp = model_.greenCluster0MatUp().data()(nn, 0);
                        greenMatsubaraUp_(nn, 0) = tmp - tmp * tmp * Maverage_.at(nn) / (model_.beta() * signMeas_);
                }
        }

      private:
        Model_t model_;
        std::shared_ptr<ISDataCTINT> datactint_;
        size_t totalNMeas_;
        size_t NMat_;
        double signMeas_;
        double expOrder_;
        ClusterMatrixCD_t greenMatsubaraUp_;

        std::vector<cd_t> Maverage_;
        ClusterMatrix_t Maveraged_;
        size_t NMeas_;
}; // namespace Obs

} // namespace Obs
} // namespace Markov
