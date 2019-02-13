

#include "Utilities.hpp"
#include "Model.hpp"

namespace SelfCon
{

using Model_t = Models::Model_2D;

class SelfConsistency
{

  public:
    const size_t Nc = 1;

    SelfConsistency(const Json &jj, const ClusterMatrixCD_t &greenImpurity, const ClusterMatrixCD_t &hybridization) : model_(jj),
                                                                                                                      greenImpurity_(greenImpurity),
                                                                                                                      hybridization_(hybridization),
                                                                                                                      selfEnergy_(),
                                                                                                                      hybNext_()
    {

        const size_t NMat = greenImpurity_.n_rows;

        selfEnergy_.resize(NMat, Nc);
        for (size_t nn = 0; nn < NMat; nn++)
        {
            const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_(nn, 0) = -1.0 / greenImpurity_(nn, 0) + zz - hybridization_(nn, 0);
        }

        std::cout << "In Selfonsistency constructor " << std::endl;
        Save("self", selfEnergy_);
        std::cout << "In Selfonsistency constructor, after save selfenery " << std::endl;
    }

    void DoSCGrid()
    {

        std::cout << "In Selfonsistency DOSC " << std::endl;
        const size_t NMat = greenImpurity_.n_rows;
        const size_t NKPTS = 200;
        const double NKPTS_D = 200.0;

        ClusterMatrixCD_t gImpUpPrime(NMat, Nc);
        gImpUpPrime.zeros();
        ClusterMatrixCD_t hybridizationNext(NMat, Nc);

        for (size_t nn = 0; nn < NMat; nn++)
        {
            const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            for (size_t kxindex = 0; kxindex < NKPTS; kxindex++)
            {
                const double kx = -M_PI + 2.0 * M_PI * static_cast<double>(kxindex) / NKPTS_D;
                for (size_t kyindex = 0; kyindex < NKPTS; kyindex++)
                {
                    const double ky = -M_PI + 2.0 * M_PI * static_cast<double>(kyindex) / NKPTS_D;
                    const double eps0k = -2.0 * (std::cos(kx) + std::cos(ky));
                    gImpUpPrime(nn, 0) += 1.0 / (NKPTS_D * NKPTS_D) * 1.0 / (zz - eps0k - selfEnergy_(nn, 0));
                }
            }
            hybridizationNext(nn, 0) = -1.0 / gImpUpPrime(nn, 0) - selfEnergy_(nn, 0) + zz;
        }

        Save("hybNext", hybridizationNext);
        hybNext_ = hybridizationNext;
    }

    void Save(const std::string &fname, const ClusterMatrixCD_t &green)
    {

        const size_t NMat = green.n_rows;

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        for (size_t nn = 0; nn < NMat; nn++)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / model_.beta();
            fout << iwn << " ";

            fout << green(nn, 0).real()
                 << " "
                 << green(nn, 0).imag()
                 << " ";
            fout << "\n";
        }

        fout.close();

        green.save(fname + std::string(".arma"), arma::arma_ascii);
    }

  private:
    Model_t model_;
    ClusterMatrixCD_t greenImpurity_;
    ClusterMatrixCD_t hybridization_;
    ClusterMatrixCD_t selfEnergy_;
    ClusterMatrixCD_t hybNext_;
};

} // namespace SelfCon