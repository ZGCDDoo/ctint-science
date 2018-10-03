#pragma once
#include "MarkovChain.hpp"

namespace MC
{
using Clock_t = std::chrono::high_resolution_clock;
using MarkovChain_t = Markov::MarkovChain;

class MonteCarlo
{
  public:
    MonteCarlo(const std::shared_ptr<MarkovChain_t> &markovchainPtr, const Json &jj) : markovchainPtr_(markovchainPtr),
                                                                                       updatesMeas_(jj["UPDATESMEAS"]),
                                                                                       thermalization_(jj["THERMALIZATION"]),
                                                                                       totalNMeas_(jj["TOTALNMEAS"]),
                                                                                       cleanUpdate_(jj["CLEANUPDATE"]),
                                                                                       NMeas_(0),
                                                                                       NCleanUpdates_(0),
                                                                                       updatesProposed_(0)
    {
    }

    ~MonteCarlo(){};

    void RunMonteCarlo()
    {
        auto startTime = Clock_t::now();
        std::cout << std::string("Start Thermalization, ") << std::endl;

        while (NMeas_ < thermalization_)
        {
            updatesProposed_++;
            markovchainPtr_->DoStep();
            if (updatesProposed_ % updatesMeas_ == 0)
            {
                NMeas_++;
            }
        }

        auto endTime = Clock_t::now();
        std::cout << "End Thermalization, Time " + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()) << std::endl;

        NMeas_ = 0;
        updatesProposed_ = 0;
        startTime = Clock_t::now();
        std::cout << std::string("Start Measurements, ") << std::endl;

        while (NMeas_ < totalNMeas_)
        {
            markovchainPtr_->DoStep(); //One simple sweep
            updatesProposed_++;

            if (updatesProposed_ % updatesMeas_ == 0)
            {
                markovchainPtr_->Measure();
                if (NMeas_ % cleanUpdate_ == 0 && NMeas_ != 0)
                {
                    markovchainPtr_->CleanUpdate();
                    NCleanUpdates_++;
                }
                NMeas_++;
            }
        }

        endTime = Clock_t::now();
        std::cout << std::string("End Measurements, Time ") + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()) << std::endl;
        markovchainPtr_->Save();
        return;
    }

    //Getters
    size_t NMeas() const { return NMeas_; };
    size_t NCleanUpdates() const { return NCleanUpdates_; };
    size_t updatesProposed() const { return updatesProposed_; };

  private:
    //attributes
    const std::shared_ptr<MarkovChain_t> markovchainPtr_;
    const size_t updatesMeas_;
    const size_t thermalization_;
    const size_t totalNMeas_;
    const size_t cleanUpdate_;

    size_t NMeas_;
    size_t NCleanUpdates_;
    size_t updatesProposed_;
};
} // namespace MC
