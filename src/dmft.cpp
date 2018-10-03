
#include "Includes/MonteCarloNew.hpp"
#include "Includes/SelfConsistency.hpp"
#include "Includes/MarkovChain.hpp"
#include "Includes/FS.hpp"
#include "Includes/GreenMat.hpp"
#include "Includes/GreenTau.hpp"

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
    }

    const std::string paramsName = argv[1];
    const int ITER = atoi(argv[2]);
    std::string fname_params = paramsName + std::to_string(ITER) + std::string(".json");
    std::ifstream fin(fname_params);
    Json jj;
    fin >> jj;

    //init a model, to make sure all the files are present and that not all proc write to the same files

    std::cout << "Iter = " << ITER << std::endl;
    size_t seed = jj["SEED"];

    using Markov_t = Markov::MarkovChain;
    using Model_t = Models::Model_2D;
    MC::MonteCarlo monteCarloMachine(std::make_shared<Markov_t>(jj, seed), jj);
    monteCarloMachine.RunMonteCarlo();
    Model_t model(jj);
    ClusterCubeCD_t greenImpurity;
    assert(greenImpurity.load("greenUp.dat"));

    SelfCon::SelfConsistency selfcon(jj, greenImpurity, model.hybridizationMatUp());
    selfcon.DoSCGrid();
    IO::FS::PrepareNextIter(paramsName, ITER);

    return EXIT_SUCCESS;
}
