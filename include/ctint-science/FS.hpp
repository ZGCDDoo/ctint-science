#pragma once

#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <stdlib.h>
#include "json.hpp"

namespace IO
{
namespace FS
{
void WriteToFile(const size_t &iter, const double &value, const std::string &name)
{
    std::string fname = name + std::string(".dat");
    std::ofstream fout(fname, std::ios_base::out | std::ios_base::app);
    fout << iter << " " << value << std::endl;
    fout.close();
}

void PrepareNextIter(const std::string paramsName, const size_t &iter)
{

    using boost::filesystem::copy_file;
    copy_file("hybNext.arma", std::string("hyb") + std::to_string(iter + 1) + std::string(".arma"));
    copy_file("hybNext.dat", std::string("hyb") + std::to_string(iter + 1) + std::string(".dat"));
    copy_file("self.dat", std::string("self") + std::to_string(iter) + std::string(".dat"));

    copy_file("greenUp.dat", std::string("greenUp") + std::to_string(iter) + std::string(".dat"));

    std::string fname = paramsName + std::to_string(iter) + std::string(".json");
    std::ifstream fin(fname);
    Json params;
    Json results;
    fin >> params;
    fin.close();
    fin.open("Obs.json");
    fin >> results;
    fin.close();

    for (Json::iterator it = results.begin(); it != results.end(); ++it)
    {
        WriteToFile(iter, it.value(), it.key());
    }

    params["SEED"] = rand() / 2;
    params["HybFile"] = std::string("hyb") + std::to_string(iter + 1);

    std::ofstream fout(std::string("params") + std::to_string(iter + 1) + std::string(".json"));
    fout << std::setw(4) << params << std::endl;
    fout.close();
}
} // namespace FS
} // namespace IO