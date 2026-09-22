#include "Record.hpp"
#include "io.hpp"

#include <fstream>
#include <iostream>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    std::ifstream fileIn;
    fileIn.open("ML_AB", std::ifstream::in);
    Record dataset;
    io::readMlab(dataset, fileIn);
    fileIn.close();

    data::MemoryCounter f;
    traverse(dataset, f);

    //std::cout << data::printRecordByteSizes();
    std::cout << data::printMemoryUsage(f.memInfo);

    std::ofstream fileOut;
    fileOut.open("ML_AB.json");
    io::writeJson(dataset, fileOut);
    fileOut.close();

    return 0;
}
