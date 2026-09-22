#include "BasisFunctions.hpp"
#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "Frame.hpp"
#include "IoHandlerML_FF.hpp"
#include "Kernel.hpp"
#include "KernelPolynomial.hpp"
#include "MlMPI.hpp"
#include "Predictor.hpp"
#include "ShmemArray.hpp"
#include "SpillingFactor.hpp"
#include "Structure.hpp"
#include "Timer.hpp"
#include "constants.hpp"
#include "cutoff.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using namespace vaspml;

int main(int argc, char** argv)
{

    using CT = vaspml::math::CutoffType;

    if (argc < 3)
    {
        std::cout << str("USAGE: %s <POSCAR> <ML_FF>\n", argv[0]);
        throw std::runtime_error("ERROR: Wrong number of arguments.");
    }
    std::shared_ptr<MlMPI> mpiMain;
    mpiMain = std::make_shared<MlMPI>();
    mpiMain->make_WorldComm();
    std::string poscar_name = argv[1];
    std::string forceFieldName = argv[2];

    IoHandlerML_FF forceFieldData(forceFieldName, mpiMain);
    forceFieldData.mainReader();
    forceFieldData.convertUnitsToAtomicUnits();

    //    std::cout << "Setting up basis functions " << std::endl;
    std::map<std::string, std::shared_ptr<BasisFunctions>> basisFunctions;

    basisFunctions["2-body"] = std::make_shared<BasisFunctionsRadialSpline>(
        static_cast<CT>(std::get<Int>(forceFieldData["SHS2-2-body-cutoff-type"])),
        std::get<Real>(forceFieldData["SHS2-2-body-cutoff"]),
        std::get<Real>(forceFieldData["SHS2-2-body-smearing-param"]),
        std::get<Int>(forceFieldData["SHS2-2-body-number-radial-grid-points"]),
        0,
        std::get<Int>(forceFieldData["SHS2-2-body-max-radial-funcs"]),
        BasisFunctionType::bodyOrder2);
    basisFunctions["3-body"] = std::make_shared<BasisFunctionsAngular>(
        static_cast<CT>(std::get<Int>(forceFieldData["SHS3-3-body-cutoff-type"])),
        std::get<Real>(forceFieldData["SHS3-3-body-cutoff"]),
        std::get<Real>(forceFieldData["SHS3-3-body-smearing-param"]),
        std::get<Int>(forceFieldData["SHS3-3-body-number-radial-grid-points"]),
        std::get<Int>(forceFieldData["SHS3-3-body-max-angular-number"]),
        std::get<Int>(forceFieldData["SHS3-3-body-max-radial-funcs"]),
        BasisFunctionType::bodyOrder3);

    VASPML_PROFILING_START("update");
    Frame singleStructure;
    singleStructure.init(forceFieldData, *basisFunctions["3-body"]);
    singleStructure.update(poscar_name, basisFunctions);
    KernelPolynomial polyKernel;
    polyKernel.init(forceFieldData, mpiMain);

    polyKernel.updatePolynomialKernel(singleStructure);

    Predictor predictor;
    predictor.init(forceFieldData, mpiMain);

    predictor.update(polyKernel);

    Real volume = polyKernel.get_descriptorCollection()
                      .getDescriptor("SHS2-2-body")
                      .get_neighborList_ptr()
                      ->get_latticeVolume();
    predictor.compute_atomicForces(polyKernel.get_descriptorCollection());
    predictor.compute_stressTensor(polyKernel.get_descriptorCollection(), volume);
    VASPML_PROFILING_STOP("update");

    std::cout << str("ENERGY %24.16E\n", predictor.get_totalEnergy() * constants::EUNIT);
    Real pressure = 0.0;
    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << "STRESS ";
        for (std::size_t j = 0; j < 3; j++)
        {
            std::cout << str(" %24.16E", predictor.get_totalStressTensor(i, j) * constants::SUNIT);
        }
        pressure += predictor.get_totalStressTensor(i, i);
        std::cout << "\n";
    }
    pressure *= constants::SUNIT / 3.0;
    std::cout << str("PRESSURE %24.16E\n", pressure);

    for (std::size_t atom = 0; atom < polyKernel.get_descriptorCollection()
                                          .getDescriptor("SHS2-2-body")
                                          .get_neighborList_ptr()
                                          ->get_nAtoms();
         atom++)
    {
        auto& [x, y, z] = predictor.get_atomicForces(atom);
        std::cout << str("FORCE %-6zu %24.16E %24.16E %24.16E\n",
                         atom,
                         x * constants::FUNIT,
                         y * constants::FUNIT,
                         z * constants::FUNIT);
    }

    if (mpiMain->get_rank() == 0) { VASPML_PROFILING_WRITE(); }

    {
        SpillingFactor spilli(forceFieldData, mpiMain);
        spilli.computeSpillingFactor(*polyKernel.get_kernelMatrix(),
                                     *polyKernel.get_nAtomsType(),
                                     *singleStructure.get_typeMap());
        spilli.writeToScreen();
    }

    MPI_Finalize();
    return 0;
}
