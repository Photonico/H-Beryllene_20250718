#include "BasisFunctions.hpp"
#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "DescriptorSHS2.hpp"
#include "Frame.hpp"
#include "IoHandlerML_FF.hpp"
#include "Kernel.hpp"
#include "KernelPolynomial.hpp"
#include "MlMPI.hpp"
#include "ParallelEnvironemt.hpp"
#include "Predictor.hpp"
#include "ShmemArray.hpp"
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
    std::string poscar_name = argv[1];
    std::string forceFieldName = argv[2];

    IoHandlerML_FF forceFieldData(forceFieldName);
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

    Structure                               struc;
    std::shared_ptr<NearestNeighborNSquare> neighborLists =
        std::make_shared<NearestNeighborNSquare>(forceFieldData["SHS2-2-body-cutoff"].cget<Real>(),
                                                 true,
                                                 false);
    struc.readPoscar(poscar_name);
    neighborLists->computeNearestNeighborsDirectCoordinates(struc);
    std::shared_ptr<DescriptorSHS2> descriptors =
        std::make_shared<DescriptorSHS2>(forceFieldData["SHS2-2-body-weight"].cget<Real>(),
                                         forceFieldData["SHS2-2-body-is-normalized"].cget<bool>(),
                                         nullptr,
                                         nullptr,
                                         nullptr,
                                         nullptr
                                         //,ExecutionPolicy::gpuStdLib );
        );
    std::cout << "START " << std::endl;

    /*******************************************************************************************
     * classical approach
     *******************************************************************************************/
    //    descriptors->updatePairCoefficients( neighborLists, basisFunctions[ "2-body" ] );
    //    VASPML_PROFILING_START( "Compute descriptors Classical Direct" );
    //    for ( Size_t i = 0; i < 500; i++ )
    //    {
    //        descriptors->updatePairCoefficients( neighborLists, basisFunctions[ "2-body" ] );
    //        descriptors->computeVaspCoefficientsFromPairCoefficients();
    //    }
    //    VASPML_PROFILING_STOP( "Compute descriptors Classical Direct" );
    /*******************************************************************************************
    *******************************************************************************************/

    //    descriptors->set_basisFunctions( basisFunctions[ "2-body" ] );
    //    VASPML_PROFILING_START( "Compute descriptors Direct" );
    //    for ( Size_t i = 0; i < 500; i++ )
    //    descriptors->updateVaspCoefficients( neighborLists );
    //    VASPML_PROFILING_STOP( "Compute descriptors Direct" );

    //    const auto& desc = descriptors -> get_descriptor();
    //
    //    auto file = file_io::openFileO( "Descriptor.dat" );
    //    for ( const auto& x : desc ){
    //        for ( const auto& y : x ){
    //           file << y << std::endl;
    //       }
    //    }

    std::shared_ptr<Frame> singleStructure = std::make_shared<Frame>();
    singleStructure->init(forceFieldData, *basisFunctions["3-body"], ExecutionPolicy::gpuStdLib);
    //    singleStructure->init( forceFieldData, *basisFunctions["3-body"] );
    singleStructure->update(poscar_name, basisFunctions);
    //    singleStructure->set_basisFunctions( basisFunctions );

    //    VASPML_PROFILING_START( "Compute descriptors" );
    //    for ( Size_t i = 0; i < 100; i++ )
    //    singleStructure->update();
    //    VASPML_PROFILING_STOP( "Compute descriptors" );

    std::shared_ptr<KernelPolynomial> polyKernel =
        std::make_shared<KernelPolynomial>(forceFieldData,
                                           nullptr,
                                           nullptr,
                                           nullptr,
                                           ExecutionPolicy::gpuStdLib);
    VASPML_PROFILING_START("update");
    polyKernel->updatePolynomialKernel(*singleStructure);
    VASPML_PROFILING_STOP("update");

    std::cout << "PREDICTOR CREATION" << std::endl;
    std::shared_ptr<Predictor> predictor = std::make_shared<Predictor>(forceFieldData,
                                                                       nullptr,
                                                                       nullptr,
                                                                       nullptr,
                                                                       nullptr,
                                                                       ExecutionPolicy::gpuStdLib);
    predictor->update(*polyKernel);
    std::cout << "PREDICTOR UPDATE DONE" << std::endl;

    Real volume = polyKernel->get_descriptorCollection()
                      .getDescriptor("SHS2-2-body")
                      .get_neighborList_ptr()
                      ->get_latticeVolume();
    VASPML_PROFILING_START("Atomic forces");
    predictor->compute_atomicForces(polyKernel->get_descriptorCollection());
    VASPML_PROFILING_STOP("Atomic forces");
    predictor->compute_stressTensor(polyKernel->get_descriptorCollection(), volume);

    std::cout << str("ENERGY %24.16E\n", predictor->get_totalEnergy() * constants::EUNIT);
    Real pressure = 0.0;
    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << "STRESS ";
        for (std::size_t j = 0; j < 3; j++)
        {
            std::cout << str(" %24.16E", predictor->get_totalStressTensor(i, j) * constants::SUNIT);
        }
        pressure += predictor->get_totalStressTensor(i, i);
        std::cout << "\n";
    }
    pressure *= constants::SUNIT / 3.0;
    std::cout << str("PRESSURE %24.16E\n", pressure);

    for (std::size_t atom = 0; atom < polyKernel->get_descriptorCollection()
                                          .getDescriptor("SHS2-2-body")
                                          .get_neighborList_ptr()
                                          ->get_nAtoms();
         atom++)
    {
        auto& [x, y, z] = predictor->get_atomicForces(atom);
        std::cout << str("FORCE %-6zu %24.16E %24.16E %24.16E\n",
                         atom,
                         x * constants::FUNIT,
                         y * constants::FUNIT,
                         z * constants::FUNIT);
    }
    VASPML_PROFILING_WRITE();

    return 0;
}
