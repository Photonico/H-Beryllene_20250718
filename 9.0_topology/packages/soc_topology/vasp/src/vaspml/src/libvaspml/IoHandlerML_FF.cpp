#include "IoHandlerML_FF.hpp"

#include "SemanticVersion.hpp"
#include "constants.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace vaspml;

// constructor
IoHandlerML_FF::IoHandlerML_FF(const std::string& fname, const std::shared_ptr<MlMPI>& mpiIn) :
    fileName(fname),
    mpiMain(mpiIn)
{
#ifdef use_shmem
    bool readAllCoresOn = false;
#else
    bool readAllCoresOn = true;
#endif
    rankReadsOn = false;
    if (mpiMain == nullptr) rankReadsOn = true;
    else if (mpiMain->get_rank() == 0 || readAllCoresOn) rankReadsOn = true;
}

void IoHandlerML_FF::init(const std::string& fname, const std::shared_ptr<MlMPI>& mpiIn)
{
    fileName = fname;
    mpiMain = mpiIn;

#ifdef use_shmem
    bool readAllCoresOn = false;
#else
    bool readAllCoresOn = true;
#endif
    rankReadsOn = false;
    if (mpiMain == nullptr) rankReadsOn = true;
    else if (mpiMain->get_rank() == 0 || readAllCoresOn) rankReadsOn = true;
}

void IoHandlerML_FF::mainReader()
{
    // open ML_FF file
    input_stream = file_io::openFileI(fileName, std::ifstream::binary);

    // read header
    read_header();

    // read version number
    Int major = 0;
    Int minor = 0;
    Int patch = 0;
    read_values::io_scalar(input_stream, major);
    read_values::io_scalar(input_stream, minor);
    read_values::io_scalar(input_stream, patch);
    versionIn = SemanticVersion(major, minor, patch);

    // check version number of header and binary
    check_version_number_header_and_binary();

    if (versionIn > versionMax)
    {
        global_scope::tutor.error("The given force field file is too new " + versionIn.toString()
                                  + "), this vaspml distribution only supports"
                                    "files up to version "
                                  + versionMax.toString() + ".");
        //throw std::runtime_error(
        //      "The given force field file is too new " + versionIn.toString()
        //       + "), this vaspml distribution only supports"
        //       "files up to version "
        //        + versionMax.toString() + "."
        //      );
    }

    // Read ML_DESC_TYPE - descriptor type
    data["ML_DESC_TYPE"] = (Int)(0);
    if (versionIn >= SemanticVersion(0, 2, 1))
    {
        read_values::io_scalar(input_stream, data["ML_DESC_TYPE"].get<Int>());
    }
    const Int descriptorType = data["ML_DESC_TYPE"].cget<Int>();

    // Read number of types
    data["number-types"] = (Int)(0);
    read_values::io_scalar(input_stream, data["number-types"].get<Int>());
    const Int ntypes = data["number-types"].cget<Int>();

    // Read element type names into array of strings
    data["types"] = std::make_shared<Vec1String>();
    read_values::io_array_binary(input_stream, data["types"].dget<ShVec1String>(), ntypes);
    string_tools::trimVector(data["types"].dget<ShVec1String>());

    // Read reference atom energies (units: eV per default)
    data["ML_EATOM_REF"] = std::make_shared<Vec1Real>();
    read_values::io_array_binary(input_stream, data["ML_EATOM_REF"].dget<ShVec1Real>(), ntypes);

    // Read atomic masses
    data["atomic-masses"] = std::make_shared<Vec1Real>();
    read_values::io_array_binary(input_stream, data["atomic-masses"].dget<ShVec1Real>(), ntypes);

    // Read ML_IWEIGHT - type of weighting
    data["ML_IWEIGHT"] = (Int)(0);
    read_values::io_scalar(input_stream, data["ML_IWEIGHT"].get<Int>());

    // Read ML_WTOTEN - weight for energies
    data["ML_WTOTEN"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["ML_WTOTEN"].get<Real>());

    // Read ML_WTIFOR - weight for forces
    data["ML_WTIFOR"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["ML_WTIFOR"].get<Real>());

    // Read ML_WTSIF - weight for forces
    data["ML_WTSIF"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["ML_WTSIF"].get<Real>());

    // Read root mean square error in energy (eV atom^-1)
    data["rmse-total-energy"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["rmse-total-energy"].get<Real>());

    // Read root mean square error in forces (eV Angst^-1)
    data["rmse-forces"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["rmse-forces"].get<Real>());

    // Read root mean square error in stresses (kB)
    data["rmse-stresses"] = (Real)(0.0);
    read_values::io_scalar(input_stream, data["rmse-stresses"].get<Real>());

    // Read many body switch - always true
    data["many-body"] = false;
    read_values::io_scalar(input_stream, data["many-body"].get<bool>());

    // Although this is always true we need to check here
    if (data["many-body"].cget<bool>())
    {
        // Read evil self interaction on/off switch
        data["ML_LSIC"] = false;
        read_values::io_scalar(input_stream, data["ML_LSIC"].get<bool>());

        // Read ML_LSUPERVEC - evil super vector switch
        data["ML_LSUPERVEC"] = false;
        read_values::io_scalar(input_stream, data["ML_LSUPERVEC"].get<bool>());

        // Read ML_W1 - weight for two body descriptor
        data["ML_W1"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_W1"].get<Real>());

        // Read ML_W2 - weight for three body descriptor
        // (this doesn't exist in this version in fortran version)
        data["ML_W2"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_W2"].get<Real>());

        // Read ML_ICUT1 - type of cutoff function for two body descriptor
        data["ML_ICUT1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_ICUT1"].get<Int>());

        // Read ML_ICUT2 - type of cutoff function for three body descriptor
        data["ML_ICUT2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_ICUT2"].get<Int>());

        // Read ML_RCUT1 - cutoff radius for two body descriptor (Angstrom)
        data["ML_RCUT1"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_RCUT1"].get<Real>());

        // Read ML_RCUT2 - cutoff radius for three body descriptor (Angstrom)
        data["ML_RCUT2"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_RCUT2"].get<Real>());

        // Read ML_IBROAD1 - broadening type for two body descriptor
        data["ML_IBROAD1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_IBROAD1"].get<Int>());

        // Read ML_IBROAD2 - broadening type for three body descriptor
        data["ML_IBROAD2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_IBROAD2"].get<Int>());

        // Read ML_SION1 - broadening value for two body descriptor
        data["ML_SION1"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_SION1"].get<Real>());

        // Read ML_SION2 - broadening value for three body descriptor
        data["ML_SION2"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["ML_SION2"].get<Real>());

        // Read ML_MRB1 - Maximum number of radial functions for two body descriptor
        data["ML_MRB1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_MRB1"].get<Int>());

        // Read ML_MRB2 - Maximum number of radial functions for three body descriptor
        data["ML_MRB2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_MRB2"].get<Int>());

        if (versionIn >= SemanticVersion(0, 2, 2))
        {
            // Read ML_BASIS_TYPE2 - type of basis function
            data["ML_BASIS_TYPE2"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_BASIS_TYPE2"].get<Int>());
            const Int basisType = data["ML_BASIS_TYPE2"].cget<Int>();

            // Read ML_MRB2_MAX - maximum number of summation for basis_type 4 and 5
            if (basisType == 4 || basisType == 5)
            {
                data["ML_MRB2_MAX"] = (Int)(0);
                read_values::io_scalar(input_stream, data["ML_MRB2_MAX"].get<Int>());
            }

            // Reading variables for dual descriptor
            if (descriptorType == 2 || descriptorType == 4)
            {
                // Read ML_MRB2_RED - Maximum number of radial functions for dual descriptor
                data["ML_MRB2_RED"] = (Int)(0);
                read_values::io_scalar(input_stream, data["ML_MRB2_RED"].get<Int>());

                // Read ML_BASIS_TYPE2_RED - type of basis function for dual descriptor
                data["ML_BASIS_TYPE2_RED"] = (Int)(0);
                read_values::io_scalar<Int>(input_stream, data["ML_BASIS_TYPE2_RED"].get<Int>());
                const Int basisTypeRed = data["ML_BASIS_TYPE2_RED"].cget<Int>();

                // Read ML_MRB2_MAX_RED - maximum number of summation for basis_type 4 and 5
                if (basisTypeRed == 4 || basisTypeRed == 5)
                {
                    data["ML_MRB2_MAX_RED"] = (Int)(0);
                    read_values::io_scalar(input_stream, data["ML_MRB2_MAX_RED"].get<Int>());
                }
            }
        }

        // Read ML_NR1 - Number of grid points for radial integration for two body descriptor
        data["ML_NR1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_NR1"].get<Int>());

        // Read ML_NR2 - Number of grid points for radial integration for three body descriptor
        data["ML_NR2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_NR2"].get<Int>());

        // Read ML_MSPL1 - Number of grid points for spline interpolation for two body descriptor
        data["ML_MSPL1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_MSPL1"].get<Int>());

        // Read ML_MSPL2 - Number of grid points for spline interpolation for three body descriptor
        data["ML_MSPL2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_MSPL2"].get<Int>());

        // Read ML_LMAX1 - Maximum l for spherical harmonics for two body descriptor
        data["ML_LMAX1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_LMAX1"].get<Int>());

        // Read ML_LMAX2 - Maximum l for spherical harmonics for three body descriptor
        data["ML_LMAX2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_LMAX2"].get<Int>());

        // Reading variables for dual descriptor
        if (descriptorType == 2 || descriptorType == 4)
        {
            // Read ML_LMAX2_RED - Maximum l for spherical harmonics for dual descriptor
            data["ML_LMAX2_RED"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_LMAX2_RED"].get<Int>());

            // Read ML_DESC_RATIO_DUAL - Factor of first descriptor within dual descriptor
            data["ML_DESC_RATIO_DUAL"] = (Real)(0.0);
            read_values::io_scalar(input_stream, data["ML_DESC_RATIO_DUAL"].get<Real>());

            // Read factor of second descriptor within dual descriptor
            data["desc_ratio_dual_second"] = (Real)(0.0);
            read_values::io_scalar(input_stream, data["desc_ratio_dual_second"].get<Real>());

            // amplification factor for three body descriptor SHS3
            data["ML_DESC_FACTORE_TESTE"] = (Real)(0.0);
            read_values::io_scalar(input_stream, data["ML_DESC_FACTORE_TESTE"].get<Real>());
        }

        // Read ML_NHYP1 - Exponential coefficient for two body descriptor
        data["ML_NHYP1"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_NHYP1"].get<Int>());

        // Read ML_NHYP2 - Exponential coefficient for three body descriptor
        data["ML_NHYP2"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_NHYP2"].get<Int>());

        // Read ML_LNORM1 - Normalization of descriptor for two body descriptor
        data["ML_LNORM1"] = false;
        read_values::io_scalar(input_stream, data["ML_LNORM1"].get<bool>());

        // Read ML_LNORM2 - Normalization of descriptor for three body descriptor
        data["ML_LNORM2"] = false;
        read_values::io_scalar(input_stream, data["ML_LNORM2"].get<bool>());

        // Read ML_LWINDOW1 - Use of window function for two body descriptor
        data["ML_LWINDOW1"] = false;
        read_values::io_scalar(input_stream, data["ML_LWINDOW1"].get<bool>());

        // Read ML_LWINDOW2 - Use of window function for three body descriptor
        data["ML_LWINDOW2"] = false;
        read_values::io_scalar(input_stream, data["ML_LWINDOW2"].get<bool>());

        // Optional: Read ML_IWINDOW1  - Type of window function for two body descriptor
        if (data["ML_LWINDOW1"].cget<bool>())
        {
            data["ML_IWINDOW1"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_IWINDOW1"].get<Int>());
        }

        // Optional: Read ML_IWINDOW2  - Type of window function for three body descriptor
        if (data["ML_LWINDOW2"].cget<bool>())
        {
            data["ML_IWINDOW2"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_IWINDOW2"].get<Int>());
        }

        // Read ML_LAFILT2 - Angular filtering on/off switch
        data["ML_LAFILT2"] = false;
        read_values::io_scalar(input_stream, data["ML_LAFILT2"].get<bool>());

        // Optional: Only done if angular filtering is on
        if (data["ML_LAFILT2"].cget<bool>())
        {
            // ML_IAFILT2: Angular filtering type
            data["ML_IAFILT2"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_IAFILT2"].get<Int>());

            // ML_AFILT2: Angular filtering parameter for type 2
            if (data["ML_IAFILT2"].cget<Int>() == 2)
            {
                data["ML_AFILT2"] = (Real)(0.0);
                read_values::io_scalar(input_stream, data["ML_AFILT2"].get<Real>());
            }
        }

        // Read ML_LMETRIC1 - Metric function for two body descriptor
        data["ML_LMETRIC1"] = false;
        read_values::io_scalar(input_stream, data["ML_LMETRIC1"].get<bool>());

        // Read ML_LMETRIC2 - Metric function for three body descriptor
        data["ML_LMETRIC2"] = false;
        read_values::io_scalar(input_stream, data["ML_LMETRIC2"].get<bool>());

        // Optional: Only done if two body metric is true
        if (data["ML_LMETRIC1"].cget<bool>())
        {
            // ML_NMETRIC1: Number of metric functions
            data["ML_NMETRIC1"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_NMETRIC1"].get<Int>());

            // ML_RMETRIC1: Metric factor
            data["ML_RMETRIC1"] = (Real)(0.0);
            read_values::io_scalar(input_stream, data["ML_RMETRIC1"].get<Real>());
        }

        // Optional: Only done if thre body metric is true
        if (data["ML_LMETRIC2"].cget<bool>())
        {
            // ML_NMETRIC2: Number of metric functions
            data["ML_NMETRIC2"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_NMETRIC2"].get<Int>());

            // ML_RMETRIC2: Metric factor
            data["ML_RMETRIC2"] = (Real)(0.0);
            read_values::io_scalar(input_stream, data["ML_RMETRIC2"].get<Real>());
        }

        // Read ML_LVARTRAN1 - Variable transformation for two body descriptor
        data["ML_LVARTRAN1"] = false;
        read_values::io_scalar(input_stream, data["ML_LVARTRAN1"].get<bool>());

        // Read ML_LVARTRAN2 - Variable transformation for three body descriptor
        data["ML_LVARTRAN2"] = false;
        read_values::io_scalar(input_stream, data["ML_LVARTRAN2"].get<bool>());

        // Optional: Only done if two body variable transformation is true
        if (data["ML_LVARTRAN1"].cget<bool>())
        {
            // ML_NVARTRAN1: Number of variable transforms
            data["ML_NVARTRAN1"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_NVARTRAN1"].get<Int>());
        }

        // Optional: Only done if three body variable transformation is true
        if (data["ML_LVARTRAN2"].cget<bool>())
        {
            // ML_NVARTRAN2: Number of variable transforms
            data["ML_NVARTRAN2"] = (Int)(0);
            read_values::io_scalar(input_stream, data["ML_NVARTRAN2"].get<Int>());
        }

        // The following is only done if angular descriptor exists
        if (data["ML_W2"].cget<Real>() > 0.0)
        {
            // Read number of radial functions for each angular channel
            // nrb2
            data["SHS3-3-body-number-radial-basis"] = std::make_shared<Vec1Int>();
            read_values::io_array_binary(input_stream,
                                         data["SHS3-3-body-number-radial-basis"].dget<ShVec1Int>(),
                                         data["ML_LMAX2"].get<Int>() + 1);

            // Reading variables for dual descriptor
            if (descriptorType == 2 || descriptorType == 4)
            {
                // Read number of radial functions for each angular channel
                data["SHS3-3-body-reduced-number-radial-basis"] = std::make_shared<Vec1Int>();
                read_values::io_array_binary(
                    input_stream,
                    data["SHS3-3-body-reduced-number-radial-basis"].dget<ShVec1Int>(),
                    data["ML_LMAX2_RED"].get<Int>() + 1);
            }

            // Make array for number of descriptors per type
            // ndesc_per_type2
            data["SHS3-3-body-number-descriptors-per-type"] = std::make_shared<Vec1Int>();

            // Make a 2d vector and resize its outer dimension to number of element types
            // descriptor_list2
            data["SHS3-3-body-descriptor-list"] = std::make_shared<Vec2Int>();
            data["SHS3-3-body-descriptor-list"].dget<ShVec2Int>().resize(ntypes);

            // Same as above for dual descriptor
            if (descriptorType == 2 || descriptorType == 4)
            {
                data["ndesc_per_type2_dual"] = std::make_shared<Vec1Int>();
                data["descriptor_list2_dual"] = std::make_shared<Vec2Int>();
                data["descriptor_list2_dual"].dget<ShVec2Int>().resize(ntypes);
            }

            // In the following everything is read in for each type separately
            for (Int itype = 0; itype < ntypes; ++itype)
            {
                // Read number of active desriptors
                read_values::io_array_binary(
                    input_stream,
                    data["SHS3-3-body-number-descriptors-per-type"].dget<ShVec1Int>(),
                    1);

                // Read active descriptors
                read_values::io_array_binary(
                    input_stream,
                    data["SHS3-3-body-descriptor-list"].dget<ShVec2Int>()[itype],
                    data["SHS3-3-body-number-descriptors-per-type"].dcget<ShVec1Int>()[itype]);

                // Reading variables for dual descriptor, same as above
                if (descriptorType == 2 || descriptorType == 4)
                {
                    read_values::io_array_binary(input_stream,
                                                 data["ndesc_per_type2_dual"].dget<ShVec1Int>(),
                                                 1);

                    read_values::io_array_binary(
                        input_stream,
                        data["descriptor_list2_dual"].dget<ShVec2Int>()[itype],
                        data["ndesc_per_type2_dual"].dcget<ShVec1Int>()[itype]);
                }
            } // closing number-types

            // Reading stuff related to SIC correction, we leave this
            // for the moment empty since fast version won't need SIC
            if (data["ML_LSIC"].cget<bool>())
            {
                global_scope::tutor.error(
                    "ERROR: This force field file does not support the fast "
                    "prediction mode. Please do a refitting with ML_MODE=REFIT.");
                //throw std::runtime_error(
                //      "ERROR: This force field file does not support the fast "
                //         "prediction mode. Please do a refitting with ML_MODE=REFIT."
                //      );
            }

        } // closing ML_W2

        // Read number of local reference configurations (basis sets in old code) per type
        // nlc
        data["number-local-reference-configs"] = std::make_shared<Vec1Int>();
        read_values::io_array_binary(input_stream,
                                     data["number-local-reference-configs"].dget<ShVec1Int>(),
                                     ntypes);
        const Vec1Int nlrc = data["number-local-reference-configs"].dcget<ShVec1Int>();

        // Read ML_ISCALE_TOTEN - type to scale total energy
        data["ML_ISCALE_TOTEN"] = (Int)(0);
        read_values::io_scalar(input_stream, data["ML_ISCALE_TOTEN"].get<Int>());

        // Read average energy per atom (eV atom^-1)
        data["average-energy-per-atom"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["average-energy-per-atom"].get<Real>());

        // Regression coefficients: Make a 2d vector as above and loop
        // over types (ev/atom)
        data["regression-coeff"] = std::make_shared<Vec2Real>();
        if (rankReadsOn)
        {
            // only read on node==0 zero since the array might be stored
            // in shared memory
            data["regression-coeff"].dget<ShVec2Real>().resize(ntypes);
            for (Int itype = 0; itype < ntypes; ++itype)
            {
                read_values::io_array_binary(input_stream,
                                             data["regression-coeff"].dget<ShVec2Real>()[itype],
                                             nlrc[itype]);
            }
        }

        // The following is only done if radial descriptor exists
        if (data["ML_W1"].cget<Real>() > 0.0)
        {
            // Number of descriptors for radial descriptor. It is the
            // same for each type hence it is only saved in a scalar.
            // ndesc_per_type1
            data["SHS2-2-body-number-descriptors-per-type"] = data["ML_MRB1"].cget<Int>() * ntypes;

            /// TODO: It would be ideal to read it in as -
            /// std::shared_ptr<std::vector<std::shared_ptr<std::vector<Real>>>>
            /// or in shared memory container directly
            /// rad_desc
            data["SHS2-2-body-reference-configs"] = std::make_shared<Vec2Real>();
            if (rankReadsOn)
            {
                // only read on node==0 zero since the array might be stored
                // in shared memory
                data["SHS2-2-body-reference-configs"].dget<ShVec2Real>().resize(ntypes);
                for (Int itype = 0; itype < ntypes; ++itype)
                {
                    read_values::io_array_binary(
                        input_stream,
                        data["SHS2-2-body-reference-configs"].dget<ShVec2Real>()[itype],
                        nlrc[itype] * data["SHS2-2-body-number-descriptors-per-type"].cget<Int>());
                }
            }
        }

        data["SHS3-3-body-reference-configs"] = std::make_shared<Vec2Real>();
        // The following is only done if angular descriptor exists
        if (data["ML_W2"].cget<Real>() > 0.0)
        {
            /// TODO: It would be ideal to read it in as -
            /// std::shared_ptr<std::vector<std::shared_ptr<std::vector<Real>>>>
            /// or in shared memory container directly
            if (rankReadsOn)
            {
                // only read on node==0 zero since the array might be stored
                // in shared memory
                data["SHS3-3-body-reference-configs"].dget<ShVec2Real>().resize(ntypes);
                for (Int itype = 0; itype < ntypes; ++itype)
                {
                    read_values::io_array_binary(
                        input_stream,
                        data["SHS3-3-body-reference-configs"].dget<ShVec2Real>()[itype],
                        nlrc[itype]
                            * data["SHS3-3-body-number-descriptors-per-type"]
                                  .dcget<ShVec1Int>()[itype]);
                }
            }
            // Reading variables for dual descriptor, same as above
            if (descriptorType == 2 || descriptorType == 4)
            {
                /// TODO: It would be ideal to read it in as -
                /// std::shared_ptr<std::vector<std::shared_ptr<std::vector<Real>>>>
                /// or in shared memory container directly
                data["ang_desc_dual"] = std::make_shared<Vec2Real>();
                data["ang_desc_dual"].dget<ShVec2Real>().resize(ntypes);

                for (Int itype = 0; itype < ntypes; ++itype)
                {
                    read_values::io_array_binary(
                        input_stream,
                        data["ang_desc_dual"].dget<ShVec2Real>()[itype],
                        nlrc[itype] * data["ndesc_per_type2_dual"].dcget<ShVec1Int>()[itype]);
                }
            }
        } // closing ML_W2

        // Read normalized noise parameter (-) not to be confused with ML_SIGV0
        data["sigv"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["sigv"].get<Real>());

        // Read normalized precision (unceartanty) parameter (-) not to be confused with ML_SIGV0
        data["sigw"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["sigw"].get<Real>());

        // Read variance in energies per atom in training data (eV atom^-1)
        data["variance-training-energies"] = (Real)(0.0);
        read_values::io_scalar(input_stream, data["variance-training-energies"].get<Real>());

        // Read variance in forces in training data (eV Angst^-1) (3 cartesian directions)
        data["variance-training-force"] = std::make_shared<Vec1Real>();
        read_values::io_array_binary(input_stream,
                                     data["variance-training-force"].dget<ShVec1Real>(),
                                     3);

        // Read variance in stresses in training data (eV Angst^-1) (6 independent terms)
        data["variance-training-stress"] = std::make_shared<Vec1Real>();
        read_values::io_array_binary(input_stream,
                                     data["variance-training-stress"].dget<ShVec1Real>(),
                                     6);

        // ML_LFAST is present here since 0.2.4.
        data["ML_LFAST"] = true;
        if (versionIn >= SemanticVersion(0, 2, 4))
        {
            read_values::io_scalar(input_stream, data["ML_LFAST"].get<bool>());
            if (!data["ML_LFAST"].cget<bool>())
            {
                global_scope::tutor.error(
                    "ERROR: This force field file does not support the fast "
                    "prediction mode. Please do a refitting with ML_MODE=REFIT.");
            }
        }

        if (rankReadsOn)
        {
            data["inverse-cov-matrix"] = std::make_shared<Vec2Real>();
            data["inverse-cov-matrix"].dget<ShVec2Real>().resize(ntypes);
            if (versionIn >= SemanticVersion(0, 2, 3) && data["ML_LFAST"].cget<bool>())
            {
                for (std::size_t type0 = 0; type0 < (std::size_t)ntypes; type0++)
                {
                    for (std::size_t locRef0 = 0;
                         locRef0 < (std::size_t)data["number-local-reference-configs"]
                                       .dcget<ShVec1Int>()[type0];
                         locRef0++)
                    {
                        for (std::size_t locRef1 = 0;
                             locRef1 < (std::size_t)data["number-local-reference-configs"]
                                           .dcget<ShVec1Int>()[type0];
                             locRef1++)
                        {
                            Real value;
                            read_values::io_scalar(input_stream, value);
                            data["inverse-cov-matrix"].dget<ShVec2Real>()[type0].push_back(value);
                        }
                    }
                }
            }
        }
    } // closing if LMLMB

    // close ML_FF file
    input_stream.close();

    // TODO: Here must come a post-read sanity check. For example, to disable newer features which
    // are not supported by older ML_FF versions. Put this into a separate function.
    if (versionIn < SemanticVersion(0, 2, 4))
    {
        // TODO: Disable spilling factor, warn user.
    }

    return;
}

const MultiTypeMap& IoHandlerML_FF::get_data(void) const
{
    return data;
}

// this only reads header of ML_FF into a buffer
void IoHandlerML_FF::read_header()
{
    char header_buffer[4096];
    input_stream.getline(header_buffer, sizeof(header_buffer));
    header = header_buffer;

    // no header found
    if (header.substr(0, 5) != "ML_FF")
    {
        global_scope::tutor.error(
            "ERROR: Unknown file format, input file does not seem to contain a "
            "valid header.");
        //throw std::runtime_error(
        //      "ERROR: Unknown file format, input file does not seem to contain a "
        //          "valid header."
        //      );
    }
}

void IoHandlerML_FF::check_version_number_header_and_binary()
{
    // Version string should start first six characters, search for space afterwards
    // (e.g. "ML_FF 1.2.3 binary").
    std::size_t i = header.substr(6).find(" ");
    if (i == 0 && i == std::string::npos)
    {
        throw std::runtime_error("ERROR: Invalid header, could not find version string.");
    }
    SemanticVersion versionHeader(header.substr(6, 6 + i));

    // Throw error if version numbers do not agree.
    if (versionHeader != versionIn)
    {
        throw std::runtime_error("ERROR: Binary version number (" + versionIn.toString()
                                 + ") and header version number (" + versionHeader.toString()
                                 + ") do not match.");
    }

    return;
}

const MultiType& IoHandlerML_FF::operator[](const std::string& key) const
{
    VASPML_DEBUG_L1(
        if (data.count(translator.get_InputTag(key)) == 0)
        {
            throw std::runtime_error(
                "ERROR in const MultiType& IoHandlerML_FF::operator[] ( const std::string& key ) \n"
                "      key does not exist in dictionary data");
        }
    );
    return data.at(translator.get_InputTag(key));
}

const MultiType& IoHandlerML_FF::operator()(const std::string& key) const
{
    VASPML_DEBUG_L1(
        if (data.count(key) == 0)
        {
            throw std::runtime_error(
                "ERROR in const MultiType& IoHandlerML_FF::operator() ( const std::string& key ) \n"
                "      key does not exist in dictionary data");
        }
    );
    return data.at(key);
}

const TagTranslator& IoHandlerML_FF::get_translator(void) const
{
    return translator;
}

std::map<String, Real> IoHandlerML_FF::generate_weightsMap(void) const
{
    std::map<String, Real> weights;
    for (const String& key : constants::descriptorKeyList)
    {
        std::string tag = key + "-weight";
        weights[key] = data.at(translator.get_InputTag(tag)).cget<Real>();
    }
    return weights;
}

std::map<String, ShVec1Int> IoHandlerML_FF::generate_featureSpaceMap(void) const
{
    std::map<String, ShVec1Int> featureSpaceSize;
    std::map<String, ShVec1Int> featureSpaceSize_check;
    const Vec1Int&              locRef =
        data.at(translator.get_InputTag("number-local-reference-configs")).dcget<ShVec1Int>();
    for (const String& key : constants::descriptorKeyList)
    {
        std::string tag = key + "-number-descriptors-per-type";
        if (std::get_if<Int>(&data.at(translator.get_InputTag(tag))))
        {
            featureSpaceSize[key] = std::make_shared<Vec1Int>(
                math::vectorTimesScalar(locRef, data.at(translator.get_InputTag(tag)).cget<Int>()));
        }
        else if (std::get_if<ShVec1Int>(&data.at(translator.get_InputTag(tag))))
        {
            featureSpaceSize[key] = std::make_shared<Vec1Int>(
                math::elementwiseProduct(locRef,
                                         data.at(translator.get_InputTag(tag)).dcget<ShVec1Int>()));
        }
    }
    return featureSpaceSize;
}

std::map<String, ShVec1Int> IoHandlerML_FF::get_numberDescriptorsMap(void) const
{
    std::map<String, ShVec1Int> numberDescriptors;
    numberDescriptors["SHS3-3-body"] = std::make_shared<Vec1Int>(
        data.at(translator.get_InputTag("SHS3-3-body-number-descriptors-per-type"))
            .dcget<ShVec1Int>());
    Int twoBody =
        data.at(translator.get_InputTag("SHS2-2-body-number-descriptors-per-type")).cget<Int>();
    Vec1Int twoBodyVector(numberDescriptors["SHS3-3-body"]->size());
    for (std::size_t i = 0; i < twoBodyVector.size(); i++) { twoBodyVector[i] = twoBody; }
    numberDescriptors["SHS2-2-body"] = std::make_shared<Vec1Int>(twoBodyVector);
    return numberDescriptors;
}

Real IoHandlerML_FF::maxCutoffRadius() const
{
    return std::max((*this)["SHS2-2-body-cutoff"].cget<Real>(),
                    (*this)["SHS3-3-body-cutoff"].cget<Real>());
}

void IoHandlerML_FF::convertUnitsToAtomicUnits(void)
{

    // EATOM_REF
    math::vectorTimesScalarNoCopy(
        data[translator.get_InputTag("reference-energies")].dget<ShVec1Real>(),
        1.0 / constants::EUNIT);

    // FF%WTOTEN
    data[translator.get_InputTag("energy-weight")].get<Real>() /= constants::EUNIT;
    // FF%WTIFOR
    data[translator.get_InputTag("force-weight")].get<Real>() /= constants::FUNIT;
    // FF%WTSIF
    data[translator.get_InputTag("stress-weight")].get<Real>() /= constants::SUNIT;

    // FF%STOTEN
    data[translator.get_InputTag("rmse-total-energy")].get<Real>() /= constants::EUNIT;
    // FF%STIFOR
    data[translator.get_InputTag("rmse-forces")].get<Real>() /= constants::FUNIT;
    // FF%STSIF
    data[translator.get_InputTag("rmse-stresses")].get<Real>() /= constants::SUNIT;
    // FFM%RCUT1
    data[translator.get_InputTag("SHS2-2-body-cutoff")].get<Real>() /= constants::AUTOA;
    // FFM%RCUT2
    data[translator.get_InputTag("SHS3-3-body-cutoff")].get<Real>() /= constants::AUTOA;
    // FFM%SION1
    data[translator.get_InputTag("SHS2-2-body-smearing-param")].get<Real>() /= constants::AUTOA;
    // FFM%SION2
    data[translator.get_InputTag("SHS3-3-body-smearing-param")].get<Real>() /= constants::AUTOA;
    // FFM%TOTENAV
    data[translator.get_InputTag("average-energy-per-atom")].get<Real>() /= constants::EUNIT;
    // FFM%WMAT
    if (rankReadsOn)
    {
        math::vectorTimesScalarNoCopy(
            data[translator.get_InputTag("regression-coeff")].dget<ShVec2Real>(),
            (Real)1.0 / constants::EUNIT);
    }
    // FFM%SIG(1)
    data[translator.get_InputTag("variance-training-energies")].get<Real>() /= constants::EUNIT;
    // FFM%SIG(2:4)
    math::vectorTimesScalarNoCopy(
        data[translator.get_InputTag("variance-training-force")].dget<ShVec1Real>(),
        1.0 / constants::FUNIT);
    // FFM%SIG(5:10)
    math::vectorTimesScalarNoCopy(
        data[translator.get_InputTag("variance-training-stress")].dget<ShVec1Real>(),
        1.0 / constants::SUNIT);
}
