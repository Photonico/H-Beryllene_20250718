#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE io_ML_FF

#include "IoHandlerML_FF.hpp"
#include "boost_helpers.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <boost/test/unit_test.hpp>
#include <vector>

using namespace vaspml;

Real const tolerance = 10.0 * std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(InitializeWrite_ML_FF_and_read)
{

    // initialize filename and class for reading
    std::string file_name = "ML_FF";

    // first make an output file stream
    std::ofstream output_stream;
    output_stream = file_io::openFileO(file_name, std::ofstream::binary);

    // write header
    std::string letter;
    letter = "ML_FF 0.2.1 binary";
    output_stream << letter;
    letter = " ";
    for (Int i = 18; i < 4095; i++) { output_stream << letter; }
    letter = "\n";
    output_stream << letter;

    // from here on we set data for writing and write
    Int v1 = 0;
    Int v2 = 2;
    Int v3 = 1;
    read_values::io_scalar(output_stream, v1);
    read_values::io_scalar(output_stream, v2);
    read_values::io_scalar(output_stream, v3);

    Int ML_DESC_TYPE = 0;
    read_values::io_scalar(output_stream, ML_DESC_TYPE);

    Int ntype = 2;
    read_values::io_scalar(output_stream, ntype);

    std::vector<std::string> types;
    types.push_back("Li");
    types.push_back("H");
    read_values::io_array_binary(output_stream, types, ntype);

    std::vector<Real> reference_energies;
    reference_energies.push_back(0.0);
    reference_energies.push_back(0.0);
    read_values::io_array_binary(output_stream, reference_energies, ntype);

    std::vector<Real> atomic_masses;
    atomic_masses.push_back(7.01);
    atomic_masses.push_back(1.0);
    read_values::io_array_binary(output_stream, atomic_masses, ntype);

    Int ML_IWEIGHT = 3;
    read_values::io_scalar(output_stream, ML_IWEIGHT);

    Real ML_WTOTEN = 1.0;
    read_values::io_scalar(output_stream, ML_WTOTEN);

    Real ML_WTIFOR = 1.0;
    read_values::io_scalar(output_stream, ML_WTIFOR);

    Real ML_WTSIF = 1.0;
    read_values::io_scalar(output_stream, ML_WTSIF);

    Real mse_toten = 2.13557005E-04;
    read_values::io_scalar(output_stream, mse_toten);

    Real mse_forces = 6.52952959E-03;
    read_values::io_scalar(output_stream, mse_forces);

    Real mse_stresses = 3.92703803E-01;
    read_values::io_scalar(output_stream, mse_stresses);

    bool many_body = true;
    read_values::io_scalar(output_stream, many_body);

    bool ML_LSIC = false;
    read_values::io_scalar(output_stream, ML_LSIC);

    bool ML_LSUPERVEC = true;
    read_values::io_scalar(output_stream, ML_LSUPERVEC);

    Real ML_W1 = 0.1;
    read_values::io_scalar(output_stream, ML_W1);

    Real ML_W2 = 0.9;
    read_values::io_scalar(output_stream, ML_W2);

    Int ML_ICUT1 = 1;
    read_values::io_scalar(output_stream, ML_ICUT1);

    Int ML_ICUT2 = 1;
    read_values::io_scalar(output_stream, ML_ICUT2);

    Real ML_RCUT1 = 8.0;
    read_values::io_scalar(output_stream, ML_RCUT1);

    Real ML_RCUT2 = 5.0;
    read_values::io_scalar(output_stream, ML_RCUT2);

    Int ML_IBROAD1 = 2;
    read_values::io_scalar(output_stream, ML_IBROAD1);

    Int ML_IBROAD2 = 2;
    read_values::io_scalar(output_stream, ML_IBROAD2);

    Real ML_SION1 = 0.5;
    read_values::io_scalar(output_stream, ML_SION1);

    Real ML_SION2 = 0.5;
    read_values::io_scalar(output_stream, ML_SION2);

    Int ML_MRB1 = 3;
    read_values::io_scalar(output_stream, ML_MRB1);

    Int ML_MRB2 = 8;
    read_values::io_scalar(output_stream, ML_MRB2);

    Int ML_NR1 = 1000;
    read_values::io_scalar(output_stream, ML_NR1);

    Int ML_NR2 = 1000;
    read_values::io_scalar(output_stream, ML_NR2);

    Int ML_MSPL1 = 1000;
    read_values::io_scalar(output_stream, ML_MSPL1);

    Int ML_MSPL2 = 1000;
    read_values::io_scalar(output_stream, ML_MSPL2);

    Int ML_LMAX1 = 0;
    read_values::io_scalar(output_stream, ML_LMAX1);

    Int ML_LMAX2 = 3;
    read_values::io_scalar(output_stream, ML_LMAX2);

    Int ML_NHYP1 = 4;
    read_values::io_scalar(output_stream, ML_NHYP1);

    Int ML_NHYP2 = 4;
    read_values::io_scalar(output_stream, ML_NHYP2);

    bool ML_LNORM1 = true;
    read_values::io_scalar(output_stream, ML_LNORM1);

    bool ML_LNORM2 = true;
    read_values::io_scalar(output_stream, ML_LNORM2);

    bool ML_LWINDOW1 = false;
    read_values::io_scalar(output_stream, ML_LWINDOW1);

    bool ML_LWINDOW2 = false;
    read_values::io_scalar(output_stream, ML_LWINDOW2);

    bool ML_LAFILT2 = true;
    read_values::io_scalar(output_stream, ML_LAFILT2);

    Int ML_IAFILT2 = 2;
    read_values::io_scalar(output_stream, ML_IAFILT2);

    Real ML_AFILT2 = 2.00000E-03;
    read_values::io_scalar(output_stream, ML_AFILT2);

    bool ML_LMETRIC1 = false;
    read_values::io_scalar(output_stream, ML_LMETRIC1);

    bool ML_LMETRIC2 = false;
    read_values::io_scalar(output_stream, ML_LMETRIC2);

    bool ML_LVARTRAN1 = false;
    read_values::io_scalar(output_stream, ML_LVARTRAN1);

    bool ML_LVARTRAN2 = false;
    read_values::io_scalar(output_stream, ML_LVARTRAN2);

    Vec1Int nrb2;
    nrb2.push_back(8);
    nrb2.push_back(7);
    nrb2.push_back(7);
    nrb2.push_back(6);
    read_values::io_array_binary(output_stream, nrb2, ML_LMAX2 + 1);

    Vec1Int ndesc_per_type2;
    ndesc_per_type2.push_back(2);
    ndesc_per_type2.push_back(2);

    Vec2Int descriptor_list2;
    descriptor_list2.resize(ntype);
    descriptor_list2[0].push_back(1);
    descriptor_list2[0].push_back(2);
    descriptor_list2[1].push_back(1);
    descriptor_list2[1].push_back(2);

    for (Int itype = 0; itype < ntype; itype++)
    {
        read_values::io_scalar(output_stream, ndesc_per_type2[itype]);

        read_values::io_array_binary(output_stream,
                                     descriptor_list2[itype],
                                     ndesc_per_type2[itype]);
    }

    Vec1Int nlc;
    nlc.push_back(3);
    nlc.push_back(3);
    read_values::io_array_binary(output_stream, nlc, ntype);

    Int ML_ISCALE_TOTEN = 2;
    read_values::io_scalar(output_stream, ML_ISCALE_TOTEN);

    Real av_en_per_atom = -2.9298355732899264;
    read_values::io_scalar(output_stream, av_en_per_atom);

    Vec2Real reg_coeff;
    reg_coeff.resize(ntype);
    reg_coeff[0].push_back(-686.69536154642321);
    reg_coeff[0].push_back(-865.10562024210083);
    reg_coeff[0].push_back(990.78098090149274);
    reg_coeff[1].push_back(224.18618969343643);
    reg_coeff[1].push_back(2337.4403403029360);
    reg_coeff[1].push_back(2042.5068756847654);
    for (Int itype = 0; itype < ntype; itype++)
    {
        read_values::io_array_binary(output_stream, reg_coeff[itype], nlc[itype]);
    }

    Int ndesc_per_type1 = ML_MRB1 * ntype;

    Vec2Real rad_desc;
    rad_desc.resize(ntype);
    rad_desc[0].push_back(0.22978615754061188);
    rad_desc[0].push_back(2.0579589580524006E-002);
    rad_desc[0].push_back(-5.1934144181894383E-002);
    rad_desc[0].push_back(-4.0534776586418166E-002);
    rad_desc[0].push_back(-3.9764923963422301E-002);
    rad_desc[0].push_back(5.1997688066030756E-003);
    rad_desc[0].push_back(3.3764626192374560E-002);
    rad_desc[0].push_back(1.5768234758278408E-002);
    rad_desc[0].push_back(-6.1622963101823905E-003);
    rad_desc[0].push_back(-1.0582779739814221E-002);
    rad_desc[0].push_back(-3.2312113434266438E-003);
    rad_desc[0].push_back(2.0250885854626332E-003);
    rad_desc[0].push_back(0.24615861497264371);
    rad_desc[0].push_back(5.1751610814095086E-002);
    rad_desc[0].push_back(-1.2726067106775785E-002);
    rad_desc[0].push_back(3.2243806386104829E-003);
    rad_desc[0].push_back(-7.0299097133974366E-003);
    rad_desc[0].push_back(-5.2878912400089809E-002);
    rad_desc[1].push_back(0.24229135702269469);
    rad_desc[1].push_back(5.0947900071248615E-002);
    rad_desc[1].push_back(-1.4243153389865208E-002);
    rad_desc[1].push_back(9.2006451934847542E-004);
    rad_desc[1].push_back(-4.5633726331625679E-003);
    rad_desc[1].push_back(-4.4319637106942984E-002);
    rad_desc[1].push_back(-4.1588438603423747E-002);
    rad_desc[1].push_back(7.7065567125702995E-003);
    rad_desc[1].push_back(1.0940168151592565E-002);
    rad_desc[1].push_back(7.5800866065998578E-003);
    rad_desc[1].push_back(1.0122071147858201E-002);
    rad_desc[1].push_back(2.2944268013535798E-004);
    rad_desc[1].push_back(0.22621575291956064);
    rad_desc[1].push_back(2.0144970692080818E-002);
    rad_desc[1].push_back(-5.3198956098150160E-002);
    rad_desc[1].push_back(-4.5252310975261820E-002);
    rad_desc[1].push_back(-4.8150459295592619E-002);
    rad_desc[1].push_back(8.5785445491588889E-003);
    for (Int itype = 0; itype < ntype; itype++)
    {
        read_values::io_array_binary(output_stream, rad_desc[itype], nlc[itype] * ndesc_per_type1);
    }

    Vec2Real ang_desc;
    ang_desc.resize(ntype);
    ang_desc[0].push_back(0.25525218849981657);
    ang_desc[0].push_back(-0.12754185410327940);
    ang_desc[0].push_back(-0.15287888624545887);
    ang_desc[0].push_back(7.2229350179085619E-002);
    ang_desc[0].push_back(4.8612164015436415E-002);
    ang_desc[0].push_back(-3.0738888156802379E-002);
    ang_desc[1].push_back(0.45471692940332076);
    ang_desc[1].push_back(0.12156331310372255);
    ang_desc[1].push_back(-0.15602082512783186);
    ang_desc[1].push_back(-0.19639025564030410);
    ang_desc[1].push_back(-2.0856127255728975E-002);
    ang_desc[1].push_back(7.2422878262464158E-002);
    for (Int itype = 0; itype < ntype; itype++)
    {
        read_values::io_array_binary(output_stream,
                                     ang_desc[itype],
                                     nlc[itype] * ndesc_per_type2[itype]);
    }

    Real sigv = 1.0;
    read_values::io_scalar(output_stream, sigv);

    Real sigw = 1e-07;
    read_values::io_scalar(output_stream, sigw);

    Real train_var_energy = 2.3315240295174282E-002;
    read_values::io_scalar(output_stream, train_var_energy);

    Vec1Real train_var_force;
    train_var_force.push_back(6.0799749030095119E-003);
    train_var_force.push_back(7.8956271398709435E-003);
    train_var_force.push_back(5.0997646192182319E-003);
    read_values::io_array_binary(output_stream, train_var_force, 3);

    Vec1Real train_var_stress;
    train_var_stress.push_back(16.036727642191934);
    train_var_stress.push_back(4.3523261400672242);
    train_var_stress.push_back(4.4873194025469640);
    train_var_stress.push_back(15.676315403732437);
    train_var_stress.push_back(7.3886773570760589);
    train_var_stress.push_back(16.924959396332085);
    read_values::io_array_binary(output_stream, train_var_stress, 6);

    output_stream.close();
    //
    ////---------------------------------------------------------------------------------
    //
    // read everything we just wrote
    std::ifstream input_stream;
    input_stream = file_io::openFileI(file_name, std::ifstream::binary);

    IoHandlerML_FF ml_ff_handler(file_name);
    ml_ff_handler.mainReader();

    input_stream.close();

    Int  ntype_new = std::get<Int>(ml_ff_handler["number-types"]);
    Real train_var_energy_new = ml_ff_handler.get_scalar<Real>("variance-training-energies");
    Int  ML_MRB1_new = ml_ff_handler.get_scalar<Int>("SHS2-2-body-max-radial-funcs");
    Int  ndesc_per_type1_new = ML_MRB1_new * ntype_new;

    const Vec1String& types_new = ml_ff_handler.get_array<String>("types");
    const Vec2Real& rad_desc_new = ml_ff_handler.get_2Darray<Real>("SHS2-2-body-reference-configs");
    const Vec2Real& ang_desc_new = ml_ff_handler.get_2Darray<Real>("SHS3-3-body-reference-configs");
    const Vec2Real& reg_coeff_new = ml_ff_handler.get_2Darray<Real>("regression-coeff");
    const Vec1Int&  nlc_new = ml_ff_handler.get_array<Int>("number-local-reference-configs");
    const Vec1Int&  ndesc_per_type2_new =
        ml_ff_handler.get_array<Int>("SHS3-3-body-number-descriptors-per-type");
    const Vec2Int& descriptor_list2_new =
        ml_ff_handler.get_2Darray<Int>("SHS3-3-body-descriptor-list");

    BOOST_REQUIRE_EQUAL(ntype_new, ntype);
    BOOST_REQUIRE_EQUAL(ndesc_per_type1_new, ndesc_per_type1);
    // Maybe CLOSE is better?
    BOOST_REQUIRE_CLOSE(train_var_energy_new, train_var_energy, tolerance);

    REQUIRE_EQUAL_COLLECTIONS(types_new, types, "types");
    REQUIRE_EQUAL_COLLECTIONS(nlc_new, nlc, "number-local-reference-configs");
    REQUIRE_EQUAL_COLLECTIONS(ndesc_per_type2_new, ndesc_per_type2, "ndesc_per_type2_new");
    REQUIRE_EQUAL_COLLECTIONS_2D(descriptor_list2_new,
                                 descriptor_list2,
                                 "SHS3-3-body-descriptor-list");
    // Maybe CLOSE is better?
    REQUIRE_CLOSE_COLLECTIONS_2D(rad_desc_new,
                                 rad_desc,
                                 tolerance,
                                 "SHS2-2-body-reference-configs");
    REQUIRE_CLOSE_COLLECTIONS_2D(ang_desc_new,
                                 ang_desc,
                                 tolerance,
                                 "SHS3-3-body-reference-configs");
    REQUIRE_CLOSE_COLLECTIONS_2D(reg_coeff_new, reg_coeff, tolerance, "regression-coeff");
}

BOOST_AUTO_TEST_SUITE_END()
