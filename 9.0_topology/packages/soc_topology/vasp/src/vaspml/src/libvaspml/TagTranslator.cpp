#include "TagTranslator.hpp"
#include "Tutor.hpp"
#include "debug.hpp"

#include <stdexcept>

using namespace vaspml;

TagTranslator::TagTranslator(void)
{
    mapToInputTags["descriptor-type"] = "ML_DESC_TYPE";
    mapToInputTags["number-types"] = "number-types";
    mapToInputTags["types"] = "types";
    mapToInputTags["reference-energies"] = "ML_EATOM_REF";
    mapToInputTags["atomic-masses"] = "atomic-masses";
    mapToInputTags["data-weighing-scheme"] = "ML_IWEIGHT";
    mapToInputTags["energy-weight"] = "ML_WTOTEN";
    mapToInputTags["force-weight"] = "ML_WTIFOR";
    mapToInputTags["stress-weight"] = "ML_WTSIF";
    mapToInputTags["rmse-total-energy"] = "rmse-total-energy";
    mapToInputTags["rmse-forces"] = "rmse-forces";
    mapToInputTags["rmse-stresses"] = "rmse-stresses";
    // not used is a tag which's value is always true
    mapToInputTags["many-body"] = "many-body";
    // not used in run mode fast
    mapToInputTags["SIC-on-off"] = "ML_LSIC";
    // supervector variale will always be true for execution mode
    mapToInputTags["supervector-on-off"] = "ML_LSUPERVEC";
    mapToInputTags["SHS2-2-body-weight"] = "ML_W1";
    mapToInputTags["SHS3-3-body-weight"] = "ML_W2";
    mapToInputTags["SHS2-2-body-cutoff-type"] = "ML_ICUT1";
    mapToInputTags["SHS3-3-body-cutoff-type"] = "ML_ICUT2";
    mapToInputTags["SHS2-2-body-cutoff"] = "ML_RCUT1";
    mapToInputTags["SHS3-3-body-cutoff"] = "ML_RCUT2";
    mapToInputTags["SHS2-2-body-smearing-type"] = "ML_IBROAD1";
    mapToInputTags["SHS3-3-body-smearing-type"] = "ML_IBROAD2";
    mapToInputTags["SHS2-2-body-smearing-param"] = "ML_SION1";
    mapToInputTags["SHS3-3-body-smearing-param"] = "ML_SION2";
    mapToInputTags["SHS2-2-body-max-radial-funcs"] = "ML_MRB1";
    mapToInputTags["SHS3-3-body-max-radial-funcs"] = "ML_MRB2";
    // not used. Could be used as switch for basis functions in future
    mapToInputTags["SHS3-3-body-basis-func-type"] = "ML_BASIS_TYPE2";
    // not used
    mapToInputTags["ML_MRB2_MAX"] = "ML_MRB2_MAX";
    // not used
    mapToInputTags["SHS3-3-body-reduced-max-radial-funcs"] = "ML_MRB2_RED";
    // not used. Could be used as switch for basis functions in future
    mapToInputTags["SHS3-3-body-reduced-basis-type"] = "ML_BASIS_TYPE2_RED";
    // not used
    mapToInputTags["ML_MRB2_MAX_RED"] = "ML_MRB2_MAX_RED";
    mapToInputTags["SHS2-2-body-number-radial-grid-points"] = "ML_NR1";
    mapToInputTags["SHS3-3-body-number-radial-grid-points"] = "ML_NR2";
    // not used. instead ML_NR1 is used
    mapToInputTags["SHS2-2-body-number-spline-grid-points"] = "ML_MSPL1";
    // not used. instead ML_NR2 is used
    mapToInputTags["SHS3-3-body-number-spline-grid-points"] = "ML_MSPL2";
    // not used because the assiigned value is always zero
    mapToInputTags["SHS2-2-body-max-angular-number"] = "ML_LMAX1";
    mapToInputTags["SHS3-3-body-max-angular-number"] = "ML_LMAX2";
    // not used yet
    mapToInputTags["SHS3-3-body-reduced-max-angular-number"] = "ML_LMAX2_RED";
    mapToInputTags["ML_DESC_RATIO_DUAL"] = "ML_DESC_RATIO_DUAL";
    mapToInputTags["desc_ratio_dual_second"] = "desc_ratio_dual_second";
    mapToInputTags["ML_DESC_FACTORE_TESTE"] = "ML_DESC_FACTORE_TESTE";
    mapToInputTags["SHS2-2-body-exponent"] = "ML_NHYP1";
    mapToInputTags["SHS3-3-body-exponent"] = "ML_NHYP2";
    mapToInputTags["SHS2-2-body-is-normalized"] = "ML_LNORM1";
    mapToInputTags["SHS3-3-body-is-normalized"] = "ML_LNORM2";
    // not used
    mapToInputTags["ML_LWINDOW1"] = "ML_LWINDOW1";
    // not used
    mapToInputTags["ML_LWINDOW2"] = "ML_LWINDOW2";
    // not used
    mapToInputTags["ML_IWINDOW1"] = "ML_IWINDOW1";
    // not used
    mapToInputTags["ML_IWINDOW2"] = "ML_IWINDOW2";
    mapToInputTags["SHS3-3-body-angular-filter-on"] = "ML_LAFILT2";
    mapToInputTags["SHS3-3-body-angular-filter-type"] = "ML_IAFILT2";
    mapToInputTags["SHS3-3-body-angular-filter-scale"] = "ML_AFILT2";
    // not used
    mapToInputTags["ML_LMETRIC1"] = "ML_LMETRIC1";
    // not used
    mapToInputTags["ML_LMETRIC2"] = "ML_LMETRIC2";
    // not used
    mapToInputTags["ML_NMETRIC1"] = "ML_NMETRIC1";
    // not used
    mapToInputTags["ML_RMETRIC1"] = "ML_RMETRIC1";
    // not used
    mapToInputTags["ML_NMETRIC2"] = "ML_NMETRIC2";
    // not used
    mapToInputTags["ML_RMETRIC2"] = "ML_RMETRIC2";
    // not used
    mapToInputTags["ML_LVARTRAN1"] = "ML_LVARTRAN1";
    // not used
    mapToInputTags["ML_LVARTRAN2"] = "ML_LVARTRAN2";
    // not used
    mapToInputTags["ML_NVARTRAN1"] = "ML_NVARTRAN1";
    // not used
    mapToInputTags["ML_NVARTRAN2"] = "ML_NVARTRAN2";
    // not used
    mapToInputTags["SHS3-3-body-number-radial-basis"] = "SHS3-3-body-number-radial-basis";
    // not used
    mapToInputTags["SHS3-3-body-reduced-number-radial-basis"] =
        "SHS3-3-body-reduced-number-radial-basis";
    mapToInputTags["SHS2-2-body-number-descriptors-per-type"] =
        "SHS2-2-body-number-descriptors-per-type";
    mapToInputTags["SHS3-3-body-number-descriptors-per-type"] =
        "SHS3-3-body-number-descriptors-per-type";
    mapToInputTags["SHS3-3-body-descriptor-list"] = "SHS3-3-body-descriptor-list";
    // not used
    mapToInputTags["ndesc_per_type2_dual"] = "ndesc_per_type2_dual";
    // not used
    mapToInputTags["descriptor_list2_dual"] = "descriptor_list2_dual";
    mapToInputTags["number-local-reference-configs"] = "number-local-reference-configs";
    mapToInputTags["SHS2-2-body-reference-configs"] = "SHS2-2-body-reference-configs";
    mapToInputTags["SHS3-3-body-reference-configs"] = "SHS3-3-body-reference-configs";
    mapToInputTags["scale-type-total-energy"] = "ML_ISCALE_TOTEN";
    mapToInputTags["average-energy-per-atom"] = "average-energy-per-atom";
    mapToInputTags["regression-coeff"] = "regression-coeff";
    // not used
    mapToInputTags["sigma-v"] = "sigv";
    // not used
    mapToInputTags["sigma-w"] = "sigw";
    // not used
    mapToInputTags["fast"] = "ML_LFAST";
    mapToInputTags["variance-training-energies"] = "variance-training-energies";
    mapToInputTags["variance-training-force"] = "variance-training-force";
    mapToInputTags["variance-training-stress"] = "variance-training-stress";
    mapToInputTags["inverse-cov-matrix"] = "inverse-cov-matrix";

    for (const auto& [key, value] : mapToInputTags) { mapToCodeTag[value] = key; }
}

std::string TagTranslator::get_InputTag(const std::string& tag) const
{
    VASPML_DEBUG_L1(
        if (mapToInputTags.count(tag) == 0)
        {
            throw std::runtime_error(
                "ERROR in std::string TagTranslator::get_InputTag( const std::string& tag ) \n"
                "      tag "
                + tag + " does not exist in dictionary data");
        }
    );
    return mapToInputTags.at(tag);
}

std::string TagTranslator::get_CodeTag(const std::string& tag) const
{
    VASPML_DEBUG_L1(
        if (mapToCodeTag.count(tag) == 0)
        {
            throw std::runtime_error(
                "ERROR in std::string TagTranslator::get_CodeTag( const std::string& tag ) \n"
                "      tag "
                + tag + " does not exist in dictionary data");
        }
    );
    return mapToCodeTag.at(tag);
}
