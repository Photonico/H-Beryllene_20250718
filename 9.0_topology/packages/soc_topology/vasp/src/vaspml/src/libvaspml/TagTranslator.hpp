#ifndef TAGTRANSLATOR_HPP
#define TAGTRANSLATOR_HPP

#include <map>
#include <string>

namespace vaspml
{

/*******************************************************************************************
 * @class TagTranslator
 * @brief this class gives the possibility to convert input tags to code tags
 *
 * the class gives the possibility to convert input tags which are defined in
 * the IoHandlerML_FF. The values that exist also in the INCAR file of vasp
 * are assigned to the INCAR tags. For the code those tags might not be of practical
 * use because they are not systematic. The map gives the possibility to
 * convert the INCAR tag to a systematic tag which can be used in loops over
 * dictionaries. For values wich are not defined in the INCAR file but are stored
 * in the ML_FF file the map will give a unit mapping, converting a tg to istelf.
 * All tags that are defined in the IoHandlerML_FF can be found here.
 *
 * @note in the current code only the systematic VASPml tags are allowed to be used
 *
 * <table>
 * <caption id="multi_row">Tag comparison table</caption>
 * <tr><th> VASPml tag <th> traditional tag <th> is in use
 * <tr><th> descriptor-type <th> ML_DESC_TYPE <th> used
 * <tr><th> number-types <th> number-types <th> used
 * <tr><th> types <th> types <th> used
 * <tr><th> reference-energies <th> ML_EATOM_REF <th> used
 * <tr><th> atomic-masses <th> atomic-masses <th> used
 * <tr><th> data-weighing-scheme <th> ML_IWEIGHT <th> used
 * <tr><th> energy-weight <th> ML_WTOTEN <th> used
 * <tr><th> force-weight <th> ML_WTIFOR <th> used
 * <tr><th> stress-weight <th> ML_WTSIF <th> used
 * <tr><th> rmse-total-energy <th> rmse-total-energy <th> used
 * <tr><th> rmse-forces <th> rmse-forces <th> used
 * <tr><th> rmse-stresses <th> rmse-stresses <th> used
 * <tr><th> many-body <th> many-body <th> not used
 * <tr><th> SIC-on-off <th> ML_LSIC <th> not used
 * <tr><th> supervector-on-off <th> ML_LSUPERVEC <th> used
 * <tr><th> SHS2-2-body-weight <th> ML_W1 <th> used
 * <tr><th> SHS3-3-body-weight <th> ML_W2 <th> used
 * <tr><th> SHS2-2-body-cutoff-type <th> ML_ICUT1 <th> used
 * <tr><th> SHS3-3-body-cutoff-type <th> ML_ICUT2 <th> used
 * <tr><th> SHS2-2-body-cutoff <th> ML_RCUT1 <th> used
 * <tr><th> SHS3-3-body-cutoff <th> ML_RCUT2 <th> used
 * <tr><th> SHS2-2-body-smearing-type <th> ML_IBROAD1 <th> used
 * <tr><th> SHS3-3-body-smearing-type <th> ML_IBROAD2 <th> used
 * <tr><th> SHS2-2-body-smearing-param <th> ML_SION1 <th> used
 * <tr><th> SHS3-3-body-smearing-param <th> ML_SION2 <th> used
 * <tr><th> SHS2-2-body-max-radial-funcs <th> ML_MRB1 <th> used
 * <tr><th> SHS3-3-body-max-radial-funcs <th> ML_MRB2 <th> used
 * <tr><th> SHS3-3-body-basis-func-type <th> ML_BASIS_TYPE2 <th> not used
 * <tr><th> ML_MRB2_MAX <th> ML_MRB2_MAX <th> not used
 * <tr><th> SHS3-3-body-reduced-max-radial-funcs <th> ML_MRB2_RED <th> not used
 * <tr><th> SHS3-3-body-reduced-basis-type <th> ML_BASIS_TYPE2_RED <th> not used
 * <tr><th> ML_MRB2_MAX_RED <th> ML_MRB2_MAX_RED <th> not used
 * <tr><th> SHS2-2-body-number-radial-grid-points <th> ML_NR1 <th> used
 * <tr><th> SHS3-3-body-number-radial-grid-points <th> ML_NR2 <th> used
 * <tr><th> SHS2-2-body-number-spline-grid-points <th> ML_MSPL1 <th> not used
 * <tr><th> SHS3-3-body-number-spline-grid-points <th> ML_MSPL2 <th> not used
 * <tr><th> SHS2-2-body-max-angular-number <th> ML_LMAX1 <th> not used
 * <tr><th> SHS3-3-body-max-angular-number <th> ML_LMAX2 <th> used
 * <tr><th> SHS3-3-body-reduced-max-angular-number <th> ML_LMAX2_RED <th> not used
 * <tr><th> ML_DESC_RATIO_DUAL <th> ML_DESC_RATIO_DUAL <th> used
 * <tr><th> desc_ratio_dual_second <th> desc_ratio_dual_second <th> used
 * <tr><th> ML_DESC_FACTORE_TESTE <th> ML_DESC_FACTORE_TESTE <th> used
 * <tr><th> SHS2-2-body-exponent <th> ML_NHYP1 <th> used
 * <tr><th> SHS3-3-body-exponent <th> ML_NHYP2 <th> used
 * <tr><th> SHS2-2-body-is-normalized <th> ML_LNORM1 <th> used
 * <tr><th> SHS3-3-body-is-normalized <th> ML_LNORM2 <th> used
 * <tr><th> ML_LWINDOW1 <th> ML_LWINDOW1 <th> not used
 * <tr><th> ML_LWINDOW2 <th> ML_LWINDOW2 <th> not used
 * <tr><th> ML_IWINDOW1 <th> ML_IWINDOW1 <th> not used
 * <tr><th> ML_IWINDOW2 <th> ML_IWINDOW2 <th> not used
 * <tr><th> SHS3-3-body-angular-filter-on <th> ML_LAFILT2 <th> used
 * <tr><th> SHS3-3-body-angular-filter-type <th> ML_IAFILT2 <th> used
 * <tr><th> SHS3-3-body-angular-filter-scale <th> ML_AFILT2 <th> used
 * <tr><th> ML_LMETRIC1 <th> ML_LMETRIC1 <th> not used
 * <tr><th> ML_LMETRIC2 <th> ML_LMETRIC2 <th> not used
 * <tr><th> ML_NMETRIC1 <th> ML_NMETRIC1 <th> not used
 * <tr><th> ML_RMETRIC1 <th> ML_RMETRIC1 <th> not used
 * <tr><th> ML_NMETRIC2 <th> ML_NMETRIC2 <th> not used
 * <tr><th> ML_RMETRIC2 <th> ML_RMETRIC2 <th> not used
 * <tr><th> ML_LVARTRAN1 <th> ML_LVARTRAN1 <th> not used
 * <tr><th> ML_LVARTRAN2 <th> ML_LVARTRAN2 <th> not used
 * <tr><th> ML_NVARTRAN1 <th> ML_NVARTRAN1 <th> not used
 * <tr><th> ML_NVARTRAN2 <th> ML_NVARTRAN2 <th> not used
 * <tr><th> SHS3-3-body-number-radial-basis <th> SHS3-3-body-number-radial-basis <th> not used
 * <tr><th> SHS3-3-body-reduced-number-radial-basis <th> SHS3-3-body-reduced-number-radial-basis
 * <th> not used
 * <tr><th> SHS2-2-body-number-descriptors-per-type <th> SHS2-2-body-number-descriptors-per-type
 * <th> used
 * <tr><th> SHS3-3-body-number-descriptors-per-type <th> SHS3-3-body-number-descriptors-per-type
 * <th> used
 * <tr><th> SHS3-3-body-descriptor-list <th> SHS3-3-body-descriptor-list <th> used
 * <tr><th> ndesc_per_type2_dual <th> ndesc_per_type2_dual <th> not used
 * <tr><th> descriptor_list2_dual <th> descriptor_list2_dual <th> not used
 * <tr><th> number-local-reference-configs <th> number-local-reference-configs <th> used
 * <tr><th> SHS2-2-body-reference-configs <th> SHS2-2-body-reference-configs <th> used
 * <tr><th> SHS3-3-body-reference-configs <th> SHS3-3-body-reference-configs <th> used
 * <tr><th> scale-type-total-energy <th> ML_ISCALE_TOTEN <th> used
 * <tr><th> average-energy-per-atom <th> average-energy-per-atom <th> used
 * <tr><th> regression-coeff <th> regression-coeff <th> used
 * <tr><th> sigma-v <th> sigv <th> not used
 * <tr><th> sigma-w <th> sigw <th> not used
 * <tr><th> variance-training-energies <th> variance-training-energies <th> used
 * <tr><th> variance-training-force <th> variance-training-force <th> used
 * <tr><th> variance-training-stress <th> variance-training-stress <th> used
 * </table>
 *******************************************************************************************/
class TagTranslator
{
  public:
    TagTranslator(void);
    /**
     * function returns the INCAR tag when given code tag
     *
     *@param tag input tag which is converted to INCAR tag from code tag
     */
    std::string get_InputTag(const std::string& tag) const;
    /**
     * function returns the code tag when given INCAR tag
     *
     *@param tag input tag which is converted to code tag from INCAR tag
     */
    std::string get_CodeTag(const std::string& tag) const;

  private:
    /**
     *map to INCAR tag. stores map from code tag INCAR tag
     */
    std::map<std::string, std::string> mapToInputTags;
    /**
     *map to Code tag. stores map from INCAR tag to code tag
     */
    std::map<std::string, std::string> mapToCodeTag;
};

} // namespace vaspml
#endif
