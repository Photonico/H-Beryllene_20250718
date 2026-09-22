#ifndef SMARTENUM_HPP
#define SMARTENUM_HPP

#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace vaspml
{

/*================================================================================================+
 | NOTE: Steps to add new SmartEnums here:
 | 1. Add forward declaration right below.
 | 2. Add enum to metafunction isSmart below.
 | 3. Add declaration of static "names" member at the end of this file.
 +================================================================================================*/
// Forward declare all enums that should be treated as smart enum.
// Add aditional enums here:
namespace math
{
enum class CutoffType;
}
enum class ItemIndex : int8_t;
enum class DescriptorStorage;
enum class TotalEnergyType;
// Define metafunction to check whether a given enum is a smart enum.
// Only if listed here the smart enum features work.
// First, the general case that defines the negative outcome.
// Add additional enums here with || (OR) combination:
// clang-format off
template<typename T, bool smart =
    std::is_same_v<T, math::CutoffType> ||
    std::is_same_v<T, ItemIndex> ||
    std::is_same_v<T, DescriptorStorage> ||
    std::is_same_v<T, TotalEnergyType>
>
// clang-format on
struct isSmartEnum : std::false_type
{};
// Now, the partial specialization if the above expression is true.
template<typename T>
struct isSmartEnum<T, true> : std::true_type
{};

/** Provides static functions to convert enums to strings and vice-versa.
 *
 * Strings corresponding to enums need to be defined first in a map, see implementation for
 * math::CutoffType. In addition, the metafunction isSmartEnum needs to be updated (requires a
 * forward declaration in the header of SmartEnum. Also, an explicit specialization of the "names"
 * static member needs to be added to the header.
 */
template<typename T>
class SmartEnum
{
    static_assert(isSmartEnum<T>::value, "Template argument is not a smart enum.");

  public:
    static T toEnum(std::string stringToConvert)
    {
        for (const auto& n : names)
        {
            if (n.second == stringToConvert) return n.first;
        }
        throw std::runtime_error("ERROR: Provided string \"" + stringToConvert
                                 + "\" does not correspond to an enum.");
    };

    static std::string toString(T enumToConvert) { return names.at(enumToConvert); };

    static std::vector<T> listEnums()
    {
        std::vector<T> keys;
        for (const auto& n : names) keys.push_back(n.first);
        return keys;
    };

    static std::vector<std::string> listNames()
    {
        std::vector<std::string> values;
        for (const auto& n : names) values.push_back(n.second);
        return values;
    };

  private:
    typedef std::map<T, std::string> EnumMap;
    static EnumMap                   names;
};

template<typename T, typename std::enable_if_t<isSmartEnum<T>::value>* = nullptr>
std::ostream& operator<<(std::ostream& os, const T& enumToStream)
{
    os << SmartEnum<T>::toString(enumToStream);
    return os;
}

template<typename T, typename std::enable_if_t<isSmartEnum<T>::value>* = nullptr>
std::string toString(const T& enumToConvert)
{
    std::ostringstream ss;
    ss << SmartEnum<T>::toString(enumToConvert);
    return ss.str();
}

// Explicit specialization of static data member "names".
template<>
SmartEnum<math::CutoffType>::EnumMap SmartEnum<math::CutoffType>::names;
template<>
SmartEnum<ItemIndex>::EnumMap SmartEnum<ItemIndex>::names;
template<>
SmartEnum<DescriptorStorage>::EnumMap SmartEnum<DescriptorStorage>::names;
template<>
SmartEnum<TotalEnergyType>::EnumMap SmartEnum<TotalEnergyType>::names;

} //namespace vaspml

#endif
