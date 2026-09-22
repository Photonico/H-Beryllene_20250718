#ifndef IO_HPP
#define IO_HPP

#include "Record.hpp"
#include "utils.hpp"

namespace vaspml::io
{

void checkLine(const String& line, const String& search, const Int n);

Int readLine(std::ifstream& input, std::string& buffer, Int skip = 0, bool expectEOF = false);

template<typename T>
void writeJsonScalarItem(std::ostream& output, T item)
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real  >) output << str("%24.16E", item);
    else if constexpr (std::is_same_v<T, Int   >) output << str("%d", item);
    else if constexpr (std::is_same_v<T, String>) output << '"' + item + '"';
    else if constexpr (std::is_same_v<T, bool  >) output << (item ? "true" : "false");
    else output << String{"null"};
    // clang-format on

    return;
}

void writeJson(const Record& record, std::ostream& output, Int level = 0, Int indentSpaces = 2);

template<typename T>
T readMlabScalar(String& line)
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real  >) return stod(line);
    else if constexpr (std::is_same_v<T, Int   >) return stoi(line);
    else if constexpr (std::is_same_v<T, String>) return string_tools::trim(line);
    // clang-format on
}

template<typename T>
T readMlabVector(std::ifstream& input, Int& n, bool expectEOF = false)
{
    String line;
    T      result;
    using scalarT = typename T::value_type;
    using namespace string_tools;
    while (readLine(input, line, 0, expectEOF) > 0)
    {
        // Increase line counter.
        n++;
        // Skip empty lines (actually this should not happen).
        if (line.size() < 1) continue;
        // If first character is one of the separator characters we reached the end of the vector.
        if (checkForSubstring("*-=", line.substr(0, 1))) break;
        // Split current line and transform elements into desired type.
        for (auto element : splitString(trim(line), " "))
        {
            result.push_back(readMlabScalar<scalarT>(element));
        }
    }

    return result;
}

void readMlab(Record& record, std::ifstream& input);

} // namespace vaspml::io

#endif
