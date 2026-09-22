#include "utils.hpp"

#include <cstdarg>
#include <iostream>
#include <iterator>
#include <stdexcept>

namespace vaspml
{

std::string str(const char* format, ...)
{
    char buffer[STR_MAX];

    va_list args;
    va_start(args, format);
    int n = vsnprintf(buffer, STR_MAX, format, args);
    if (n < 0)
    {
        std::cerr << "WARNING: There was an encoding error, maybe some output was lost "
                     "(format string: \""
                  << format << "\").\n";
    }
    // n is the number of characters which should have been written (without null terminator). If n
    // is STR_MAX, already one character is lost.
    else if (static_cast<std::size_t>(n) >= STR_MAX)
    {
        std::cerr << "WARNING: Write buffer was too small, " << std::to_string(n - STR_MAX + 1)
                  << " characters were lost.\n";
    }
    va_end(args);

    // Construct resulting string from buffer contents.
    std::string result(buffer);

    return result;
}

} // namespace vaspml

namespace vaspml::string_tools
{

std::string ltrim(const String& s)
{
    std::size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const String& s)
{
    std::size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const String& s)
{
    return rtrim(ltrim(s));
}

void trimVector(Vec1String& stringVector)
{
    std::transform(stringVector.begin(),
                   stringVector.end(),
                   stringVector.begin(),
                   [](String& str) { return trim(str); });
}

std::string makeLowerCase(const String& in_str)
{
    std::string out_str = in_str;
    std::transform(out_str.begin(),
                   out_str.end(),
                   out_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return out_str;
}

std::string makeUpperCase(const String& in_str)
{
    std::string out_str = in_str;
    std::transform(out_str.begin(),
                   out_str.end(),
                   out_str.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return out_str;
}

bool checkForSubstring(const String& in_string, const String& sub_string)
{

    if (in_string.find(sub_string) != std::string::npos) return true;
    else return false;
}

Vec1String splitString(const String& inString, const String& delimiter, bool removeEmpty)
{
    Vec1String  parts;
    std::size_t startPos = 0;
    std::size_t endPos = 0;

    while ((endPos = inString.find(delimiter, startPos)) != std::string::npos)
    {
        parts.push_back(inString.substr(startPos, endPos - startPos));
        startPos = endPos + delimiter.length();
    }
    // Add the last token.
    parts.push_back(inString.substr(startPos));

    if (removeEmpty)
    {
        auto it = parts.begin();
        while (it != parts.end())
        {
            if (it->empty()) it = parts.erase(it);
            else ++it;
        }
    }

    return parts;
}

} // namespace vaspml::string_tools

namespace vaspml::vector_tools
{

std::vector<std::size_t> generateIntSequence(std::size_t start, std::size_t end)
{
    // assert( ( end <= start ) &&
    //         "ERROR in function VectorTools::generateIntSequence. end <= start " );
    std::vector<std::size_t> indexArray(end - start);
    std::iota(indexArray.begin(), indexArray.end(), start);
    return indexArray;
}

} // namespace vaspml::vector_tools

namespace vaspml::file_io
{

/// second argument declares whether it is a binary file or not
std::ifstream openFileI(const std::string fname, const std::ifstream::openmode file_mode)
{
    std::ifstream input_stream;

    input_stream.open(fname, file_mode);
    if (!input_stream)
    {
        std::cout << "Opening infile stream failed " << std::endl;
        std::cout << "Check your file name " << fname << std::endl;
        throw;
    }
    return input_stream;
}

std::ofstream openFileO(const std::string fname, const std::ofstream::openmode file_mode)
{
    std::ofstream output_stream;

    output_stream.open(fname, file_mode);
    if (!output_stream)
    {
        std::cout << "Opening outfile stream failed " << std::endl;
        std::cout << "Check your file name " << fname << std::endl;
        throw;
    }
    return output_stream;
}

} // namespace vaspml::file_io
