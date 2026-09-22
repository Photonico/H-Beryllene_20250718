#ifndef UTILS_HPP
#define UTILS_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include <iostream>

namespace vaspml
{

/// Buffer length for vsnprintf() call in str() function.
const std::size_t STR_MAX = 1024;
/***************************************************************************************************
 * Alternative to `printf()`, use format specifiers as usual but return string.
 *
 * Under the hood this uses `vsnprintf()` to avoid buffer overflows.
 *
 * @param format Format specifier forwarded to `vsnprintf()`.
 * @return String constructed from output of `vsnprintf()`.
 *
 * @warning As this function is calling `vsnprintf()` it inherits its problems: if the format
 * specifier does not match with actual arguments, undefined behavior may occur. Also an
 * intermediate fixed length buffer (size vaspml::STR_MAX) is used to construct the string result.
 * Hence, if the output exceeds the buffer, some characters may be lost.
 *
 * @note Originally it was intended to pass the format as a std::string. However, the NEC compiler
 * (5.0.1) does not give correct results (contains garbage strings) when `format` is passed as
 * std::string and `format.c_str()` is used as argument for `vsnprintf()`. Even when the contents of
 * `format` are copied first to a char array which is passed to `vsnprintf()` this does not help.
 *
 * @todo Test newer versions of the NEC compiler and eventually change the argument type back to
 * `std::string`.
 **************************************************************************************************/
std::string str(const char* format, ...);

} // namespace vaspml

namespace vaspml::string_tools
{

/// define a symbol for white-space
const std::string WHITESPACE = " \n\r\t\f\v";

/** left trim an input string
 * @param s -> input string to trim will not change
 */
std::string ltrim(const String& s);

/** right trim an input string
 * @param s -> input string to trim will not change
 */
std::string rtrim(const String& s);

/** apply left and right trim to remove trailing
 *  and ending white spaces
 * @param s -> input string to trim will not change
 */
std::string trim(const String& s);

/*******************************************************************************************
 * trimming the elements of string vector element wise.
 * @param stringVector the leading and trailing spaces of this vector will be removed
 *******************************************************************************************/
void trimVector(Vec1String& stringVector);

/** extract single value from string
 * @param data -> string containing single value to extract
 */
template<typename T>
T extractValue(const String& data)
{
    T value;
    if constexpr (std::is_same<T, int>::value) { value = std::stoi(trim(data).c_str()); }
    else if constexpr (std::is_same<T, float>::value) { value = std::stof(trim(data).c_str()); }
    else if constexpr (std::is_same<T, double>::value) { value = std::stod(trim(data).c_str()); }
    else if constexpr (std::is_same<T, std::string>::value) { value = trim(data).c_str(); }
    return value;
}

/** extract Vector from a string separated by spaces
 * @param data -> input string containing numbers separated by space
 */
template<typename T>
std::vector<T> extractData(const String& input)
{
    std::string       temp;
    std::vector<T>    data;
    std::stringstream strStream(input.c_str());
    while (strStream >> temp)
    {
        if constexpr (std::is_same<T, int>::value) { data.push_back(std::stoi(temp)); }
        else if constexpr (std::is_same<T, float>::value) { data.push_back(std::stof(temp)); }
        else if constexpr (std::is_same<T, double>::value) { data.push_back(std::stod(temp)); }
        else if constexpr (std::is_same<T, std::string>::value) { data.push_back(temp); }
    }
    return std::vector<T>(data.begin(), data.end());
}
/*
 * convert a string to lower case letters ( copy of string will be made )
 * @param in_str input string to be converted to lower case
 *
 */
std::string makeLowerCase(const String& in_str);
/*
 * convert a string to upper case letters ( copy of string will be made )
 * @param in_str input string to be converted to upper case case
 *
 */
std::string makeUpperCase(const String& in_str);
/*
 *
 * check if a string contains a substring ( function if case sensitive )
 * @param in_string string which will be searched to contain
 * @param sub_string the subsrting for which will be searched to be conatined
 * in in_string
 *
 */
bool checkForSubstring(const String& in_string, const String& sub_string);

/***************************************************************************************************
 * Split string at given delimeter string and return vector of string parts.
 *
 * @param inString Input string which will be split. Will not be destroyed.
 * @param delimiter String at which the input string will be split (can be also multiple
 *                  characters).
 * @param removeEmpty Removes empty string parts from the return vector at the end of the function.
 *
 * @note If `removeEmpty` is set to `false`, multiple consecutive delimiters may result in unwanted
 *       empty return vector entries. For example, splitting `a__b` (two underscores between
 *       `a` and `b`) with underscore as delimiter, results in a vector containing `"a"`, `""`
 *       and `"b"`!
 ***************************************************************************************************/
Vec1String splitString(const String& inString, const String& delimiter, bool removeEmpty = true);
} // namespace vaspml::string_tools

namespace vaspml::vector_tools
{
/** sort vector input_vector and return the indexes in sorted order
 *  input_vector will not be destroyed
 * @param input_vector  -->  vector containing numerical T to sort
 * @param sorted_indx -> contains the indices of the vector in sorted order
 */
template<typename T>
std::vector<std::size_t> argSort(const std::vector<T>& input_vector)
{
    std::vector<std::size_t> sorted_indx(input_vector.size());
    std::iota(sorted_indx.begin(), sorted_indx.end(), 0);

    // sort the sorted_indx array based on lambda funtcion
    // based on input_vector
    std::sort(sorted_indx.begin(),
              sorted_indx.end(),
              // def lambda function that defines the order
              [&input_vector](std::size_t indx1, std::size_t indx2) -> bool
              { return input_vector[indx1] < input_vector[indx2]; });
    return sorted_indx;
}
/**
 *check if two vectors are equal
 *
 *if the size of the vectors does not match they are said to be not equal
 *
 *@param vec1 first vector for comparison
 *@param vec2 second vector for comprison
 *
 *@note the tmeplate type T of the functions has to implement elementwise
 *comparison operator
 */
template<typename T>
bool checkEqual(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    if (vec1.size() != vec2.size()) { return false; }
    return std::equal(vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
}

template<typename T>
std::vector<std::size_t> argSortSlice(const std::vector<T>& input_vector,
                                      const std::size_t     start,
                                      const std::size_t     end)
{
    std::vector<std::size_t> sorted_indx(end - start);
    std::iota(sorted_indx.begin(), sorted_indx.end(), start);
    // sort the sorted_indx array based on lambda funtcion
    // based on input_vector
    std::sort(sorted_indx.begin(),
              sorted_indx.end(),
              // def lambda function that defines the order
              [&input_vector](std::size_t indx1, std::size_t indx2) -> bool
              { return input_vector[indx1] < input_vector[indx2]; });
    return sorted_indx;
}

/** return slice of a vector, note that the copy constructor will be called, can be expensive
 * @param array -> array from which a slice will be extracted
 * @param indx1 -> start indx from which slice starts
 * @param indx2 -> exclusive indx at which slice ends
 */
template<typename T>
std::vector<T> sliceVector(const std::vector<T>& array, std::size_t indx1, std::size_t indx2)
{
    return std::vector<T>(array.begin() + indx1, array.begin() + indx2);
}
/**
 * allocate a 1d array, all elements are set to zero
 * @param array -> input std::vector<T> to resize
 * @param Int n -> allocation size
 */
template<typename T>
void allocate_vector(std::vector<T>& array, const std::size_t n)
{
    array.resize(n);
}
/**
 * allocate a 2d array, all elements are set to zero
 * @param array  ->  input std::vector<std::vector<T>> to resize
 * @param Int n0 ->  allocation size for first dimension
 * @param Int n1 ->  allocation size for second dimension
 */
template<typename T>
void allocate_vector(std::vector<std::vector<T>>& array, const std::size_t n0, const std::size_t n1)
{
    array.resize(n0);
    for (std::size_t i = 0; i < n0; i++) { array[i].resize(n1); }
}
/**
 * allocate a 3d array, all elements are set to zero
 * @param array  ->  input std::vector<std::vector<T>> to resize
 * @param Int n0 ->  allocation size for first dimension
 * @param Int n1 ->  allocation size for second dimension
 * @param Int n2 ->  allocation size for second dimension
 */
template<typename T>
void allocate_vector(std::vector<std::vector<std::vector<T>>>& array,
                     const std::size_t                         n0,
                     const std::size_t                         n1,
                     const std::size_t                         n2)
{
    array.resize(n0);
    for (std::size_t i = 0; i < n0; i++)
    {
        array[i].resize(n1);
        for (std::size_t j = 0; j < n1; j++) { array[i][j].resize(n2); }
    }
}
/*******************************************************************************************
 * allocate a std::vector<std::vector<T>> with varaible lengths for every entry
 *
 * @param array will be resized to the desired lengths on output
 * @param sizes contains the array lengths for every entry in the first dimension of the array
 *
 * array[ 0 ].resize( sizes[0] ), array[1].resize( sizes[1] ) ....
 *******************************************************************************************/
template<typename T>
void allocate_vector(std::vector<std::vector<T>>& array, const Vec1Int& sizes)
{
    array.resize(sizes.size());
    for (std::size_t i = 0; i < sizes.size(); i++) { array[i].resize(sizes[i]); }
}
/**
 * generate a list of integer starting at start and ending at end
 * @param start start integer
 * @param end integer
 */
std::vector<std::size_t> generateIntSequence(std::size_t start, std::size_t end);

/**
 * generate a list of integer starting at start and ending at end
 * @param start start integer
 * @param end integer
 */
std::vector<std::size_t> generateIntSequence(std::size_t start, std::size_t end);

/*
 * get unique elements of vector without destruction
 *@param input_vector -> vector from which unique elements are taken
 *@return unique elements of vector
 */
template<typename T>
std::vector<T> get_unique(const std::vector<T>& input_vector)
{
    std::vector<T>                    output_vector = input_vector;
    typename std::vector<T>::iterator ip;
    std::sort(output_vector.begin(), output_vector.end());
    ip = std::unique(output_vector.begin(), output_vector.end());
    output_vector.resize(std::distance(output_vector.begin(), ip));
    return output_vector;
}

/*
 * get unique elements of vector with destruction of input vector
 *@param input_vector -> vector from which unique elements are taken
 *                       on output vector contains unique elements
 */
template<typename T>
void get_unique(std::vector<T>& vector)
{
    typename std::vector<T>::iterator ip;
    std::sort(vector.begin(), vector.end());
    ip = std::unique(vector.begin(), vector.end());
    vector.resize(std::distance(vector.begin(), ip));
}

/**
 *
 * count the number of occurences of
 *@param element in the std::vector
 *@param vec_in
 */
template<typename T>
T count_element(const std::vector<T>& vec_in, const T element)
{
    return count(vec_in.begin(), vec_in.end(), element);
}

/**
 * function return the maximum value in your std::vector
 *@param vec input vector from which to extract max value
 */
template<typename T>
T get_max(const std::vector<T>& vec)
{
    return *max_element(vec.begin(), vec.end());
}

/**
 * function return the maximum value in your std::vector
 *@param vec input vector from which to extract max value
 */
template<typename T>
T get_min(const std::vector<T>& vec)
{
    return *min_element(vec.begin(), vec.end());
}

/*******************************************************************************************
 * summing all elements of e vector
 * @param input vector which has to  be summed over
 *******************************************************************************************/
template<typename T>
T sum(const std::vector<T>& input)
{
    return std::accumulate(input.begin(), input.end(), (T)0);
}

/*******************************************************************************************
 * making a sum over a vector slice
 *
 * @param input input vector over which the sum is taken
 * @param start index of vector can at least be zero and maximal end -1
 * @param end can at least be start + 1 and maximimal input.size()-1
 *******************************************************************************************/
template<typename T>
T sum(const std::vector<T>& input, const std::size_t start, const std::size_t end)
{
    //VASPML_DEBUG_L1(
    //      if ( end < start ){
    //         tutor::bug( "ERROR: sum( const std::vector<T>& input, const std::size_t start, const
    //         std::size_t end )\n"
    //                     "       end index has to be larger than start index" );
    //      }
    //);
    return std::accumulate(input.begin() + start, input.end() + end, T::value_type(0));
}

/*******************************************************************************************
 * Makes the Hadamard product of two input vectors
 * @param vec1 input vector one which will be multiplied
 * @param vec2 input vector two which will be multiplied with vec1
 *
 * @return result \f[=vec1[i]*vec2[i]\f]
 *******************************************************************************************/
template<typename T>
std::vector<T> elementWiseProduct(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    //VASPML_DEBUG_L1(
    //      if ( vec1.size() != vec2.size() ){
    //         tutor::bug( "std::vector<T> elementWiseProduct( const std::vector<T>& vec1, const
    //         std::vector<T>& vec2 )\n"
    //                     "       size of vec1 and vec2 do not agree. vec1.size() != vec2.size()"
    //                     );
    //      }
    //);

    std::vector<T> result(vec1.size());
    std::transform(vec1.begin(),
                   vec1.end(),
                   vec2.begin(),
                   result.begin(),
                   [](const T& a, const T& b) { return a * b; });
    return result;
}

} // namespace vaspml::vector_tools

namespace vaspml::file_io
{

/** open file for input reading
 * @param fname -> string containing filename
 */
std::ifstream openFileI(const std::string             fname,
                        const std::ifstream::openmode file_mode = std::ifstream::in);
/** open file for input writing
 * @param fname -> string containing filename
 */
std::ofstream openFileO(const std::string             fname,
                        const std::ofstream::openmode file_mode = std::ifstream::out);

} // namespace vaspml::file_io

namespace vaspml::read_values
{

/// scalar reader
template<typename T, typename V>
void io_scalar(V& io_stream, T& scalar)
{

    // Do reading if input stream
    if constexpr (std::is_same<V, std::ifstream>::value)
    {

        if constexpr (std::is_integral<T>::value && std::is_same<T, bool>::value)
        {
            Int int_value;
            io_stream.read((char*)&int_value, sizeof(int_value));
            scalar = T(int_value);
        }
        else { io_stream.read((char*)&scalar, sizeof(scalar)); }
    }
    // Do writing
    else
    {
        if constexpr (std::is_integral<T>::value && std::is_same<T, bool>::value)
        {
            Int int_value = (Int)scalar;
            io_stream.write((char*)&int_value, sizeof(int_value));
        }
        else { io_stream.write((char*)&scalar, sizeof(scalar)); }
    }
}

/// array reader
template<typename T, typename V>
void io_array_binary(V& io_stream, std::vector<T>& array, const Int& number_entries)
{
    // Do reading if input stream
    if constexpr (std::is_same<V, std::ifstream>::value)
    {
        if constexpr (std::is_same<T, std::string>::value) // value here is member variable of
                                                           // is_same!!!
        {
            // Character reading
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                char value[3];

                io_stream.read((char*)&value, sizeof(value) - 1);
                value[2] = '\0';

                array.push_back(value);
            }
        }
        else
        {
            // Non-character reading
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                T value;

                io_stream.read((char*)&value, sizeof(value));
                array.push_back(value);
            }
        }
    }
    // Otherwise do writing
    else
    {
        if constexpr (std::is_same<T, std::string>::value) // value here is member variable of
                                                           // is_same!!!
        {
            // Character writing
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                io_stream.write(array[i_entry].c_str(), 2);
            }
        }
        else
        {
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                io_stream.write((char*)&array[i_entry], sizeof(array[i_entry]));
            }
        }
    }
}

} // namespace vaspml::read_values

#endif
