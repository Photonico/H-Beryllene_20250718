#ifndef RECORD_HPP
#define RECORD_HPP

#include "types.hpp"

#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace vaspml
{

/*================================================================================================+
 |
 | SECTION Record, Item
 |
 +================================================================================================*/

// Forward declare Record, required for types and pointer in Item.
struct Record;

using ShRec = std::shared_ptr<Record>;
using ShCRec = std::shared_ptr<const Record>;
using Vec1ShRec = std::vector<ShRec>;
using Vec1ShCRec = std::vector<ShCRec>;

/**************************************************************************************************
 * Variant specialization which can represent all data types typically used in vaspml.
 *
 * For details, refer to Record.
 **************************************************************************************************/
using Item = std::variant<Real,
                          Int,
                          String,
                          bool,
                          ShRec,
                          Vec1Real,
                          Vec1Int,
                          Vec1String,
                          Vec1ShRec,
                          Vec2Real,
                          Vec2Int>;

enum class ItemIndex : int8_t
{
    // clang-format off
    REAL       = 1,
    INT        = 2,
    STRING     = 3,
    BOOL       = 4,
    SHREC      = 5,
    VEC1REAL   = 6,
    VEC1INT    = 7,
    VEC1STRING = 8,
    VEC1SHREC  = 9,
    VEC2REAL   = 10,
    VEC2INT    = 11
    // clang-format on
};

template<typename T>
constexpr std::false_type alwaysFalse{};
template<typename T>
constexpr ItemIndex itemIndex()
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real      >) return ItemIndex::REAL;
    else if constexpr (std::is_same_v<T, Int       >) return ItemIndex::INT;
    else if constexpr (std::is_same_v<T, String    >) return ItemIndex::STRING;
    else if constexpr (std::is_same_v<T, bool      >) return ItemIndex::BOOL;
    else if constexpr (std::is_same_v<T, ShRec     >) return ItemIndex::SHREC;
    else if constexpr (std::is_same_v<T, Vec1Real  >) return ItemIndex::VEC1REAL;
    else if constexpr (std::is_same_v<T, Vec1Int   >) return ItemIndex::VEC1INT;
    else if constexpr (std::is_same_v<T, Vec1String>) return ItemIndex::VEC1STRING;
    else if constexpr (std::is_same_v<T, Vec1ShRec >) return ItemIndex::VEC1SHREC;
    else if constexpr (std::is_same_v<T, Vec2Real  >) return ItemIndex::VEC2REAL;
    else if constexpr (std::is_same_v<T, Vec2Int   >) return ItemIndex::VEC2INT;
    else static_assert(alwaysFalse<T>, "Provided type is not an Index type.");
    // clang-format on
}

using RecordMap = std::unordered_map<String, Item>;

/// Metafunction to determine if a type is a std::vector of std::shared_ptr, general case (false).
template<typename T>
struct isVectorSharedPtr : std::false_type
{};
/// Specialization if type is actually a std::vector of std::shared_ptr (true).
template<typename T>
struct isVectorSharedPtr<std::vector<std::shared_ptr<T>>> : std::true_type
{};

/// Metafunction to determine if a type is a std::vector, general case (false).
template<typename T>
struct isVector : std::false_type
{};
/// Specialization if type is actually a std::vector (true).
template<typename T>
struct isVector<std::vector<T>> : std::true_type
{};

// clang-format off
/// Metafunction to determine if a type is a scalar type in the Item variant, general case (false).
template<typename T, bool scalar =
    std::is_same_v<T, Real>   ||
    std::is_same_v<T, Int>    ||
    std::is_same_v<T, String> ||
    std::is_same_v<T, bool>
>
// clang-format on
struct isScalar : std::false_type
{};
/// Specialization if type is actually a scalar type in the Item variant (true).
template<typename T>
struct isScalar<T, true> : std::true_type
{};

/**************************************************************************************************
 * Generic data structure used to store most data of vaspml.
 *
 * A Record contains a dictionary from string keys to vaspml::Item which is a specialization of
 * `std::variant`. Item contains only the common types we use in vaspml (e.g. Vec1Real, Vec2Int).
 * Two special entries are ShRec and Vec1ShRec, i.e., a single or array of smart pointers to another
 * instance of a Record. In this way records can be nested. Basically this data structure can be
 * thought of an in-memory representation of a JSON file (with limitations regarding array nesting).
 *
 * The Record class adds to this dictionary functionality to safely read and write access the data.
 * To add a new Record entry use the bracket notation equal to `std::map`:
 *
 * ~~~{.cpp}
 * Record myRecord;
 * myRecord["description"] = String{"This is my first Record entry"};
 * myRecord["age"] = 99;
 * myRecord["size"] = 2.00;
 * ~~~
 *
 * @warning Always wrap literal strings with `String{...}`, otherwise the compiler will try to
 * convert them to integer or bool type! It is best practice to be as explicit as possible to avoid
 * automatic type determination issues. For example, use `2.0` or `Real(2)`, otherwise `2` will be
 * stored in the Item as integer.
 *
 * Depending on
 * the underlying type use the following functions to retrieve data:
 *
 * Type        |  get()  |  cget() |  dget() | dcget() | vcget() |
 * ------------|:-------:|:-------:|:-------:|:-------:|:-------:|
 * Real        |    x    |    x    |         |         |         |
 * Int         |    x    |    x    |         |         |         |
 * String      |    x    |    x    |         |         |         |
 * bool        |    x    |    x    |         |         |         |
 * ShRec       |    x    |    x    |    x    |    x    |         |
 * Vec1Real    |    x    |    x    |         |         |         |
 * Vec1Int     |    x    |    x    |         |         |         |
 * Vec1String  |    x    |    x    |         |         |         |
 * Vec1ShRec   |    x    |         |         |         |    x    |
 * Vec2Real    |    x    |    x    |         |         |         |
 * Vec2Int     |    x    |    x    |         |         |         |
 *
 * Here, the variations of the get() functions have the following meaning:
 *
 * * get() ... read/write access to data,
 * * cget() ... read-only access to data,
 * * dget() ... read/write access to dereferenced data (for `shared_ptr` content),
 * * dcget() ... read-only access to dereferenced data (for `shared_ptr` content).
 * * vcget() ... special operation for Vec1ShRec, returns vector with pointers to read-only records.
 *
 * @attention Avoid using the `std::get<>()` function directly on the dictionary values because this
 * allows to accidentially modify data which is meant to be constant.
 *
 * **Examples:**
 *
 * ~~~{.cpp}
 * // Create a record intended to hold a data set of structures.
 * Record dataset;
 *
 * // Add basic data set information.
 * dataset["record"] = String("dataset");
 * dataset["numStructures"] = 2;
 * dataset["types"] = Vec1String{"Fe", "O"};
 *
 * // Create a vector of "sub"-records for the structures.
 * dataset["structures"] = Vec1ShRec(2, std::make_shared<Record>());
 *
 * // Fill the first structure with some data.
 * Record& s1 = *dataset.get<Vec1ShRec>("structures")[0];
 * s1["lattice"] = Vec2Real{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
 * s1["positions"] = Vec1Real{0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5};
 *
 * // Fill the second structure with some data.
 * Record& s2 = *dataset.get<Vec1ShRec>("structures")[1];
 * s2["lattice"] = Vec2Real{{1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
 * s2["positions"] = Vec1Real{0.5, 0.5, 0.3, 0.3, 0.0, 0.3, 0.0, 0.7, 0.1};
 *
 * // Access data and print to screen.
 * const Int& n = dataset.cget<Int>("numStructures");
 * //    Int& m = dataset.cget<Int>("numStructures"); // Does not compile, cget() returns const!
 * std::cout << "Record contains " << n << " structures with elements \"";
 * for (auto e : dataset.cget<Vec1String>("types")) std::cout << e << " ";
 * std::cout << "\": \n";
 * Vec1ShCRec structures = dataset.vcget<Vec1ShRec>("structures"); // Read-only access.
 * for (Int i = 0; i < n; ++i)
 * {
 *     std::cout << "Structure " << i << ": positions = ";
 *     for (auto p : structures[i]->cget<Vec1Real>("positions")) std::cout << p << " ";
 *     std::cout << "\n";
 * }
 *
 * // Alternatively, just write the whole Record in JSON format.
 * io::writeJson(dataset, std::cout);
 * ~~~
 *
 * @note The getter function get() and dget() by design cannot be called for instances of `const
 * Record`, only "const"-Versions (cget(), dcget(), vcget()) can be used in this case. This helps
 * enforcing const-correctness.
 **************************************************************************************************/
struct Record
{
    /// Dictionary holing all data in this Record.
    RecordMap dict;

    /// Forwards bracket operator of dictionary to Record.
    Item& operator[](String key);
    /// Forwards at() member function of dictionary to Record.
    const Item& at(String key) const;
    /// Check if key is present in Record.
    bool contains(String key) const;
    /// Return list of Record keys.
    Vec1String keys() const;
    /// Wrapper around std::get<>().
    template<typename T>
    T& get(String key)
    {
        return std::get<T>(dict[key]);
    }
    /// For `shared_ptr`-type content only, return reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto& dget(String key)
    {
        return *std::get<T>(dict[key]);
    }
    /// For non-`shared_ptr` and non-`vector<shared_ptr>` content only, return const reference to
    /// underlying data.
    template<typename T,
             typename std::enable_if_t<std::conjunction_v<std::negation<isSharedPtr<T>>,
                                                          std::negation<isVectorSharedPtr<T>>>>* =
                 nullptr>
    T const& cget(String key) const
    {
        return std::get<T>(dict.at(key));
    }
    /**********************************************************************************************
     * For `shared_ptr`-type content only, return `shared_ptr` to const version of underlying data.
     *
     * @note Creates a new `shared_ptr`, hence increases `use_count()`. Prefer dcget<>() if
     * possible.
     **********************************************************************************************/
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto cget(String key) const
    {
        using constT = const typename std::pointer_traits<T>::element_type;
        return std::shared_ptr<constT>(std::get<T>(dict.at(key)));
    }
    /// For `shared_ptr`-type content only, return const reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    const typename std::pointer_traits<T>::element_type& dcget(String key) const
    {
        return *std::get<T>(dict.at(key));
    }
    /**********************************************************************************************
     * Special getter function for Vec1ShRec Record entries.
     *
     * This getter is used to obtain a "const"-version of an Vec1ShRec Record entry. Background: a
     * naive version of such a getter could return a "const Vec1ShRec". However, this would still
     * allow to change the Record behind the shared pointers which we do NOT want! Hence, this
     * function instead creates a new vector of shared pointers to "const"-casted Record entries.
     *
     * @return Vector of shared pointers to const Records.
     **********************************************************************************************/
    template<typename T, typename std::enable_if_t<isVectorSharedPtr<T>::value>* = nullptr>
    auto vcget(String key) const
    {
        if (!std::holds_alternative<Vec1ShRec>(dict.at(key)))
        {
            throw std::runtime_error("ERROR: Item mapped to dictionary entry \"" + key
                                     + "\" does not hold a vector of shared pointers,"
                                     + "cannot execute vcget() function.");
        }
        // Unpack types, e.g. from Vec1ShRec extract Record.
        using sharedPtrT = typename T::value_type;
        using baseT = typename std::pointer_traits<sharedPtrT>::element_type;
        // Version of template type which is const at innermost level, e.g. ShCRec.
        using constT = std::shared_ptr<const baseT>;
        T const&            in = std::get<T>(dict.at(key));
        std::vector<constT> out(in.size());
        std::copy(in.begin(), in.end(), out.begin());
        return out;
    }
};

} // namespace vaspml

namespace vaspml::data
{

/***************************************************************************************************
 * Check if a given Record is of the correct "type".
 *
 * @param record Record to check (const, will not be modified).
 * @param desiredType The desired "type" of the Record.
 * @param func The name of the function from where this is called (usually "__func__").
 **************************************************************************************************/
void checkRecordType(const Record& record, String desiredType, const char* func);

/*================================================================================================+
 |
 | SECTION RecordFunctor
 |
 +================================================================================================*/

/// Base class for all functors used to traverse() a Record.
struct RecordFunctor
{
    Int    level = 0;
    String key;

    // clang-format off
    virtual void operator()(const Record& record, bool entry);
    virtual void operator()(const String& key, const Real&       item);
    virtual void operator()(const String& key, const Int&        item);
    virtual void operator()(const String& key, const String&     item);
    virtual void operator()(const String& key, const bool&       item);
    virtual void operator()(const String& key, const Vec1Real&   item);
    virtual void operator()(const String& key, const Vec1Int&    item);
    virtual void operator()(const String& key, const Vec1String& item);
    virtual void operator()(const String& key, const Vec2Real&   item);
    virtual void operator()(const String& key, const Vec2Int&    item);
    // clang-format on

    virtual void down(const String& key, Int arrayIndex = -1);
    virtual void up();
};

/**************************************************************************************************
 * Traverse a Record and apply functor to all Items.
 *
 * @param record Record to traverse (const, will not be modified).
 * @param f A functor which applies its operator() to all Items of the traversed Record.
 * @param skipRecordArrays If true, only traverse the first entry of each Vec1ShRec Item.
 *
 * This function recursively traverses all key-Item pairs of a given Record and applies the
 * functor's operator() to all Items.
 **************************************************************************************************/
void traverse(const Record& record, RecordFunctor& f, bool skipRecordArrays = false);

/*================================================================================================+
 |
 | SECTION PackSize
 |
 +================================================================================================*/

struct PackSize : public RecordFunctor
{
    Int packSize = 0;

    template<typename T>
    void add(const String& /* key */, const T& /* item */)
    {
        if constexpr (isVector<T>::value)
        {
            using T1 = typename T::value_type;
            if constexpr (isVector<T1>::value)
            {
                //using T2 = typename T1::value_type;
            }
        }

        return;
    }

    // clang-format off
    void operator()(const Record& record, bool entry) override;
    void operator()(const String& key, const Int&        item) override;
    void operator()(const String& key, const Real&       item) override;
    void operator()(const String& key, const String&     item) override;
    void operator()(const String& key, const bool&       item) override;
    void operator()(const String& key, const Vec1Real&   item) override;
    void operator()(const String& key, const Vec1Int&    item) override;
    void operator()(const String& key, const Vec1String& item) override;
    void operator()(const String& key, const Vec2Real&   item) override;
    void operator()(const String& key, const Vec2Int&    item) override;
    // clang-format on
};

Int M_packSize(const Record& record);

/*================================================================================================+
 |
 | SECTION MemoryCounter
 |
 +================================================================================================*/

std::size_t memoryUsageMapEntry(const String& input);

struct MemoryCounter : public RecordFunctor
{
    Record memInfo;

    MemoryCounter();

    template<typename T>
    void process(const String& key, const T& item)
    {
        const Real overhead = memoryUsageMapEntry(key);
        Real       data = sizeof(T);

        if constexpr (isVector<T>::value)
        {
            using T1 = typename T::value_type;
            data += item.size() * sizeof(T1);
            if constexpr (isVector<T1>::value)
            {
                using T2 = typename T1::value_type;
                for (const auto& i : item) data += i.size() * sizeof(T2);
            }
        }

        memInfo.get<Real>("overhead") += overhead;
        memInfo.get<Real>("data") += data;

        return;
    }

    // clang-format off
    void operator()(const Record& record, bool entry) override;
    void operator()(const String& key, const Int&        item) override;
    void operator()(const String& key, const Real&       item) override;
    void operator()(const String& key, const String&     item) override;
    void operator()(const String& key, const bool&       item) override;
    void operator()(const String& key, const Vec1Real&   item) override;
    void operator()(const String& key, const Vec1Int&    item) override;
    void operator()(const String& key, const Vec1String& item) override;
    void operator()(const String& key, const Vec2Real&   item) override;
    void operator()(const String& key, const Vec2Int&    item) override;
    // clang-format on
};

std::size_t memoryUsageString(const String& input);

String printMemoryUsage(const Record& memInfo);

String printRecordByteSizes();

} // namespace vaspml::data

#endif
