#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace vaspml
{

// Basic types.
using Real = double;
using Int = int;
using Int64 = int64_t;
using UInt = unsigned int;
using UInt64 = uint64_t;
using String = std::string;
using Size_t = std::size_t;

// Vectors of basic types.
using Vec1Real = std::vector<Real>;
using Vec2Real = std::vector<std::vector<Real>>;
using Vec3Real = std::vector<std::vector<std::vector<Real>>>;

using Vec1Int = std::vector<Int>;
using Vec2Int = std::vector<std::vector<Int>>;
using Vec3Int = std::vector<std::vector<std::vector<Int>>>;

using Vec1Int64 = std::vector<Int64>;
using Vec2Int64 = std::vector<std::vector<Int64>>;
using Vec3Int64 = std::vector<std::vector<std::vector<Int64>>>;

using Vec1UInt = std::vector<UInt>;
using Vec2UInt = std::vector<std::vector<UInt>>;
using Vec3UInt = std::vector<std::vector<std::vector<UInt>>>;

using Vec1UInt64 = std::vector<UInt64>;
using Vec2UInt64 = std::vector<std::vector<UInt64>>;
using Vec3UInt64 = std::vector<std::vector<std::vector<UInt64>>>;

using Vec1Size_t = std::vector<size_t>;
using Vec2Size_t = std::vector<std::vector<size_t>>;
using Vec3Size_t = std::vector<std::vector<std::vector<size_t>>>;

using Vec1String = std::vector<String>;
using Vec2String = std::vector<std::vector<String>>;
using Vec3String = std::vector<std::vector<std::vector<String>>>;

// Shared pointers to vectors of basic types.
using ShVec1Real = std::shared_ptr<Vec1Real>;
using ShVec2Real = std::shared_ptr<Vec2Real>;
using ShVec3Real = std::shared_ptr<Vec3Real>;

using ShVec1Int = std::shared_ptr<Vec1Int>;
using ShVec2Int = std::shared_ptr<Vec2Int>;
using ShVec3Int = std::shared_ptr<Vec3Int>;

using ShVec1Int64 = std::shared_ptr<Vec1Int64>;
using ShVec2Int64 = std::shared_ptr<Vec2Int64>;
using ShVec3Int64 = std::shared_ptr<Vec3Int64>;

using ShVec1UInt = std::shared_ptr<Vec1UInt>;
using ShVec2UInt = std::shared_ptr<Vec2UInt>;
using ShVec3UInt = std::shared_ptr<Vec3UInt>;

using ShVec1UInt64 = std::shared_ptr<Vec1UInt64>;
using ShVec2UInt64 = std::shared_ptr<Vec2UInt64>;
using ShVec3UInt64 = std::shared_ptr<Vec3UInt64>;

using ShVec1String = std::shared_ptr<Vec1String>;
using ShVec2String = std::shared_ptr<Vec2String>;
using ShVec3String = std::shared_ptr<Vec3String>;

// Const content versions of shared pointers used in MultiType.
using ShCVec1Real = std::shared_ptr<const Vec1Real>;
using ShCVec1Int = std::shared_ptr<const Vec1Int>;
using ShCVec1String = std::shared_ptr<const Vec1String>;
using ShCVec2Real = std::shared_ptr<const Vec2Real>;
using ShCVec2Int = std::shared_ptr<const Vec2Int>;

/// Metafunction to determine if a type is a shared_ptr, general case (false).
template<typename T>
struct isSharedPtr : std::false_type
{};
/// Specialization if type is actually a shared_ptr (true).
template<typename T>
struct isSharedPtr<std::shared_ptr<T>> : std::true_type
{};

// clang-format off
/***************************************************************************************************
 * Variant type which can represent all data types typically used in vaspml.
 *
 * @attention Inside vaspml do not use `std::get<>()` to access data members directly because this
 * allows to accidentially modify data which is meant to be constant. Instead, use
 * one of the getter functions:
 * * get<>() ... read/write access to data,
 * * cget<>() ... read-only access to data,
 * * dget<>() ... read/write access to dereferenced data (for `shared_ptr` content),
 * * dcget<>() ... read-only access to dereferenced data (for `shared_ptr` content).
 *
 * For non-`shared_ptr`-type content only get<>() and cget<>() are available. Choose cget<>() if you
 * do not intend to change the value, e.g.
 *
 *     MultiType a = 3;
 *     ...
 *     const Real& b = a.cget<Real>();
 *     b = 4; // Compile error!
 *
 * For `shared_ptr`-type content all four getters are available, here are some typical use cases:
 *
 * - Sharing ownership with some larger function or a class instance which may modify the contents:
 *
 *       MultiType v = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3});
 *       ...
 *       ShVec1Int b = v.get<ShVec1Int>();
 *       b->push_back(4);
 *
 *   `b` increases the `use_count()` of the shared pointer and may be used to modify the underlying
 *   vector. Basically, get<>() just wraps around `std::get()`.
 *
 * - Pass on a non-modifiable (const) reference of the underlying vector for reading only:
 *
 *       MultiType v = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3});
 *       ...
 *       const Vec1Int& b = v.dcget<ShVec1Int>();
 *       std::cout << b[2]; // OK!
 *       b.push_back(4);    // Compile error!
 *
 *
 * If the MultiType instance is `const` (e.g. when passed as `const MultiType&`) only the getters
 * cget<>() and dcget<>() are available, i.e., the underlying data is read-only. The "const-ness" is
 * passed on to the vector for `shared_ptr`-type contents:
 *
 *
 *     MultiType v = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3});
 *     ...
 *     void func(const MultiType& varg)
 *     {
 *         Vec1Int& b = v.dget<ShVec1Int>();        // Compile error!
 *         const Vec1Int& c = v.dcget<ShVec1Int>(); // OK!
 *         ...
 *     }
 *
 **************************************************************************************************/
class MultiType : public std::variant<Real,
                                      Int,
                                      String,
                                      bool,
                                      ShVec1Real,
                                      ShVec1Int,
                                      ShVec1String,
                                      ShVec2Real,
                                      ShVec2Int>
// clang-format on
{
    using variant::variant;

  public:
    /// Wrapper around `std::get()`, data is read/write accessible, use with care!
    template<typename T>
    T& get()
    {
        return std::get<T>(*this);
    }
    /// For `shared_ptr`-type content only, return reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto& dget()
    {
        return *std::get<T>(*this);
    }
    /***********************************************************************************************
     * For `shared_ptr`-type content only, return `shared_ptr` to const version of underlying data.
     *
     * @note Creates a new `shared_ptr`, hence increases `use_count()`. Prefer dcget<>() if
     * possible.
     **********************************************************************************************/
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto cget() const
    {
        using constT = const typename std::pointer_traits<T>::element_type;
        return std::shared_ptr<constT>(std::get<T>(*this));
    }
    /// For non-`shared_ptr` content only, return const reference to underlying data.
    template<typename T, typename std::enable_if_t<std::negation_v<isSharedPtr<T>>>* = nullptr>
    auto& cget() const
    {
        return std::get<T>(*this);
    }
    /// For `shared_ptr`-type content only, return const reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    const typename std::pointer_traits<T>::element_type& dcget() const
    {
        return *std::get<T>(*this);
    }
};

/// A map from string to MultiType to store arbitrary data.
using MultiTypeMap = std::map<String, MultiType>;

/***************************************************************************************************
 * Helper function which either returns input or creates new shared pointer to object.
 *
 * @param in Input shared pointer or `nullptr`.
 *
 * If `in` is the null pointer a new object and a `shared_ptr` to it is created with
 * std::make_shared(). Otherwise, the input pointer is just returned. Use in constructor for classes
 * which accept their memory to be managed from outside, e.g.:
 *
 * @code
 * MyClass::MyClass(inputArray)
 * {
 *     this->array = vaspml::assignOrMakeShared(inputArray)
 * }
 * @endcode
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
T assignOrMakeShared(T in)
{
    if (in == nullptr) return std::make_shared<typename std::pointer_traits<T>::element_type>();
    else return in;
}

} //namespace vaspml

#endif
