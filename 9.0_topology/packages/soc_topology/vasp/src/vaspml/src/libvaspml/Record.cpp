#include "Record.hpp"
#include "SmartEnum.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>

using namespace vaspml;
using namespace vaspml::data;

/*================================================================================================+
 |
 | SECTION Record, Item
 |
 +================================================================================================*/

template<>
SmartEnum<ItemIndex>::EnumMap SmartEnum<ItemIndex>::names = {
    {ItemIndex::REAL,       "Real"      },
    {ItemIndex::INT,        "Int"       },
    {ItemIndex::STRING,     "String"    },
    {ItemIndex::BOOL,       "bool"      },
    {ItemIndex::SHREC,      "ShRec"     },
    {ItemIndex::VEC1REAL,   "Vec1Real"  },
    {ItemIndex::VEC1INT,    "Vec1Int"   },
    {ItemIndex::VEC1STRING, "Vec1String"},
    {ItemIndex::VEC1SHREC,  "Vec1ShRec" },
    {ItemIndex::VEC2REAL,   "Vec2Real"  },
    {ItemIndex::VEC2INT,    "Vec2Int"   },
};

Item& Record::operator[](String key)
{
    return dict[key];
}

const Item& Record::at(String key) const
{
    return dict.at(key);
}

bool Record::contains(String key) const
{
    return dict.find(key) != dict.end();
}

Vec1String Record::keys() const
{
    Vec1String vkeys(dict.size());
    std::transform(dict.begin(),
                   dict.end(),
                   vkeys.begin(),
                   [](auto const& entry) { return entry.first; });
    return vkeys;
}

void vaspml::data::checkRecordType(const Record& record, String desiredType, const char* func)
{
    if (!record.contains("record"))
    {
        throw std::runtime_error("ERROR: Record used in function \"" + String{func}
                                 + "\" is not typed (i.e. does not have a \"record\" field).");
    }
    const String& actualType = record.cget<String>("record");
    if (actualType != desiredType)
    {
        throw std::runtime_error("ERROR: Record used in function \"" + String{func}
                                 + "\" does not have desired type \"" + desiredType + "\", got \""
                                 + actualType + "\" instead.");
    }

    return;
}

/*================================================================================================+
 |
 | SECTION RecordFunctor
 |
 +================================================================================================*/

// clang-format off
void RecordFunctor::operator()(const Record& /* record */, bool /* entry */) {}
void RecordFunctor::operator()(const String& /* key */, const Real&       /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Int&        /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const String&     /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const bool&       /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1Real&   /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1Int&    /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1String& /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec2Real&   /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec2Int&    /* item */) {}
// clang-format on

void RecordFunctor::down(const String& key, Int /* arrayIndex */)
{
    level++;
    this->key = key;
    //std::cout << "Record \"" << key << "\"";
    //if (arrayIndex >= 0) std::cout << "[" << arrayIndex << "]";
    //std::cout << " (level " << level << ")\n";

    return;
}

void RecordFunctor::up()
{
    level--;
    //std::cout << "End Record \"" << key << "\"\n";

    return;
}

void vaspml::data::traverse(const Record& record, RecordFunctor& f, bool skipRecordArrays)
{
    f(record, true);
    const Vec1String keys = record.keys();
    for (auto k = keys.cbegin(); k != keys.cend(); ++k)
    {
        std::visit(
            [&](auto&& arg)
            {
                using T = std::decay_t<decltype(arg)>;
                //******************
                // Item types: ShRec
                //******************
                if constexpr (std::is_same_v<T, ShRec>)
                {
                    f.down(*k);
                    traverse(*arg, f);
                    f.up();
                }
                //**********************
                // Item types: Vec1ShRec
                //**********************
                else if constexpr (std::is_same_v<T, Vec1ShRec>)
                {
                    auto                      begin = arg.cbegin();
                    Vec1ShRec::const_iterator end;
                    if (skipRecordArrays) end = arg.cbegin() + 1;
                    else end = arg.cend();
                    for (auto a1 = begin; a1 != end; ++a1)
                    {
                        f.down(*k, std::distance(begin, a1));
                        traverse(**a1, f);
                        f.up();
                    }
                }
                //******************
                // Actual data types
                //******************
                else f(*k, arg);
            },
            record.at(*k));
    }
    f(record, false);

    return;
}

/*================================================================================================+
 |
 | SECTION PackSize
 |
 +================================================================================================*/

void PackSize::operator()(const Record& /* record */, bool /* entry */)
{

    return;
}

// clang-format off
void PackSize::operator()(const String& key, const Int&        item) { add(key, item); }
void PackSize::operator()(const String& key, const Real&       item) { add(key, item); }
void PackSize::operator()(const String& key, const String&     item) { add(key, item); }
void PackSize::operator()(const String& key, const bool&       item) { add(key, item); }
void PackSize::operator()(const String& key, const Vec1Real&   item) { add(key, item); }
void PackSize::operator()(const String& key, const Vec1Int&    item) { add(key, item); }
void PackSize::operator()(const String& key, const Vec1String& item) { add(key, item); }
void PackSize::operator()(const String& key, const Vec2Real&   item) { add(key, item); }
void PackSize::operator()(const String& key, const Vec2Int&    item) { add(key, item); }
// clang-format on

Int vaspml::data::M_packSize(const Record& record)
{
    PackSize f;
    traverse(record, f);

    return f.packSize;
}

/*================================================================================================+
 |
 | SECTION MemoryCounter
 |
 +================================================================================================*/

std::size_t vaspml::data::memoryUsageString(const String& input)
{
    // Check if input string is smaller than SSO buffer. If yes, only the the stack memory is
    // required. If no, add stack and heap memory together.
    if (input.length() <= String{}.capacity()) return sizeof(String);
    else return sizeof(String) + input.length() * sizeof(String::value_type);
}

std::size_t vaspml::data::memoryUsageMapEntry(const String& input)
{
    constexpr std::size_t mem = sizeof(Item)       // Stack size of std::variant specialization.
                              + 2 * sizeof(Item*); // Pointers for double-linked list (assumption).

    return mem + memoryUsageString(input);
}

MemoryCounter::MemoryCounter()
{
    memInfo["record"] = String{"memInfo"};
    memInfo["overhead"] = Real(0.0);
    memInfo["data"] = Real(0.0);
}

void MemoryCounter::operator()(const Record& record, bool entry)
{
    if (entry) memInfo.get<Real>("overhead") += sizeof(record);

    return;
}

// clang-format off
void MemoryCounter::operator()(const String& key, const Int&        item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Real&       item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const String&     item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const bool&       item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1Real&   item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1Int&    item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1String& item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec2Real&   item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec2Int&    item) { process(key, item); }
// clang-format on

String vaspml::data::printMemoryUsage(const Record& memInfo)
{
    String result;

    checkRecordType(memInfo, "memInfo", __func__);

    Vec1String keys = {"Memory group", "overhead", "data", "total"};
    Real       total = 0.0;
    // Compute format specifier for key column.
    const String keyFormat = str("%%-%zus",
                                 std::max_element(keys.begin(),
                                                  keys.end(),
                                                  [](const String& a, const String& b)
                                                  { return a.length() < b.length(); })
                                     ->length());

    // Prepare table separator line.
    const String separator = str(("|-" + keyFormat + "-|-%16s-|-%12s-|-%8s-|\n").c_str(),
                                 String(keys.at(0).length(), '-').c_str(),
                                 String(16, '-').c_str(),
                                 String(12, '-').c_str(),
                                 String(8, '-').c_str());

    // Table header.
    result += "\n Memory usage estimation\n";
    result += " -----------------------\n";
    result += str(("| " + keyFormat + " | %16s | %12s | %8s |\n").c_str(),
                  keys.at(0).c_str(),
                  "kB",
                  "MB",
                  "GB");
    keys.erase(keys.begin());
    result += separator;

    // Fill table with all key-value pairs.
    for (auto k : keys)
    {
        Real mem;
        if (k == "total")
        {
            mem = total;
            result += separator;
        }
        else
        {
            mem = memInfo.cget<Real>(k);
            total += mem;
        }
        result += str(("| " + keyFormat + " | %16.3lf | %12.2lf | %8.1lf |\n").c_str(),
                      k.c_str(),
                      mem / 1.0E3,
                      mem / 1.0E6,
                      mem / 1.0E9);
    }
    result += "\n";

    return result;
}

String vaspml::data::printRecordByteSizes()
{
    String result;

    result += str("sizeof(Record)     = %zu\n", sizeof(Record));
    result += str("sizeof(RecordMap)  = %zu\n\n", sizeof(RecordMap));

    result += str("sizeof(Item)       = %zu\n", sizeof(Item));
    result += str("sizeof(Real)       = %zu\n", sizeof(Real));
    result += str("sizeof(Int)        = %zu\n", sizeof(Int));
    result += str("sizeof(String)     = %zu\n", sizeof(String));
    result += str("sizeof(bool)       = %zu\n", sizeof(bool));
    result += str("sizeof(ShRec)      = %zu\n", sizeof(ShRec));
    result += str("sizeof(Vec1Real)   = %zu\n", sizeof(Vec1Real));
    result += str("sizeof(Vec1Int)    = %zu\n", sizeof(Vec1Int));
    result += str("sizeof(Vec1String) = %zu\n", sizeof(Vec1String));
    result += str("sizeof(Vec1ShRec)  = %zu\n", sizeof(Vec1ShRec));
    result += str("sizeof(Vec2Real)   = %zu\n", sizeof(Vec2Real));
    result += str("sizeof(Vec2Int)    = %zu\n", sizeof(Vec2Int));

    return result;
}
