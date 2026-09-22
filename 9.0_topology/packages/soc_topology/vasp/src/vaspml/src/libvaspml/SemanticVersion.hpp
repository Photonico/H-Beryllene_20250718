#ifndef SEMANTICVERSION_HPP
#define SEMANTICVERSION_HPP

#include <string>

namespace vaspml
{

/***************************************************************************************************
 * Implementation of semantic versioning rules (version 2.0.0).
 *
 * For detailed rules description see https://semver.org/.
 * @note This implementation does not recognize dot-separated identifiers of pre-release and build
 * metadata strings as separate tokens. Hence, for comparison only the full strings are used.
 **************************************************************************************************/
class SemanticVersion
{
  public:
    SemanticVersion() = default;
    explicit SemanticVersion(int major, int minor = 0, int patch = 0, std::string suffix = "");
    explicit SemanticVersion(std::string version);
    SemanticVersion& update(int major, int minor = 0, int patch = 0, std::string suffix = "");
    SemanticVersion& update(std::string version);
    int              major() const;
    int              minor() const;
    int              patch() const;
    std::string      suffix() const;
    std::string      prerelease() const;
    std::string      build() const;
    std::string      kind() const;
    std::string      toString() const;
    bool             operator<(const SemanticVersion& rhs) const;
    bool             operator==(const SemanticVersion& rhs) const;
    bool             operator!=(const SemanticVersion& rhs) const;
    bool             operator>(const SemanticVersion& rhs) const;
    bool             operator<=(const SemanticVersion& rhs) const;
    bool             operator>=(const SemanticVersion& rhs) const;
    static void      parseVersionString(std::string  version,
                                        int&         major,
                                        int&         minor,
                                        int&         patch,
                                        std::string& suffix);

  private:
    int         major_ = 0;
    int         minor_ = 0;
    int         patch_ = 0;
    std::string suffix_ = "";
    std::string prerelease_ = "";
    std::string build_ = "";
    std::string kind_ = "release";

    void processSuffix();
    void checkValidity() const;
};

} //namespace vaspml

#endif
