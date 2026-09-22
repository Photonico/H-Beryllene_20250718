#include "SemanticVersion.hpp"

#include <cstddef>
#include <stdexcept>
//#include <charconv>
#include <string>

using namespace vaspml;

SemanticVersion::SemanticVersion(int major, int minor, int patch, std::string suffix)
{
    update(major, minor, patch, suffix);
}

SemanticVersion::SemanticVersion(std::string version)
{
    update(version);
}

SemanticVersion& SemanticVersion::update(int major, int minor, int patch, std::string suffix)
{
    major_ = major;
    minor_ = minor;
    patch_ = patch;
    suffix_ = suffix;

    processSuffix();
    checkValidity();

    return *this;
}

SemanticVersion& SemanticVersion::update(std::string version)
{
    parseVersionString(version, major_, minor_, patch_, suffix_);
    update(major_, minor_, patch_, suffix_);

    return *this;
}

std::string SemanticVersion::toString() const
{
    return std::to_string(major_) + '.' + std::to_string(minor_) + '.' + std::to_string(patch_)
         + suffix_;
}

int SemanticVersion::major() const
{
    return major_;
}

int SemanticVersion::minor() const
{
    return minor_;
}

int SemanticVersion::patch() const
{
    return patch_;
}

std::string SemanticVersion::suffix() const
{
    return suffix_;
}

std::string SemanticVersion::prerelease() const
{
    return prerelease_;
}

std::string SemanticVersion::build() const
{
    return build_;
}

std::string SemanticVersion::kind() const
{
    return kind_;
}

bool SemanticVersion::operator<(const SemanticVersion& rhs) const
{
    // clang-format off
    if      (major_ < rhs.major_) return true;
    else if (major_ > rhs.major_) return false;
    if      (minor_ < rhs.minor_) return true;
    else if (minor_ > rhs.minor_) return false;
    if      (patch_ < rhs.patch_) return true;
    else if (patch_ > rhs.patch_) return false;
    // clang-format on
    // An empty pre-release string indicates a higher version.
    if (prerelease_.empty() != rhs.prerelease_.empty())
    {
        if (prerelease_.empty()) return false;
        else return true;
    }
    // Otherwise, compare lexicographically.
    if (prerelease_ < rhs.prerelease_) return true;
    else if (prerelease_ > rhs.prerelease_) return false;

    return false;
}

bool SemanticVersion::operator==(const SemanticVersion& rhs) const
{
    if (major_ != rhs.major_) return false;
    if (minor_ != rhs.minor_) return false;
    if (patch_ != rhs.patch_) return false;
    if (prerelease_ != rhs.prerelease_) return false;

    return true;
}

bool SemanticVersion::operator!=(SemanticVersion const& rhs) const
{
    return !((*this) == rhs);
}
bool SemanticVersion::operator>(SemanticVersion const& rhs) const
{
    return rhs < (*this);
}
bool SemanticVersion::operator<=(SemanticVersion const& rhs) const
{
    return !((*this) > rhs);
}
bool SemanticVersion::operator>=(SemanticVersion const& rhs) const
{
    return !((*this) < rhs);
}

void SemanticVersion::parseVersionString(std::string  version,
                                         int&         major,
                                         int&         minor,
                                         int&         patch,
                                         std::string& suffix)
{
    // Initialize to defaults.
    major = 0;
    minor = 0;
    patch = 0;
    suffix = "";

    if (version.empty())
    {
        throw std::runtime_error("ERROR: Empty version string passed to SemanticVersion.");
    }

    // Major version number.
    //**********************
    std::size_t i = version.find('.');
    // Single number version, e.g. just "2", will convert to "2.0.0".
    if (i == std::string::npos)
    {
        //std::from_chars(version.data(), version.data() + version.size(), major);
        major = std::stoi(version);
        return;
    }
    // Only "." without a preceding number is invalid.
    else if (i < 1)
    {
        throw std::runtime_error("ERROR: Invalid version string (early major/minor separator)");
    }
    // Read major version number.
    //else std::from_chars(version.data(), version.data() + i, major);
    else major = std::stoi(version.substr(0, i));
    // If string ends after ".", finish parsing.
    if (i + 1 == version.size()) return;
    // Cut away already processed version string.
    version = version.substr(i + 1);

    // Minor version number.
    //**********************
    i = version.find('.');
    // Two-number version, e.g. just "2.1", will convert to "2.1.0".
    if (i == std::string::npos)
    {
        //std::from_chars(version.data(), version.data() + version.size(), minor);
        minor = std::stoi(version);
        return;
    }
    // Only "." without a preceding number is invalid.
    else if (i < 1)
    {
        throw std::runtime_error("ERROR: Invalid version string (early minor/patch separator)");
    }
    // Read minor version number.
    //else std::from_chars(version.data(), version.data() + i, minor);
    else minor = std::stoi(version.substr(0, i));
    // If string ends after ".", finish parsing.
    if (i + 1 == version.size()) return;
    // Cut away already processed version string.
    version = version.substr(i + 1);

    // Patch version number.
    //**********************
    i = std::min(version.find('-'), version.find('+'));
    // No suffix, e.g. just "2.1.3".
    if (i == std::string::npos)
    {
        //std::from_chars(version.data(), version.data() + version.size(), patch);
        patch = std::stoi(version);
        return;
    }
    // Read patch version number.
    //else std::from_chars(version.data(), version.data() + i, patch);
    else patch = std::stoi(version.substr(0, i));
    // Cut away already processed version string and store as suffix.
    suffix = version.substr(i);

    // Split suffix and use only first part if there are spaces.
    i = suffix.find(' ');
    if (i != std::string::npos) { suffix = suffix.substr(0, i); }

    return;
}

void SemanticVersion::processSuffix()
{
    // Initialize substrings.
    prerelease_ = "";
    build_ = "";

    if (suffix_.empty())
    {
        kind_ = "release";
        return;
    }

    // Split suffix and use only first part if there are spaces.
    std::size_t i = suffix_.find(' ');
    if (i != std::string::npos) { suffix_ = suffix_.substr(0, i); }

    if (suffix_.at(0) == '-')
    {
        i = suffix_.find('+');
        if (i != std::string::npos)
        {
            if (i > 1) prerelease_ = suffix_.substr(1, i - 1);
            if (i < suffix_.size() - 1) build_ = suffix_.substr(i + 1);
            kind_ = "pre-release-dev";
        }
        else
        {
            prerelease_ = suffix_.substr(1);
            kind_ = "pre-release";
        }
    }
    else if (suffix_.at(0) == '+')
    {
        build_ = suffix_.substr(1);
        kind_ = "dev";
    }
    else
    {
        throw std::runtime_error("ERROR: Invalid semantic version, "
                                 "suffix must start with \"-\" or \"+\"");
    }

    return;
}

void SemanticVersion::checkValidity() const
{
    // Check if all version numbers are greater than zero.
    std::string wrongVersion;
    if (major_ < 0) wrongVersion = "major";
    if (minor_ < 0) wrongVersion = "minor";
    if (patch_ < 0) wrongVersion = "patch";
    if (!wrongVersion.empty())
    {
        throw std::runtime_error("ERROR: Invalid semantic version number \"" + wrongVersion
                                 + "\", contains negative value.");
    }

    if (kind_.find("pre-release") != std::string::npos)
    {
        if (prerelease_.empty())
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", empty pre-release string.");
        }
    }

    if (kind_.find("dev") != std::string::npos)
    {
        if (build_.empty())
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", empty build metadata string.");
        }
        if (build_.find("+") != std::string::npos)
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", build metadata contains another \"+\" character.");
        }
    }

    return;
}
