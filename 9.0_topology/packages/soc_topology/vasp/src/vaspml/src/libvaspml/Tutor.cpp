#include "Tutor.hpp"

#include <iostream>
#include <stdexcept>

using namespace vaspml;

void Tutor::warning(const std::string message) const
{
    std::cout << headerWarning;
    std::cout << message;
    std::cout << footerWarning;
}

void Tutor::error(const std::string message) const
{
    std::cout << headerError;
    std::cout << message;
    std::cout << footerError;
    throw std::runtime_error("");
}

void Tutor::bug(const std::string message) const
{
    std::cout << headerBug;
    std::cout << message;
    std::cout << footerBug;
    throw std::runtime_error("");
}

namespace global_scope
{
Tutor tutor;
}
