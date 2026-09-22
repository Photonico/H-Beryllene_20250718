#include "Timer.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace vaspml;
using hrc = std::chrono::high_resolution_clock;

void Timer::start(const std::string& id)
{
    // Take time immediately, check if timer exists later.
    auto current = hrc::now();

    if (data.find(id) == data.end())
    {
        // New timer, take down current time.
        std::get<2>(data[id]) = (double)0;
        std::get<3>(data[id]) = (size_t)1;
        std::get<4>(data[id]) = true;
        std::get<0>(data[id]) = current;
    }
    else
    {
        // If timer was started again without being stopped first, throw an error.
        if (std::get<4>(data[id]))
        {
            throw std::runtime_error("ERROR: Timer \"" + id + "\" was started already before.");
        };
        // Restarted timer, note current time and increase counter.
        std::get<3>(data[id])++;
        std::get<4>(data[id]) = true;
        std::get<0>(data[id]) = current;
    }

    return;
}

void Timer::stop(const std::string& id)
{
    // Take time immediately, check if timer exists later.
    auto current = hrc::now();

    if (data.find(id) != data.end())
    {
        // Store current stop time.
        std::get<1>(data[id]) = current;
        // If timer was not started before it cannot be stopped.
        if (!std::get<4>(data[id]))
        {
            throw std::runtime_error("ERROR: Timer \"" + id + "\" was not started before.");
        };
        // Accumulate total time of this timer.
        std::get<2>(data[id]) +=
            std::chrono::duration<double>(std::get<1>(data[id]) - std::get<0>(data[id])).count();
        std::get<4>(data[id]) = false;
    }
    else
    {
        throw std::runtime_error("ERROR: New Timer \"" + id + "\" was stopped but never started.");
    }

    return;
}

void Timer::writeToScreen()
{
    // Get maximum total time of all timers.
    using compareType = const std::pair<std::string, TimerStorage>&;
    double maxAverageTimer =
        std::get<2>(std::max_element(data.begin(),
                                     data.end(),
                                     [](compareType a, compareType b) -> bool
                                     { return std::get<2>(a.second) < std::get<2>(b.second); })
                        ->second);

    for (const auto& idata : data)
    {
        std::cout << "TIMER \"" << idata.first << "\":"
                  << " total " << std::fixed << std::setprecision(15) << std::get<2>(idata.second)
                  << " calls " << std::setw(6) << std::get<3>(idata.second) << " average "
                  << std::fixed << std::setprecision(15)
                  << std::get<2>(idata.second) / (double)std::get<3>(idata.second) << " relative "
                  << std::fixed << std::setprecision(5)
                  << std::get<2>(idata.second) / maxAverageTimer << std::endl;
    }

    return;
}

namespace global_scope
{
Timer timer;
}
