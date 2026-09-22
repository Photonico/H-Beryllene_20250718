#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <map>
#include <string>
#include <tuple>

#if PROFILING
#define VASPML_PROFILING_START(keyword) do { \
      global_scope::timer.start( keyword ); \
   } while(0)
#else
#define VASPML_PROFILING_START(X)
#endif

#if PROFILING
#define VASPML_PROFILING_STOP(keyword) do { \
      global_scope::timer.stop( keyword ); \
   } while(0)
#else
#define VASPML_PROFILING_STOP(X)
#endif

#if PROFILING
#define VASPML_PROFILING_WRITE() do { \
      global_scope::timer.writeToScreen(); \
   } while(0)
#else
#define VASPML_PROFILING_WRITE()
#endif

namespace vaspml
{

/*******************************************************************************************
 * @class Timer
 * @brief implements a class which can be used for code timings
 *
 * This class can be used to time code segements. The timings will be averaged over all
 * calls of the function calls and will also write the absolute time.
 *
 * Usage example:\n
 * @code
 * Timer time; \n
 * time.start( "MyTimerID" ); \n
 * some code segment which will be timed \n
 * time.stop( "MyTimerID" ); \n
 * time.writeToScreen(); \n
 * @endcode
 *******************************************************************************************/
class Timer
{
  public:
    /*******************************************************************************************
     * Start timing with given ID.
     *
     * @param id Representative name for the timer, e.g., function name that is timed.
     *******************************************************************************************/
    void start(const std::string& id);
    /*******************************************************************************************
     * Stop timing with given ID.
     *
     * @param id Representative name for the timer, e.g., function name that is timed.
     *******************************************************************************************/
    void stop(const std::string& id);
    /*******************************************************************************************
     * writing the timing data to screen.
     *
     * The write to screen function will print the total time for evert set timer, the average
     * time for every set timer and the factor between the current timer and the fastest timer
     *******************************************************************************************/
    void writeToScreen();

  private:
    /*******************************************************************************************
     * Storage for start and end times, total time and call count.
     *
     * Tuple content:
     * - 0: start time
     * - 1: stop time
     * - 2: accumulated time in seconds
     * - 3: how often this timer was started
     * - 4: whether timer is currently running
     *******************************************************************************************/

  public:
    typedef std::tuple<std::chrono::time_point<std::chrono::high_resolution_clock>,
                       std::chrono::time_point<std::chrono::high_resolution_clock>,
                       double,
                       size_t,
                       bool>
        TimerStorage;
    /*******************************************************************************************
     * storing timings. The keys of the map are determined by the strings to start/stop timer
     *******************************************************************************************/
    std::map<std::string, TimerStorage> data;
};

namespace global_scope
{
inline Timer timer;
}

} //namespace vaspml

#endif
