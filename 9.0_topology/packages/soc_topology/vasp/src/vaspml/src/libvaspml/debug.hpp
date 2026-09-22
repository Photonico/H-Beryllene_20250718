#ifndef DEBUG_HPP
#define DEBUG_HPP

#if VASPML_DEBUG_LEVEL > 0
#define VASPML_DEBUG_L1(X) do { X } while(0)
#else
#define VASPML_DEBUG_L1(X)
#endif

#if VASPML_DEBUG_LEVEL > 1
#define VASPML_DEBUG_L2(X) do { X } while(0)
#else
#define VASPML_DEBUG_L2(X)
#endif

#if VASPML_DEBUG_LEVEL > 2
#define VASPML_DEBUG_L3(X) do { X } while(0)
#else
#define VASPML_DEBUG_L3(X)
#endif

#endif
