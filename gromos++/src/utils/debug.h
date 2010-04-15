/**
 * @file debug.h
 * debugging for GROMOS++
 */

#undef DEBUG

extern int debug_level;

#ifdef NDEBUG
#define DEBUG(level, s)
#else
#define DEBUG(level, s) \
  if ((level) <= ::debug_level){ \
    std::cout << "DEBUG: " << s << std::endl; \
  }

#endif	/* DEBUG_H */

