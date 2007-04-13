#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <ostream>

#ifdef DEBUG
extern std::ostream& debug;
#else

class DummyOstream {
};

template<class C>
inline DummyOstream& operator<<(DummyOstream& os,
                         const C&)
{
  return os;
}
extern DummyOstream& debug;
#endif

#endif
