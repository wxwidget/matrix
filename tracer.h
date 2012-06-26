/** brief this contains basic tracer to trace the memory allocating and deallocating
 *
 */
#ifndef ALICE_MATRIX_TRACER_H
#define ALICE_MATRIX_TRACER_H
#include <stdarg.h>
#include <iostream>

namespace alice
{
namespace matrix
{
class NoneTrace
{
public:
    //call do nothing
    static inline void Trace(std::ostream& stream, const char* format, ...)
    {
    }
    //call vprint
    static inline void Trace(const std::string& str)
    {
    }
};

class Tracer
{
    static inline void Trace(std::ostream& stream, const char* format, ...)
    {
        char buf[256];
        va_list args;
        va_start (args, format);
        vsprintf (buf, format, args);
        va_end (args);
        stream << buf;
    }
    static inline void Trace(const std::string& str)
    {
        std::cout << str;
    }
};
}
}
#endif
