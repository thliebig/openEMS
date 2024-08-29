#include <boost/predef.h>

#if BOOST_ARCH_X86
#include <xmmintrin.h>
#endif

// Disable denormal (subnormal) floating point numbers. These exceedingly
// small numbers may create a substantial overhead depending on the CPU
// (microcode assists are required os x86).
//
// TODO: Only implemented on x86. Do other CPUs like POWER, ARM have
// denormal overheads? If so, implement them too.

namespace Denormal
{
	inline void Disable();
};

inline void Denormal::Disable()
{
#if BOOST_ARCH_X86
	// read the old MXCSR setting
	unsigned int oldMXCSR = _mm_getcsr();

	// set DAZ and FZ bits (flush to zero)
	unsigned int newMXCSR = oldMXCSR | 0x8040;

	// write the new MXCSR setting to the MXCSR
	_mm_setcsr( newMXCSR );
#endif
}
