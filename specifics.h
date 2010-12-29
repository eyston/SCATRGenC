#ifndef SPECIFICS_H
#define SPECIFICS_H

#if defined(__GNUC__)
	#define __GCC__
#elif defined(_MSC_VER)
	#define __MSVC__
#else
	#error unsupported compiler.
#endif

#if defined __GCC__
	#define MM_ALIGN(x)    __attribute__((aligned(x)))

	/// force inlining.
	#define FINLINE        __attribute__((always_inline))
	/// prohibit inlining.
	#define NOINLINE       __attribute__((noinline))
	/// function attribute that says: returned pointer doesn't alias.
	#define FA_MALLOC      __attribute__((malloc))

	// pretend we have C++0X
	#define alignof(x)     __alignof__(x)

	/// branch/condition is expected to hold.
	#define likely(expr)   __builtin_expect((expr), 1)
	/// branch/condition is expected to not hold.
	#define unlikely(expr) __builtin_expect((expr), 0)

#elif defined __MSVC__
	#define MM_ALIGN(x)    __declspec(align(x))
	// pretend we have C++0X
	#define alignof(x)     __alignof(x)

	/// function attribute that says: returned pointer doesn't alias.
	#define FA_MALLOC      __declspec(restrict)
	/// force inlining.
	#define FINLINE        __inline __forceinline
	/// prohibit inlining.
	#define NOINLINE       __declspec(noinline)

	#define likely(expr)   (expr)
	#define unlikely(expr) (expr)
#endif

#define MM_ALIGN16     MM_ALIGN(16)
#define MM_ALIGN64     MM_ALIGN(64)

#endif // SPECIFICS_H
