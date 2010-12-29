#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include "specifics.h"
#include <xmmintrin.h> // for __m128, pff.
#include <cstdlib>

/*
const size_t NPLMAX = 64;
const size_t NTPMAX = 8096;
const size_t NPLMAX_V = 16;
*/

enum {
	NPLMAX   = 64,
	NTPMAX   = 8192, // 8*1024 = 8192 foo, not 8096.
	NPLMAX_V = 16,
};


struct SimulationParameters
{
	float t0, tstop, dt, tinc;
	float dtout, dtdump;
	bool lflag[6];
	char lflag_c[6];
	float rmin, rmax, rmaxu, qmin;
	bool lclose;
	char lclose_c;
	float rcrit;
	char output_file_name[256];
	char frame[4];
	char fopenstat[8];
};


struct NBodies
{
	int nbod, npl;
	float j2rp2, j4rp4;

	float *mass;
	float *rpl;

	float *xh;
	float *yh;
	float *zh;

	float *vxh;
	float *vyh;
	float *vzh;

	float *xj;
	float *yj;
	float *zj;

	float *vxj;
	float *vyj;
	float *vzj;

	void SetNBod(int p_nbod)
	{
		nbod = p_nbod;

		mass = new float[nbod];
		rpl = new float[nbod];

		xh = new float[nbod];
		yh = new float[nbod];
		zh = new float[nbod];

		vxh = new float[nbod];
		vyh = new float[nbod];
		vzh = new float[nbod];

		xj = new float[nbod];
		yj = new float[nbod];
		zj = new float[nbod];

		vxj = new float[nbod];
		vyj = new float[nbod];
		vzj = new float[nbod];
	}
};


/// round up to a power of 2.
/// @note need a compile-time constant.
template<size_t n, size_t denom> struct round_up_to_a_power_of_two {
	enum {
		value = denom > 1 ? ((n + denom - 1) & -denom)/denom : n,
	};
};

/// storage for 3 coordinates/scalars in SoA form.
/// @note not over-generalized to avoid having to refactor code too much (templating the arity of that 'vector').
template<typename T, size_t NumScalars>
	struct vec3_t {
		enum { capacity = round_up_to_a_power_of_two<NumScalars, sizeof(T)/sizeof(float)>::value };
		typedef T type;
		T x[capacity];
		T y[capacity];
		T z[capacity];
	};
/// storage for 1 coordinate/scalar, to maintain symmetry with vec3_t
/// @sa vec3_t
template<typename T, size_t NumScalars>
	struct vec1_t {
		enum { capacity = round_up_to_a_power_of_two<NumScalars, sizeof(T)/sizeof(float)>::value };
		typedef T type;
		T m[capacity];
		T operator[](size_t i) const { return m[i]; }
		T &operator[](size_t i) { return m[i]; }
	};

/// still strongly typed storage for planet related stuff.
template<typename T, size_t NumScalars>
	struct storage_planets_t {
		// using capitals to distinguish a bit.
		//FIXME: figure out & reorganize by access.
		vec1_t<T, NumScalars> MASS;
		vec3_t<T, NumScalars> H, AH, J;
		vec3_t<T, NumScalars> VJ, VH;
		vec1_t<T, NumScalars> RPL;
	};

/// still strongly typed storage for particle related stuff.
template<typename T, size_t NumScalars>
	struct storage_particles_t {
		// using capitals to distinguish a bit.
		//FIXME: figure out & reorganize by access.
		vec3_t<T, NumScalars> HT, VHT, AHT;
	};

/// degenerated storage for planets & particles.
/// @attention unions do stink, but we're selecting early (keep fingers crossed and verify).
template<size_t MaxPlanets, size_t MaxParticles>
	union MM_ALIGN64 storage_t {
		struct scalar_t {
			storage_planets_t<float, MaxPlanets> planets;
			storage_particles_t<float, MaxParticles> particles;
		};
		struct sse_t {
			storage_planets_t<__m128, MaxPlanets> planets;
			storage_particles_t<__m128, MaxParticles> particles;
		};

		scalar_t scalar;
		sse_t sse;
	};

// now, some common types.
typedef storage_planets_t<float,  NPLMAX> planets_scalar_t;
typedef storage_planets_t<__m128, NPLMAX> planets_sse_t;
typedef storage_particles_t<float,  NTPMAX> particles_scalar_t;
typedef storage_particles_t<__m128, NTPMAX> particles_sse_t;
typedef vec3_t<float, NPLMAX>  vec3_scalar_t;
typedef vec3_t<__m128, NPLMAX> vec3_sse_t;





// while we're at it, put that crap here.
namespace {
	/// automatically allocate with the proper alignment.
	template<typename T> FA_MALLOC T *allocate() {
		// ah fuck, i always forget alignof() is C++0X (and __alignof__ isn't portable at all).
		return static_cast<T*>(_mm_malloc(sizeof(T), alignof(T)));
	}

	template<typename T> void deallocate(T *p) {
		_mm_free(p);
	}
}

#endif /* STRUCTURES_H_ */
