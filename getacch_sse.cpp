#include "specifics.h"
#include <cstdio> // use C++ headers, not C (when they exist).
#include <xmmintrin.h>
#include <math.h>
#include "structures.h"
#include "getacch_sse.h"

// don't export statics, foo.
static void getacch_ir3_sse(const size_t nbod, const __m128 x[NPLMAX], const __m128 y[NPLMAX], const __m128 z[NPLMAX], __m128 ir3[NPLMAX]);

static void getacch_ah0_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX],
	const __m128 ir3h[NPLMAX], __m128 &axh0, __m128 &ayh0, __m128 &azh0);

static void getacch_ah1_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX],
	const __m128 ir3h[NPLMAX], const __m128 ir3j[NPLMAX], __m128 axh1[NPLMAX], __m128 ayh1[NPLMAX], __m128 azh1[NPLMAX]);

static void getacch_ah2_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 ir3j[NPLMAX],
	__m128 axh2[NPLMAX], __m128 ayh2[NPLMAX], __m128 azh2[NPLMAX]);

static void getacch_ah3_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh3[NPLMAX], __m128 ayh3[NPLMAX], __m128 azh3[NPLMAX]);


// for checking codegen, force emission of a function as a standalone, no matter what.
// only for testing.
// #define SCREWYOU __attribute((used, noinline, nonnull))
#define SCREWYOU

/// ad-hoc scratch buffer, avoid embiggening the stack too much.
struct scratchpad_t {
	// NPLMAX_V is NPLMAX for __m128, right?
	// correct, NPLMAX_V is NPLMAX / 4 ... but not used consistently (if at all) for function definitions
	__m128 MM_ALIGN16 ir3j[NPLMAX_V];
	__m128 MM_ALIGN16 ir3h[NPLMAX_V];
	__m128 MM_ALIGN16 etaj[NPLMAX_V];

	typedef vec3_t<__m128, NPLMAX> ah_t;
	ah_t ah1;
	ah_t ah2;
	ah_t ah3;
};

// typedef storage_planets_sse_t<__m128, NPLMAX> planets_sse_t;
// typedef vec3_t<__m128, NPLMAX> vec3_sse_t;

namespace {
	//FIXME: not inlined, perhaps we're hitting a limit of some sort (also there's 2 calls, that's a hint).
	//       investigate (for now, kludge).
	SCREWYOU
	FINLINE
	// void xxx_getacch_ir3_sse(const size_t nbod, const vec3_sse_t &v, __m128 ir3[NPLMAX]) {
	void xxx_getacch_ir3_sse(const size_t nbod, const vec3_sse_t &v, __m128 * __restrict const ir3) {
		// in case nbod's tests don't propagate up there.
		if (unlikely(nbod < 1))
			return;

		// I think this is a pretty self explanatory function.
		// Its called twice, once for h = heliocentric coordinates, once for j = jacobian coordinates
		// in both cases the suns position is (0, 0, 0) so its ir3 ends up as NaN (1 / 0).
		// This makes me think its best to pull the sun out of this vector and leave it as just planets.
		// I'll point out other places that make me think this too.

		const __m128 ones = _mm_set1_ps(1.0f);

		for(size_t i = 0; i < nbod; ++i) {
			__m128 x2 = _mm_mul_ps(v.x[i], v.x[i]);
			__m128 y2 = _mm_mul_ps(v.y[i], v.y[i]);
			__m128 z2 = _mm_mul_ps(v.z[i], v.z[i]);

			__m128 r2 = _mm_add_ps(x2, _mm_add_ps(z2, y2));

			__m128 r = _mm_sqrt_ps(r2);
			__m128 r3 = _mm_mul_ps(r2, r);
			ir3[i] = _mm_div_ps(ones, r3);
		}
	}

	SCREWYOU
	void xxx_getacch_ah0_sse(const size_t nbod,
	                         const __m128 mass[NPLMAX],
	                         const __m128 ir3h[NPLMAX],
	                         const vec3_sse_t &h,
	                         __m128 &axh0, __m128 &ayh0, __m128 &azh0 )
	{
		// in case nbod's tests don't propagate up there.
		if (unlikely(nbod < 1))
			return;

		// okay, so the whole idea of all of these functions (ah0, ah1, ah2, ah3) is to calculate different portions of the planets acceleration due to other planets
		// So if we have 4 planets and the sun, then p1, p2, p3, p4 all acceleration each other due to newtons law of gravity F = ma = (Gm1m2/r^2).
		// In this case we are specifically dealing with a system that has a dominant mass -- the sun.  This knowledge (that the sun is the dominant mass) can be exploited
		// to create much larger time steps.  Essentially you calculate the planets (p1, p2, etc) perfect orbit around the sun as if there were no other planets.  You then
		// calculate the permutations that each other planet puts on the planets orbit around the sun (so p1 orbits the sun, but is slightly modified by p2, p3, p4, etc).
		//
		// These functions (ah0, ah1, ah2, ah3) are all calculating different aspects of the planets accelerating each other.  At the end the accelerations are summed to get
		// a final acceleration.  so the acceleration of p1 = ah0 + ah1(p1) + ah2(p1) + ah3(p1).
		//
		// ah0 is special in that its the same for all planets.  So technically we could save axh0, ayh0, azh0 as a float.  The reason I return them as vectors is to make the
		// summing step easier (ah0 + ah1 + ah2, etc).  This may be dumb.
		// So at the end axh0_vector = { axh0_scalar, axh0_scalar, axh0_scalar, axh0_scalar } (a vector with 4 identical values).
		//
		// I can get an explanation from my friend of what each function accounts for.  The fortran code simply refers to them as 'terms'.

		const __m128 zeroes = _mm_set1_ps(0.0f);

		// these exist just so we can bitmask out the sun and p1.  If you look at the scalar version there is a parameter called istart that is 2.  This means we want the sun,
		// and p1 to not have any effect on the value of ah0.
		// again, this is another situation where we have to treat the sun as a special case which makes me think moving it out of the nbody vector could be smart.
		const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
		const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

		axh0 = zeroes;
		ayh0 = zeroes;
		azh0 = zeroes;

		// essentially we are pulling the first iteration out of the loop because we treat the sun and p1 special
		// since we want the sun and p1 to have no affect of the ah0 values if we set their ir3h to zero then it carries through the rest of the equal to finish as zero.
		// this is kind of ghetto (okay, real ghetto) for the sun since at this point its ir3h is NaN.  Again, the sun as a special case ...
		__m128 ir3h0 = _mm_and_ps(ir3h[0], remove_sun_and_p1);

		// okay, so if we look at scalar we get: float fac = mass[i] * ir3h[i];  which means fac for sun / p1 = 0 since ir3h is zero
		__m128 fac = _mm_mul_ps(mass[0], ir3h0);

		// xfac, yfac, zfac are fac * x, y, z ... but we also accumulate in them.  So the first iteration (i = 0) initializes their values, but in the loop we add results to them.
		__m128 xfac = _mm_mul_ps(fac, h.x[0]);
		__m128 yfac = _mm_mul_ps(fac, h.y[0]);
		__m128 zfac = _mm_mul_ps(fac, h.z[0]);

		// we just repeat for 1 to N and keep accumulationg xfac, yfac, zfac
		for(size_t i = 1; i < nbod; ++i)
		{
			fac = _mm_mul_ps(mass[i], ir3h[i]);
			xfac = _mm_add_ps(xfac, _mm_mul_ps(fac, h.x[i]));
			yfac = _mm_add_ps(yfac, _mm_mul_ps(fac, h.y[i]));
			zfac = _mm_add_ps(zfac, _mm_mul_ps(fac, h.z[i]));
		}

		// all this bullshit below is to just do a horizontal add on xfac, yfac, and zfac.
		// actually, I lied, its to do the horizontal add, get a single scalar, and then populate the axh0, ayh0, azh0 with that scalar in all positions (hence the _mm_set1_ps).

		// this version was vector -> scalar -> sum -> vector route

		//float MM_ALIGN16 mx[4], mxtot;
		//_mm_store_ps(mx, xfac);
		//mxtot = mx[0] + mx[1] + mx[2] + mx[3];
		//axh0 = _mm_set1_ps(0.0f - mxtot);

		//float MM_ALIGN16 my[4], mytot;
		//_mm_store_ps(my, yfac);
		//mytot = my[0] + my[1] + my[2] + my[3];
		//ayh0 = _mm_set1_ps(0.0f - mytot);

		//float MM_ALIGN16 mz[4], mztot;
		//_mm_store_ps(mz, zfac);
		//mztot = mz[0] + mz[1] + mz[2] + mz[3];
		//azh0 = _mm_set1_ps(0.0f - mztot);

		// What the fuck are you doing?

		// this version was doing a vector add that I found in the intel optimization guide of how to horizontally sum 4 vectors so that the sum is a single componenent of a final vector.
		// in this case the 4th vector is just zeroes, and vector 1-3 are x, y, and z.
		__m128 t1a = _mm_movelh_ps(yfac, xfac);
		__m128 t1b = _mm_movehl_ps(xfac, yfac);
		__m128 t1c = _mm_movelh_ps(zeroes, zfac);
		__m128 t1d = _mm_movehl_ps(zfac, zeroes);

		__m128 t2a = _mm_add_ps(t1a, t1b);
		__m128 t2b = _mm_add_ps(t1c, t1d);

		__m128 t3a = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(0, 2, 0, 2));
		__m128 t3b = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(1, 3, 1, 3));

		__m128 tot = _mm_add_ps(t3a, t3b); // tot = {sum(xfac), sum(yfac), sum(zfac), sum(zeroes)}

		// populate axh0 with sum(xfac) at all positions ... same for y and z
		axh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(0, 0, 0, 0));
		// make it negative.
		axh0 = _mm_sub_ps(zeroes, axh0);
		ayh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(1, 1, 1, 1));
		ayh0 = _mm_sub_ps(zeroes, ayh0);
		azh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(2, 2, 2, 2));
		azh0 = _mm_sub_ps(zeroes, azh0);
	}

	SCREWYOU
	void xxx_getacch_ah1_sse(const size_t nbod,
	                     const __m128 mass[NPLMAX],
	                     const __m128 ir3h[NPLMAX], const __m128 ir3j[NPLMAX],
	                     const vec3_sse_t &h,
	                     const vec3_sse_t &j,
	                     vec3_sse_t &ah1)
	{
		// in case nbod's tests don't propagate up there.
		if (unlikely(nbod < 1))
			return;

		// in the scalar version we use the mass of the sun for all of our calculations
		// so here we pull that out and a create a mass_sun vector ... again, another hint to make the sun its own vector instead of with the planets.
		float *mass_f = (float*) mass;
		__m128 mass_sun = _mm_set1_ps(mass_f[0]);

		// for reference h = heliocentric (sun = 0,0,0) and j = jacobian (some weird center of mass / vector type coordinate system).
		// so it looks like we are using the difference in coordinate system positions to get some value.  not an astronomer though :)
		//
		// one thing to note is that the scalar version has sun and p1 set to zero.  I don't do that.  I wonder if thats an error or it works out
		// because diff_x, diff_y and diff_z end up being zero?  I'll have to double check that.
		// if I need to add that I would use the bitmask like in ah0 (and, again, would question why the sun is in this vector).
		for(size_t i = 0; i < nbod; ++i) {
			// again, simple fix.
			// load.
			__m128 ah1j, ah1h;
			ah1j = _mm_mul_ps(j.x[i], ir3j[i]);
			ah1h = _mm_mul_ps(h.x[i], ir3h[i]);
			__m128 diff_x = _mm_sub_ps(ah1j, ah1h);
			ah1j = _mm_mul_ps(j.y[i], ir3j[i]);
			ah1h = _mm_mul_ps(h.y[i], ir3h[i]);
			__m128 diff_y = _mm_sub_ps(ah1j, ah1h);
			ah1j = _mm_mul_ps(j.z[i], ir3j[i]);
			ah1h = _mm_mul_ps(h.z[i], ir3h[i]);
			__m128 diff_z = _mm_sub_ps(ah1j, ah1h);

			// store.
			ah1.x[i] = _mm_mul_ps(mass_sun, diff_x);
			ah1.y[i] = _mm_mul_ps(mass_sun, diff_y);
			ah1.z[i] = _mm_mul_ps(mass_sun, diff_z);
		}
	}


	SCREWYOU
	void xxx_getacch_ah2_sse(const size_t nbod,
	                         const __m128 mass[NPLMAX],
	                         const __m128 ir3j[NPLMAX],
	                         const vec3_sse_t &j,
	                         vec3_sse_t &ah2,
	                         const __m128 etaj[NPLMAX])
	{
		// in case nbod's tests don't propagate up there.
		// note: it will be tested at least once for the first loop here.
		if (unlikely(nbod < 1))
			return;

		float * __restrict const etaj_s = (float*) etaj;
		const float *mass_s = (const float*) mass;

		etaj_s[0] = 0.0f;
		etaj_s[1] = mass_s[0];
		//FIXME: that's ugly.
		// okay, this is something I just gave up on.  I didn't know how to vectorize this.
		// etaj should be:
		// index		value
		// 0			0
		// 1			0
		// 2			mass of p1
		// 3			mass of p1 + mass of p2 OR mass of p2 + etaj[2]
		// 4			mass of p1 + mass of p2 + mass of p3 OR mass of p3 + etaj[3]
		// etc ...
		// so I gave up on this and just did it scalar.
		//
		// Technically this only has to be done once since mass is constant.
		for(size_t i = 2; i < nbod * 4; ++i)
			etaj_s[i] = etaj_s[i-1] + mass_s[i-1];

		const __m128 mass_sun = _mm_set1_ps(mass_s[0]);

		const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
		const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

		for(size_t i = 0; i < nbod; ++i) {
			// load
			__m128 fac = _mm_div_ps(_mm_mul_ps(_mm_mul_ps(mass[i], mass_sun), ir3j[i]), etaj[i]);
			__m128 x = _mm_mul_ps(fac, j.x[i]), y = _mm_mul_ps(fac, j.y[i]), z = _mm_mul_ps(fac, j.z[i]);
			// store
			ah2.x[i] = x;
			ah2.y[i] = y;
			ah2.z[i] = z;
		}

		// set the sun and p1 values to zero
		ah2.x[0] = _mm_and_ps(ah2.x[0], remove_sun_and_p1);
		ah2.y[0] = _mm_and_ps(ah2.y[0], remove_sun_and_p1);
		ah2.z[0] = _mm_and_ps(ah2.z[0], remove_sun_and_p1);

		float * __restrict axh2_s = (float*) ah2.x;
		float * __restrict ayh2_s = (float*) ah2.y;
		float * __restrict azh2_s = (float*) ah2.z;

		//FIXME: same problem.
		// again, this is me giving up
		// we want (for axh2):
		// index		value
		// 0			0
		// 1			0
		// 2			fac[2] * x[2]
		// 3			fac[3] * x[3] + fac[2] * x[2] OR fac[3] * x[3] + axh2[2]
		// 4			fac[4] * x[4] + fac[3] * x[3] + fac[2] * x[2] OR fac[4] * x[4] + axh2[3]
		// so again, basically accumulating at each step.
		// what I ended up doing is calculating every axh2 individually (aka axh2[N] = fac[N] * x[N]
		// I then did a scalar for loop which did the summing from 2 to N
		for(size_t i = 2; i < nbod * 4; ++i)
		{
			axh2_s[i] += axh2_s[i-1];
			ayh2_s[i] += ayh2_s[i-1];
			azh2_s[i] += azh2_s[i-1];
		}
	}


}

///FIXME: should better distinguish what's written/updated and what's not. inspect the flow ASAP.
//        perhaps simply discriminate sending vec3_t, should be enough (and already simpler).
// __attribute((used, noinline, nonnull))
void xgetacch_sse(const size_t nbod, planets_sse_t & __restrict p) {
	// __asm("int3");
	// in case nbod's tests don't propagate up there.
	if (unlikely(nbod < 1))
		return;

	scratchpad_t &scratch = *allocate<scratchpad_t>();

	const __m128 zeroes = _mm_set1_ps(0.0f);
	__m128 axh0 = zeroes, ayh0 = zeroes, azh0 = zeroes;

	// get ir3 in both coordinates -- heliocentric and jacobian
	xxx_getacch_ir3_sse(nbod, p.J, scratch.ir3j);
	xxx_getacch_ir3_sse(nbod, p.H, scratch.ir3h);

	__m128 * const __restrict mass = p.MASS.m;

	xxx_getacch_ah0_sse(nbod, mass, scratch.ir3h, p.H, axh0, ayh0, azh0);

	xxx_getacch_ah1_sse(nbod,
	                    mass, scratch.ir3h, scratch.ir3j,
	                    p.H, p.J,
	                    scratch.ah1);

	xxx_getacch_ah2_sse(nbod, mass, scratch.ir3j, p.J, scratch.ah2, scratch.etaj);

	// not touching this one with a 10 foot pole.
	// yah, I'll explain this one in detail ... you will see my madness :)
	getacch_ah3_sse(nbod, mass, p.H.x, p.H.y, p.H.z, scratch.ah3.x, scratch.ah3.y, scratch.ah3.z);

	for(size_t i = 0; i < nbod; ++i) {
		// can you explain the below?  This is good or bad?
		/*
		  obvious fix for re-ordering. perfect.
		  well, ok, could be unrolled, interleaved, but that's for later if needed.
		  opcodes are a tad too large too ;)

		  402690:       0f 28 90 00 05 00 00    movaps 0x500(%rax),%xmm2
		  402697:       0f 28 88 00 06 00 00    movaps 0x600(%rax),%xmm1
		  40269e:       0f 58 90 00 08 00 00    addps  0x800(%rax),%xmm2
		  4026a5:       0f 28 80 00 07 00 00    movaps 0x700(%rax),%xmm0
		  4026ac:       0f 58 88 00 09 00 00    addps  0x900(%rax),%xmm1
		  4026b3:       0f 58 80 00 0a 00 00    addps  0xa00(%rax),%xmm0
		  4026ba:       0f 58 90 00 02 00 00    addps  0x200(%rax),%xmm2
		  4026c1:       0f 58 88 00 03 00 00    addps  0x300(%rax),%xmm1
		  4026c8:       0f 58 80 00 04 00 00    addps  0x400(%rax),%xmm0
		  4026cf:       48 83 c0 10             add    $0x10,%rax
		  4026d3:       0f 58 d5                addps  %xmm5,%xmm2
		  4026d6:       0f 58 cc                addps  %xmm4,%xmm1
		  4026d9:       0f 58 c3                addps  %xmm3,%xmm0
		  4026dc:       41 0f 29 54 15 00       movaps %xmm2,0x0(%r13,%rdx,1)
		  4026e2:       0f 29 4c 15 00          movaps %xmm1,0x0(%rbp,%rdx,1)
		  4026e7:       0f 29 04 13             movaps %xmm0,(%rbx,%rdx,1)
		  4026eb:       48 83 c2 10             add    $0x10,%rdx
		  4026ef:       4c 39 f2                cmp    %r14,%rdx
		  4026f2:       75 9c                   jne    402690 <_Z12xgetacch_ssemPDv4_fPKS_S2_S2_S2_S2_S2_S0_S0_S0_+0x240>

		*/
		// load.
		__m128 x = _mm_add_ps(axh0, _mm_add_ps(scratch.ah1.x[i], _mm_add_ps(scratch.ah2.x[i], scratch.ah3.x[i])));
		__m128 y = _mm_add_ps(ayh0, _mm_add_ps(scratch.ah1.y[i], _mm_add_ps(scratch.ah2.y[i], scratch.ah3.y[i])));
		__m128 z = _mm_add_ps(azh0, _mm_add_ps(scratch.ah1.z[i], _mm_add_ps(scratch.ah2.z[i], scratch.ah3.z[i])));
		// store.
		p.H.x[i] = x;
		p.H.y[i] = y;
		p.H.z[i] = z;
	}

	deallocate(&scratch);
}

__attribute((used, noinline))
void getacch_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh[NPLMAX], __m128 ayh[NPLMAX], __m128 azh[NPLMAX])
{
	const __m128 zeroes = _mm_set1_ps(0.0f);

	__m128 MM_ALIGN16 ir3j[NPLMAX_V];
	__m128 MM_ALIGN16 ir3h[NPLMAX_V];

	getacch_ir3_sse(nbod, xj, yj, zj, ir3j);
	getacch_ir3_sse(nbod, xh, yh, zh, ir3h);

	__m128 axh0 = zeroes, ayh0 = zeroes, azh0 = zeroes;

	getacch_ah0_sse(nbod, mass, xh, yh, zh, ir3h, axh0, ayh0, azh0);

	__m128 MM_ALIGN16 axh1[NPLMAX_V];
	__m128 MM_ALIGN16 ayh1[NPLMAX_V];
	__m128 MM_ALIGN16 azh1[NPLMAX_V];

	getacch_ah1_sse(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h, ir3j, axh1, ayh1, azh1);

	__m128 MM_ALIGN16 axh2[NPLMAX_V];
	__m128 MM_ALIGN16 ayh2[NPLMAX_V];
	__m128 MM_ALIGN16 azh2[NPLMAX_V];

	getacch_ah2_sse(nbod, mass, xj, yj, zj, ir3j, axh2, ayh2, azh2);

	__m128 MM_ALIGN16 axh3[NPLMAX_V];
	__m128 MM_ALIGN16 ayh3[NPLMAX_V];
	__m128 MM_ALIGN16 azh3[NPLMAX_V];

	getacch_ah3_sse(nbod, mass, xh, yh, zh, axh3, ayh3, azh3);

	for(size_t i = 0; i < nbod; ++i)
	{
		axh[i] = _mm_add_ps(axh0, _mm_add_ps(axh1[i], _mm_add_ps(axh2[i], axh3[i])));
		ayh[i] = _mm_add_ps(ayh0, _mm_add_ps(ayh1[i], _mm_add_ps(ayh2[i], ayh3[i])));
		azh[i] = _mm_add_ps(azh0, _mm_add_ps(azh1[i], _mm_add_ps(azh2[i], azh3[i])));
	}

	//for(size_t i = 1; i < nbod; ++i)
	//{
	//	axh[i] = axh0 + axh1[i] + axh2[i] + axh3[i];
	//	ayh[i] = ayh0 + ayh1[i] + ayh2[i] + ayh3[i];
	//	azh[i] = azh0 + azh1[i] + azh2[i] + azh3[i];
	//}

}

static void getacch_ir3_sse(const size_t nbod, const __m128 x[NPLMAX], const __m128 y[NPLMAX], const __m128 z[NPLMAX], __m128 ir3[NPLMAX])
{
	const __m128 ones = _mm_set1_ps(1.0f);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 x2 = _mm_mul_ps(x[i], x[i]);
		__m128 y2 = _mm_mul_ps(y[i], y[i]);
		__m128 z2 = _mm_mul_ps(z[i], z[i]);

		__m128 r2 = _mm_add_ps(x2, _mm_add_ps(z2, y2));

		__m128 r = _mm_sqrt_ps(r2);
		__m128 r3 = _mm_mul_ps(r2, r);
		ir3[i] = _mm_div_ps(ones, r3);
	}
}

static void getacch_ah0_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], 
	const __m128 ir3h[NPLMAX], __m128 &axh0, __m128 &ayh0, __m128 &azh0)
{
	const __m128 zeroes = _mm_set1_ps(0.0f);
	const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

	axh0 = zeroes;
	ayh0 = zeroes;
	azh0 = zeroes;

	__m128 ir3h0 = _mm_and_ps(ir3h[0], remove_sun_and_p1);
	__m128 fac = _mm_mul_ps(mass[0], ir3h0);

	__m128 xfac = _mm_mul_ps(fac, xh[0]);
	__m128 yfac = _mm_mul_ps(fac, yh[0]);
	__m128 zfac = _mm_mul_ps(fac, zh[0]);

	for(size_t i = 1; i < nbod; ++i)
	{
		fac = _mm_mul_ps(mass[i], ir3h[i]);
		xfac = _mm_add_ps(xfac, _mm_mul_ps(fac, xh[i]));
		yfac = _mm_add_ps(yfac, _mm_mul_ps(fac, yh[i]));
		zfac = _mm_add_ps(zfac, _mm_mul_ps(fac, zh[i]));
	}

	//float MM_ALIGN16 mx[4], mxtot;
	//_mm_store_ps(mx, xfac);
	//mxtot = mx[0] + mx[1] + mx[2] + mx[3];
	//axh0 = _mm_set1_ps(0.0f - mxtot);

	//float MM_ALIGN16 my[4], mytot;
	//_mm_store_ps(my, yfac);
	//mytot = my[0] + my[1] + my[2] + my[3];
	//ayh0 = _mm_set1_ps(0.0f - mytot);

	//float MM_ALIGN16 mz[4], mztot;
	//_mm_store_ps(mz, zfac);
	//mztot = mz[0] + mz[1] + mz[2] + mz[3];
	//azh0 = _mm_set1_ps(0.0f - mztot);

	__m128 t1a = _mm_movelh_ps(yfac, xfac);
	__m128 t1b = _mm_movehl_ps(xfac, yfac);
	__m128 t1c = _mm_movelh_ps(zeroes, zfac);
	__m128 t1d = _mm_movehl_ps(zfac, zeroes);

	__m128 t2a = _mm_add_ps(t1a, t1b);
	__m128 t2b = _mm_add_ps(t1c, t1d);

	__m128 t3a = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(0, 2, 0, 2));
	__m128 t3b = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(1, 3, 1, 3));

	__m128 tot = _mm_add_ps(t3a, t3b);

	axh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(0, 0, 0, 0));
	axh0 = _mm_sub_ps(zeroes, axh0);
	ayh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(1, 1, 1, 1));
	ayh0 = _mm_sub_ps(zeroes, ayh0);
	azh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(2, 2, 2, 2));
	azh0 = _mm_sub_ps(zeroes, azh0);
}

static void getacch_ah1_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 ir3h[NPLMAX], const __m128 ir3j[NPLMAX], __m128 axh1[NPLMAX], __m128 ayh1[NPLMAX], __m128 azh1[NPLMAX])
{
	float *mass_f = (float*) mass;
	__m128 mass_sun = _mm_set1_ps(mass_f[0]);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 ah1j = _mm_mul_ps(xj[i], ir3j[i]);
		__m128 ah1h = _mm_mul_ps(xh[i], ir3h[i]);
		__m128 diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		axh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);

		ah1j = _mm_mul_ps(yj[i], ir3j[i]);
		ah1h = _mm_mul_ps(yh[i], ir3h[i]);
		diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		ayh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);

		ah1j = _mm_mul_ps(zj[i], ir3j[i]);
		ah1h = _mm_mul_ps(zh[i], ir3h[i]);
		diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		azh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);
	}
}

static void getacch_ah2_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 ir3j[NPLMAX], 
	__m128 axh2[NPLMAX], __m128 ayh2[NPLMAX], __m128 azh2[NPLMAX])
{
	float *etaj_s;
	__m128 *etaj;
	__m128 mass_sun;

	etaj_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	const float *mass_s = (float*)mass;

	etaj = (__m128*) etaj_s;

	etaj_s[0] = 0.0f;
	etaj_s[1] = mass_s[0];
	for(size_t i = 2; i < nbod * 4; ++i)
	{
		etaj_s[i] = etaj_s[i-1] + mass_s[i-1];
	}

	mass_sun = _mm_set1_ps(mass_s[0]);

	const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 fac = _mm_div_ps(_mm_mul_ps(_mm_mul_ps(mass[i], mass_sun), ir3j[i]), etaj[i]);
		axh2[i] = _mm_mul_ps(fac, xj[i]);
		ayh2[i] = _mm_mul_ps(fac, yj[i]);
		azh2[i] = _mm_mul_ps(fac, zj[i]);
	}

	axh2[0] = _mm_and_ps(axh2[0], remove_sun_and_p1);
	ayh2[0] = _mm_and_ps(ayh2[0], remove_sun_and_p1);
	azh2[0] = _mm_and_ps(azh2[0], remove_sun_and_p1);

	float *axh2_s = (float*)axh2;
	float *ayh2_s = (float*)ayh2;
	float *azh2_s = (float*)azh2;

	for(size_t i = 2; i < nbod * 4; ++i)
	{
		axh2_s[i] += axh2_s[i-1];
		ayh2_s[i] += ayh2_s[i-1];
		azh2_s[i] += azh2_s[i-1];
	}
}

static void getacch_ah3_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh3[NPLMAX], __m128 ayh3[NPLMAX], __m128 azh3[NPLMAX])
{

	// okay, this is actually the most important function in terms of scaling -- its O(n^2)
	// what we actually want to do is calculate the acceleration on pX due to every other planet (p of not X).

	// the end resust is the axh3[N] is the sum of all of the accelerations of every other planet on pN.

	__m128 zeroes = _mm_set1_ps(0.0f);
	__m128 ones = _mm_set1_ps(1.0f);


	// since we are going to accumulate  we first clear out the accelerations.
	for(size_t i = 0; i < nbod; ++i)
	{
		axh3[i] = zeroes;
		ayh3[i] = zeroes;
		azh3[i] = zeroes;
	}

	// okay, again, the sun doesn't take part in this ... god, I should move it to out of the planet vector.  I don't know why I resist.
	// so here we set the mass of the sun to zero which means the acceleration it causes on all of the planets is also zero.
	const unsigned int MM_ALIGN16 mask[4] = { 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun =_mm_load_ps((float*) mask);

	mass[0] = _mm_and_ps(mass[0], remove_sun);

		// so this is all confusing and stupid as fuck.  Let me copy / past the scalar code and comment on it.

		/*

			okay ... so for 8 planets we now have a matrix
			i/j	0	1	2	3	4	5	6	7
			0	0
			1	X	0
			2	X	X	0
			3	X	X	X	0
			4	X	X	X	X	0
			5	X	X	X	X	X	0
			6	X	X	X	X	X	X	0
			7	X	X	X	X	X	X	X	0

			So obviously we want the effect of p0 on p0 to be zero.  Repeat for every other planet.
			
			Now, we also calculate the effect of P1 on P2 AND P2 on P1 at the same time so the loop only does each pairing once.
			So these are marked with X ... the loop skips them since j = i + 1.

			Everything else we calculate.

			// skip the sun ... i = 1
			for(size_t i = 1; i < nbod - 1; ++i)
			{
				for(size_t j = i + 1; j < nbod; ++j)
				{
					// pretty easy, get the diff x, y, z
					float dx = xh[j] - xh[i];
					float dy = yh[j] - yh[i];
					float dz = zh[j] - zh[i];

					// again, easy, get the vector between them (well, squared)
					float rji2 = dx*dx + dy*dy + dz*dz;
					// inverse cube
					float irij3 = 1.0f / (rji2 * sqrtf(rji2));

					// so these values - faci, facj - are the effect of Pj on Pi (facj) and Pi on Pj (faci).
					float faci = mass[i] * irij3;
					float facj = mass[j] * irij3;

					// we now want to accumulate these accelerations (faci/facj) in the respective arrays
					axh3[j] -= faci*dx;
					ayh3[j] -= faci*dy;
					azh3[j] -= faci*dz;

					axh3[i] += facj*dx;
					ayh3[i] += facj*dy;
					azh3[i] += facj*dz;
				}
			}

			So, to me, there are two tricky parts when trying to vectorize this.

			1) vec[N] also acts on itself.  I mean, vec[N] contains 4 planets (p0, p1, p2, p3).  This means we have to do interactions:

				p0 <-> p1
				p0 <-> p2
				p0 <-> p3
				p1 <-> p2
				p1 <-> p3
				p2 <-> p3

			The way I accopmlished this is with rotating the vector by 1 and then working it on itself.  E.g.

				vec				=> p0, p1, p2, p3
				vec_rotate_1	=> p3, p0, p1, p2
				results:
								p0 <-> p3
								p1 <-> p0
								p2 <-> p1
								p3 <-> p2

				vec				=> p0, p1, p2, p3
				vec_rotate_2	=> p2, p3, p0, p1
				results:
								p0 <-> p2
								p1 <-> p3
								p2 <-> p0
								p3 <-> p1

				vec				=> p0, p1, p2, p3
				vec_rotate_2	=> p1, p2, p3, p0
				results:
								p0 <-> p1
								p1 <-> p2
								p2 <-> p3
								p3 <-> p0

			What I didn't realize at the time (and actually didn't fully realize until typing this out) was that doing it
			this way results in calculating both p2 <-> p3 AND p3 <-> p2.  This is why I only accumulate the acceration and never do a subtraction.

			This obviously also results in doing 12 interactions instead of doing 6.  I'm dumb.

			2) vec[N] acts on vec[N+1] we only do some of the interactions instead of doing all of them.  E.g.

				vec[N]			=> p0, p1, p2, p3
				vec[N+1]		=> p4, p5, p6, p7
				results:
								p0 <-> p4
								p1 <-> p5
								p2 <-> p6
								p3 <-> p7

				We still need more.  Again, I solve this with rotations:

				vec[N]			=> p0, p1, p2, p3
				vec[N+1]_rotate	=> p7, p4, p5, p6
				results:
								p0 <-> p7
								p1 <-> p4
								p2 <-> p5
								p3 <-> p6
				
				Again, repeat two more times until I end up getting: p0 <-> p4, p0 <-> p5, p0 <-> p6, p0 <-> p7, and the same for p1, p2, p3.


			***

			I have a feeling this way kind of sucks.  Another option would be to do something like

				vec[N]			=> p0, p1, p2, p3
				vec[N+1]		=> p4, p5, p6, p7

				I could move p0 horizontally:

				vec_p0			=> p0, p0, p0, p0
				vec[N+1]		=> p4, p5, p6, p7

				I then do N+2, N+3, etc until I'm done.  Finally, I horizontally sum the results to get ah3 of p0.

				Repeat for p1, p2, p3.

			This is by far the most important function, and it so far its implementation is terrible.

			I'll now walk through the implementation :)
		*/

	// we are going to go through each planet vector (nbods)!
	for(size_t i = 0; i < nbod; ++i)
	{

		// we make local copies of xh[i] as xhi (and others) because we are going to rotate them using SHUFFLE.
		__m128 xacci = zeroes, yacci = zeroes, zacci = zeroes;
		__m128 xhi = xh[i], yhi = yh[i], zhi = zh[i], massi = mass[i];

		__m128 dx, dy, dz, rji2, irji3, facj, faci;


		// okay, so this operation does the following accelerations:
		//
		// veci			=> p0, p1, p2, p3
		// vecj			=> p4, p5, p6, p7
		//
		// results		=> p0 <-> p4, p1 <-> p5, p2 <-> p6, p3 <-> p7
		//
		// results of accelerations for veci are stored in xacci, yacci, zacci
		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		// we now rotate damn near everyting for i
		// acceleration results (xacci, yacci, zacci)
		// the positions (xhi, yhi, zhi)
		// the mass (massi)

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		__m128 xhi_rot1 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 yhi_rot1 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 zhi_rot1 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 mass_rot1 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(0, 3, 2, 1));

		// now we do the first inter vector i accelerations
		// veci			=> p0, p1, p2, p3
		// vec_rot1		=> p3, p0, p1, p2

		dx = _mm_sub_ps(xhi_rot1, xh[i]);
		dy = _mm_sub_ps(yhi_rot1, yh[i]);
		dz = _mm_sub_ps(zhi_rot1, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot1, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot1;
		yhi = yhi_rot1;
		zhi = zhi_rot1;
		massi = mass_rot1;

		// now we do veci and vecj again, but this time we use the rotated veci
		//
		// veci			=> p3, p0, p1, p2
		// vecj			=> p4, p5, p6, p7
		//
		// results		=> p3 <-> p4, p0 <-> p5, p1 <-> p6, p2 <-> p7

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		// repeat the same rotations ...

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		__m128 xhi_rot2 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 yhi_rot2 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 zhi_rot2 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 mass_rot2 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(1, 0, 3, 2));

		// calculate the inter vector i accelerations

		dx = _mm_sub_ps(xhi_rot2, xh[i]);
		dy = _mm_sub_ps(yhi_rot2, yh[i]);
		dz = _mm_sub_ps(zhi_rot2, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot2, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot2;
		yhi = yhi_rot2;
		zhi = zhi_rot2;
		massi = mass_rot2;

		// doing veci rotated by 2 spots on vecj
		//
		// veci			=> p2, p3, p0, p1
		// vecj			=> p4, p5, p6, p7
		//
		// results		=> p2 <-> p4, p3 <-> p5, p0 <-> p6, p1 <-> p7

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		// rotate again!  yay!

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));


		__m128 xhi_rot3 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 yhi_rot3 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 zhi_rot3 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 mass_rot3 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(2, 1, 0, 3));

		dx = _mm_sub_ps(xhi_rot3, xh[i]);
		dy = _mm_sub_ps(yhi_rot3, yh[i]);
		dz = _mm_sub_ps(zhi_rot3, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot3, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot3;
		yhi = yhi_rot3;
		zhi = zhi_rot3;
		massi = mass_rot3;

		// last rotated veci on vecj
		//
		// veci			=> p1, p2, p3, p0
		// vecj			=> p4, p5, p6, p7
		//
		// results		=> p1 <-> p4, p2 <-> p5, p3 <-> p6, p0 <-> p7
		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		// we shuffle only the accerations to get them back into position
		//
		// what I mean is that the acceleration for p0 has gone:
		// rot 0		=> { p0,  x,  x,  x }
		// rot 1		=> {  x, p0,  x,  x }
		// rot 2		=> {  x,  x, p0,  x }
		// rot 3		=> {  x,  x,  x, p0 }
		// rot 4		=> { p0,  x,  x,  x }
		//
		// so now its back in the proper stop, we add it to the final acceleration vector
		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		axh3[i] = _mm_add_ps(axh3[i], xacci);
		ayh3[i] = _mm_add_ps(ayh3[i], yacci);
		azh3[i] = _mm_add_ps(azh3[i], zacci);

	}
}
