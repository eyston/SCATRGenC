#include <cmath>

#include "drift.h"
#include "structures.h"


// Returns the real root of cubic often found in solving kepler
// problem in universal variables
//             Input:
//                 dt            ==>  time step (real scalar)
//                 r0            ==>  Distance between `Sun' and paritcle
//                                     (real scalar)
//                 mu            ==>  Reduced mass of system (real scalar)
//                 alpha         ==>  Twice the binding energy (real scalar)
//                 u             ==>  Vel. dot radial vector (real scalar)
//             Output:
//                 s             ==>  solution of cubic eqn for the  
//                                    universal variable
//                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
static void drift_kepu_p3solve(const float dt, const float r0, const float mu, const float alpha, const float u, float &s, int &iflg)
{
	float denom, a0, a1, a2, q, r, sq2, sq, p1, p2;

	denom = (mu - alpha * r0) / 6.0f;
	a2 = 0.5f * u / denom;
	a1 = r0 / denom;
	a0 = -dt / denom;

	q = (a1 - a2 * a2 / 3.0f) / 3.0f;
	r = (a1 * a2 - 3.0f * a0) / 6.0f - powf(a2, 3) / 27.0f;
	sq2 = powf(q, 3) + powf(r, 3);

	if(sq2 > 0.0f)
	{
		sq = sqrtf(sq2);

		if((r + sq) < 0.0f)
		{
			p1 = -powf(-(r + sq), 1.0f / 3.0f);
		}
		else
		{
			p1 = powf(r + sq, 1.0f/3.0f);
		}

		if((r - sq) < 0.0f)
		{
			p2 = -powf(-(r - sq), 1.0f / 3.0f);
		}
		else
		{
			p2 = powf(r - sq, 1.0f/3.0f);
		}

		iflg = 0;
		s = p1 + p2 - a2 / 3.0f;
	}
	else
	{
		iflg = 1;
		s = 0;
	}
}



/***********************************************************************
c	                  ORBEL_SCGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sin and cos.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             sx    ==>  sin(angle)  (real scalar)
*             cx    ==>  cos(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: The HP 700 series won't return correct answers for sin
*       and cos if the angle is bigger than 3e7. We first reduce it
*       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************/

#define PI 3.14159265358979f
#define TWOPI 2 * PI
#define PIBY2 PI / 2.0f

// E: the comments sound ancient ... we can probably just replace with sinf(angle) and cosf(angle)
static void orbel_scget(const float angle, float &sx, float &cx)
{
	int nper;
	float x, PI3BY2 = 1.5f * PI;

	nper = angle / TWOPI;
	x = angle - nper * TWOPI;
	if(x < 0.0f)
	{
		x = x + TWOPI;
	}

	sx = sinf(x);
	cx = sqrtf(1.0f - sx * sx);
	if(x > PIBY2 && x < PI3BY2)
	{
		cx = -cx;
	}
}

// http://gcc.gnu.org/onlinedocs/gfortran/SIGN.html

/*
Description:
SIGN(A,B) returns the value of A with the sign of B. 
Standard:
Fortran 77 and later 
Class:
Elemental function 
Syntax:
RESULT = SIGN(A, B) 
Arguments:
A	Shall be of type INTEGER or REAL 
B	Shall be of the same type and kind as A 


Return value:
The kind of the return value is that of A and B. If B\ge 0 then the result is ABS(A), else it is -ABS(A). 
*/

// E: some fortran library function ... obviously can be improved upon :p
static float sign(float a, float b)
{
	if(b > 0.0f)
	{
		return fabs(a);
	}
	else
	{
		return -fabs(a);
	}
}

/*c*************************************************************************
c                        DRIFT_KEPU_GUESS.F
c*************************************************************************
c Initial guess for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  initial guess for the value of 
c                                    universal variable
c
c Author:  Hal Levison & Martin Duncan 
c Date:    3/12/93
c Last revision: April 6/93 */

static void drift_kepu_guess(const float dt, const float r0, const float mu, const float alpha, const float u, float &s)
{
	int iflg;
	float y, sy, cy, sigma, es;
	float x, a;
	float en, ec, e;

	if(alpha > 0.0f)
	{
		// find initial guess for elliptic motion
		if( dt / r0 < 0.4f)
		{
			s = dt / r0 - (dt * dt * u) / (2.0f * r0 * r0 * r0);
		}
		else
		{
			a = mu / alpha;
			en = sqrtf(mu / (a * a * a));
			ec = 1.0f - r0 / a;
			es = u / (en * a * a);
			e = sqrtf(ec * ec  + es * es);
			y = en * dt - es;
			orbel_scget(y, sy, cy);
			sigma = sign(1.0f, es * cy + ec * sy);
			x = y + sigma * 0.85 * e;
			s = x / sqrtf(alpha);
		}
	}
	else
	{
		// find initial guess for hyperbolic motion.
		drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg);
		if(iflg != 0)
		{
			s = dt / r0;
		}
	}
}


/*c*************************************************************************
c                        DRIFT_KEPU_STUMPFF.F
c*************************************************************************
c subroutine for the calculation of stumpff functions
c see Danby p.172  equations 6.9.15
c
c             Input:
c                 x             ==>  argument
c             Output:
c                 c0,c1,c2,c3   ==>  c's from p171-172
c                                       (real scalors)
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 2/3/93 */

static void drift_kepu_stumpff(float x, float &c0, float &c1, float &c2, float &c3)
{
	int n;
	float xm;

	n = 0;
	xm = 0.0f;
	while(fabs(x) > xm)
	{
		n = n + 1;
		x = x / 4.0f;
	}

	c2 = (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x / 182.0f) / 132.0f) / 90.0f) / 56.0f) / 30.0f) / 12.0f) / 2.0f;
	c3 = (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x * (1.0f - x / 210.0f) / 156.0f) / 110.0f) / 72.0f) / 42.0f) / 20.0f) / 6.0f;
	c1 = 1.0f - x * c3;
	c0 = 1.0f - x * c2;

	if(n != 0)
	{
		// do i=n,1,-1 ... I think this is the correct for statement
		for(int i = n; i > 0; --i)
		{
			c3 = (c2 + c0 * c3) / 4.0f;
			c2 = (c1 * c1) / 2.0f;
			c1 = c0 * c1;
			c0 = 2.0f * c0 * c0 - 1.0f;
			x = x * 4.0f;
		}
	}
}


/*c*************************************************************************
c                        DRIFT_KEPU_NEW.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using NEWTON'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 4/21/93 */

#define DANBYB 1.0e-13 // E: hmmm ... is this a double or float?  who knows!

static void drift_kepu_new(float s, const float dt, const float r0, const float mu, const float alpha, const float u, float &fp, float &c1, float &c2, float &c3, int &iflgn)
{
	int nc;
	float x, c0, ds, f, fpp, fppp, fdt;

	for(int i = 0; i < 7; ++i)
	{
		x = s * s * alpha;
		drift_kepu_stumpff(x, c0, c1, c2, c3);
		c1 = c1 * s;
		c2 = c2 * s * s;
		c3 = c3 * s * s * s;
		f = r0 * c1 + u * c2 + mu * c3 - dt;
		fp = r0 * c0 + u * c1 + mu * c2;
		fpp = (-r0 * alpha + mu) * c1 + u * c0;
		fppp = (-r0 * alpha +mu) * c0 - u * alpha * c1;
		ds = - f / fp;
		ds = - f / (fp + ds * fpp / 2.0f);
		ds = - f / (fp + ds * fpp / 2.0f + ds * ds * fppp / 6.0f);
		s = s + ds;
		fdt = f / dt;

		// quartic convergence
		if(fdt * fdt < DANBYB * DANBYB)
		{
			iflgn = 0;
		}
		// newton's method succeeded
	}

	// newton's method failed
	iflgn = 1;
}


/*c*************************************************************************
c                        DRIFT_KEPU_FCHK.F
c*************************************************************************
c Returns the value of the function f of which we are trying to find the root
c in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 r0            ==>  Distance between `Sun' and particle
c                                     (real scalar)
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 u             ==>  Vel. dot radial vector (real scalar)
c                 s             ==>  Approx. root of f 
c             Output:
c                 f             ==>  function value ( = 0 if O.K.) (integer)
c
c Author:  Martin Duncan  
c Date:    March 12/93
c Last revision: March 12/93 */

static void drift_kepu_fchk(const float dt, const float r0, const float mu, const float alpha, const float u, const float s, float &f)
{
	float x, c0, c1, c2, c3;

	x = s * s * alpha;
	drift_kepu_stumpff(x, c0, c1, c2, c3);
	c1 = c1 * s;
	c2 = c2 * s * s;
	c3 = c3 * s * s * s;
	f = r0 * c1 + u * c2 + mu * c3 - dt;
}


/*c*************************************************************************
c                        DRIFT_KEPU_LAG.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using LAGUERRE'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 4/21/93 */

#define NLAG2 400

static void drift_kepu_lag(float &s, const float dt, const float r0, const float mu, const float alpha, const float u, float &fp, float &c1, float &c2, float &c3, int &iflg)
{
	int nc, ncmax;
	float ln, x, fpp, ds, c0, f, fdt;

	int NTMP = NLAG2 + 1;

	// To get close approch needed to take lots of iterations if alpha<0
	if(alpha < 0.0f)
	{
		ncmax = NLAG2;
	}
	else
	{
		ncmax = NLAG2;
	}

	ln = 5.0f;
	// start laguere's method

	// do nc =0,ncmax
	for(nc = 0; nc < ncmax; ++nc)
	{
		x = s * s * alpha;
		drift_kepu_stumpff(x, c0, c1, c2, c3);
		c1 = c1 * s;
		c2 = c2 * s * s;
		c3 = c3 * s * s * s;
		f = r0 * c1 + u * c2 + mu * c3 - dt;
		fp = r0 * c0 + u * c1 + mu * c2;
		fpp = (-40.0f * alpha + mu) * c1 + u * c0;
		ds = - ln * f / (fp + sign(1.0f, fp) * sqrtf(fabs((ln - 1.0f) * (ln - 1.0f) * fp * fp - (ln - 1.0f) * ln * f * fpp)));
		s = s + ds;

		fdt = f/dt;

		// quartic convergence
		if(fdt * fdt < DANBYB * DANBYB)
		{
			iflg = 0;
			return;
		}
		// Laguerre's method succeeded
	}

	iflg = 2;
}


/*c*************************************************************************
c                        DRIFT_KEPU.F
c*************************************************************************
c subroutine for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflg          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison  
c Date:    2/3/93
c Last revision: 2/3/93 */


static void drift_kepu(const float dt, const float r0, const float mu, const float alpha, const float u, float &fp, float &c1, float &c2, float &c3, int &iflg)
{
	float s, st, fo, fn;

	drift_kepu_guess(dt, r0, mu, alpha, u, s);

	st = s;

	// store initial guess for possible use later in
	// laguerre's method, in case newton's method fails.
	drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);

	if(iflg != 0)
	{
		drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo);
		drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn);

		if(fabs(fo) < fabs(fn))
		{
			s = st;
		}

		drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);
	}
}


/*c********************************************************************#
c                  DRIFT_KEPMD
c********************************************************************#
c  Subroutine for solving kepler's equation in difference form for an
c  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
c  for the criteria.
c  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
c  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
c  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
c
c	Input:
c	    dm		==> increment in mean anomaly M (real*8 scalar)
c	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)
c
c       Output:
c            x          ==> solution to Kepler's difference eqn (real*8 scalar)
c            s,c        ==> sin and cosine of x (real*8 scalars)
c */

static void drift_kepmd(const float dm, const float es, const float ec, float &x, float &s, float &c)
{
	float 
		A0 = 39916800.0f, 
		A1 = 6652800.0f, 
		A2 = 332640.0f, 
		A3 = 7920.0f, 
		A4 = 110.0f;

	float dx, fac1, fac2, q, y, f, fp, fpp, fppp;

	// calculate initial guess for root
	fac1 = 1.0f / (1.0f - ec);
	q = fac1 * dm;
	fac2 = es * es * fac1 - ec / 3.0f;
	x = q * (1.0f - 0.5f * fac1 * q * (es - q * fac2));

	// excellent approx to sin and cos of x for small x
	y = x * x;
	s = x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
	c = sqrtf(1.0f - s * s);

	// computer better value for the root using quartic Newton method
	f = x - ec * s + es * (1.0f - c) - dm;
	fp = 1.0f - ec * c + es * s;
	fpp = ec * s + es * c;
	fppp = ec * c - es * s;
	dx = -f / fp;
	dx = -f / (fp + 0.5f * dx * fpp);
	dx = -f / (fp + 0.5f * dx * fpp + 0.16666666666666666f * dx * dx * fppp);
	x = x + dx;

	// excellent approx to sin and cos of x for small x
	y = x * x;
	s = x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
	c = sqrtf(1.0f - s * s);
}


/*c*************************************************************************
c                        DRIFT_DAN.F
c*************************************************************************
c This subroutine does the Danby and decides which vbles to use
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 x0,y0,z0         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx0,vy0,vz0      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dt0            ==>  time step
c             Output:
c                 x0,y0,z0         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx0,vy0,vz0      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 iflg             ==>  integer flag (zero if satisfactory)
c					      (non-zero if nonconvergence)
c
c Authors:  Hal Levison & Martin Duncan  
c Date:    2/10/93
c Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged */

static void drift_dan(const float mu, float &x0, float &y0, float &z0, float &vx0, float &vy0, float &vz0, const float dt0, int &iflg)
{
	float x,y,z,vx,vy,vz,dt;
	float f,g,fdot,c1,c2;
	float c3,gdot;
	float u,alpha,fp,r0,v0s;
	float a,asq,en;
	float dm,ec,es,esq,xkep;
	float fchk,s,c;

	// Set dt = dt0 to be sure timestep is not altered while solving
	// for new coords.
	dt = dt0;
	iflg = 0;

	r0 = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);
	v0s = vx0*vx0 + vy0 * vy0 + vz0 * vz0;
	u = x0 * vx0 + y0 * vy0 + z0 * vz0;
	alpha = 2.0f * mu / r0 - v0s;

	if(alpha > 0.0f)
	{
		a = mu / alpha;
		asq = a * a;
		en = sqrtf(mu / (a * asq));
		ec = 1.0f - r0 / a;
		es = u / (en * asq);
		esq = ec * ec + es * es;
		dm = dt * en - int(dt * en/TWOPI)*TWOPI;
		dt = dm / en;
		if(dm * dm > 0.16f || esq > 0.36f)
		{
			drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);

			if(iflg == 0)
			{
				f = 1.0f - (mu / r0) * c2;
				g = dt - mu * c3;
				fdot = -(mu / (fp * r0)) * c1;
				gdot = 1.0f - (mu / fp) * c2;

				x = x0 * f + vx0 * g;
				y = y0 * f + vy0 * g;
				z = z0 * f + vz0 * g;
				vx = x0 * fdot + vx0 * gdot;
				vy = y0 * fdot + vy0 * gdot;
				vz = z0 * fdot + vz0 * gdot;

				x0 = x;
				z0 = y;
				z0 = z;
				vx0 = vx;
				vy0 = vy;
				vz0 = vz;
			}

			return;
		}

		if(esq * dm * dm < 0.0016f)
		{
			drift_kepmd(dm, es, ec, xkep, s, c);
			fchk = (xkep - ec * s + es * (1.0f - c) - dm);

			if(fchk * fchk > DANBYB)
			{
				iflg = 1;
				return;
			}

			fp = 1.0f - ec * c + es * s;
			f = (a / r0) * (c - 1.0f) + 1.0f;
			g = dt + (s - xkep) / en;
			fdot = - (a / (r0 * fp)) * en * s;
			gdot = (c - 1.0f) / fp + 1.0f;

			x = x0 * f + vx0 * g;
			y = y0 * f + vy0 * g;
			z = z0 * f + vz0 * g;
			vx = x0 * fdot + vx0 * gdot;
			vy = y0 * fdot + vy0 * gdot;
			vz = z0 * fdot + vz0 * gdot;

			x0 = x;
			z0 = y;
			z0 = z;
			vx0 = vx;
			vy0 = vy;
			vz0 = vz;

			iflg = 0;
			return;
		}
	}

	drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);

	if(iflg == 0)
	{
		f = 1.0f - (mu / r0) * c2;
		g = dt - mu * c3;
		fdot = -(mu / (fp * r0)) * c1;
		gdot = 1.0f - (mu / fp) * c2;

		x = x0 * f + vx0 * g;
		y = y0 * f + vy0 * g;
		z = z0 * f + vz0 * g;
		vx = x0 * fdot + vx0 * gdot;
		vy = y0 * fdot + vy0 * gdot;
		vz = z0 * fdot + vz0 * gdot;

		x0 = x;
		z0 = y;
		z0 = z;
		vx0 = vx;
		vy0 = vy;
		vz0 = vz;
	}
}


/*c*************************************************************************
c                        DRIFT_ONE.F
c*************************************************************************
c This subroutine does the danby-type drift for one particle, using 
c appropriate vbles and redoing a drift if the accuracy is too poor 
c (as flagged by the integer iflg).
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mu            ==>  mass of central body (real scalar) 
c                 x,y,z         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx,vy,vz      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dt            ==>  time step
c             Output:
c                 x,y,z         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx,vy,vz      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 iflg          ==>  integer (zero for successful step)
c
c Authors:  Hal Levison & Martin Duncan 
c Date:    2/10/93
c Last revision: 2/10/93  */


static void drift_one(const float mu, float &x, float &y, float &z, float vx, float vy, float vz, const float dt, int &iflg)
{
	float dttmp;

	drift_dan(mu, x, y, z, vx, vy, vz, dt, iflg);

	if(iflg != 0)
	{
		for(int i = 0; i < 10; ++i)
		{
			dttmp = dt / 10.0f;
			drift_dan(mu, x, y, z, vx, vy, vz, dttmp, iflg);

			if(iflg != 0)
				return;
		}
	}
}

/*c*************************************************************************
c                        DRIFT_TP.F
c*************************************************************************
c This subroutine loops thorugh the TEST particles and calls the danby routine
c
c             INPUT:
c                 ntp             ==>  number of test particles (int scalar)
c                 msun            ==>  mass of the sun (real scalar)
c                 xjt,yjt,zjt     ==>  initial position in jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt  ==>  initial position in jacobi coord 
c                                      (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt              ==>  time step
c             OUTPUT:
c                 xjt,yjt,zjt     ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxjt,vyjt,vzjt  ==>  final position in jacobi coord 
c                                       (real arrays)
c
c Authors:  Hal Levison 
c Date:    2/18/93
c Last revision: */

void drift_tp(const size_t ntp, const float mass_sun, vec3_scalar_t &h, const vec3_scalar_t &vh, const float dt, int **istat)
{
	int iflg;

	for(size_t j=0; j < ntp; ++j)
	{
		if(istat[j][1] == 0)
		{
			drift_one(mass_sun, h.x[j], h.y[j], h.z[j], vh.x[j], vh.y[j], vh.z[j], dt, iflg);
			if(iflg != 0)
			{
				istat[j][1] = 1;
				istat[j][2] = -1;
			}
		}
	}
}