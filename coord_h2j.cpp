#include "coord_h2j.h"

void coord_h2j(int nbod, NBody bodies[])
{
	float *eta = new float[nbod];
	float sumx, sumy, sumz, sumvx, sumvy, sumvz;
	float capx, capy, capz, capvx, capvy, capvz;
	
	eta[0] = mass[0];
	
	for(int i=1; i<nbod; i++)
	{
		eta[i] = eta[i-1] + bodies[i].mass;
	}
	
	xj[0] = 0.0;
	yj[0] = 0.0;
	zj[0] = 0.0;
	vxj[0] = 0.0;
	vyj[0] = 0.0;
	vzj[0] = 0.0;
	
	xj[1] = xh[1];
	yj[1] = yh[1];
	zj[1] = zh[1];
	vxj[1] = vxh[1];
	vyj[1] = vyh[1];
	vzj[1] = vzh[1];
	
	sumx = mass[1]*xh[1];
	sumy = mass[1]*yh[1];
	sumz = mass[1]*zh[1];
	sumvx = mass[1]*vxh[1];
	sumvy = mass[1]*vyh[1];
	sumvz = mass[1]*vzh[1];
	
	capx = sumx/eta[1];
	capy = sumy/eta[1];
	capz = sumz/eta[1];
	capvx = sumvx/eta[1];
	capvy = sumvy/eta[1];
	capvz = sumvz/eta[1];
	
	for(int i=2; i<nbod; i++)
	{
		xj[i] = xh[i] - capx;
		yj[i] = yh[i] - capy;
		zj[i] = zh[i] - capz;
		vxj[i] = vxh[i] - capvx;
		vyj[i] = vyh[i] - capvy;
		vzj[i] = vzh[i] - capvz;
		
		if(i < (nbod-1))
		{
			sumx = sumx + mass[i]*xh[i];
			sumy = sumy + mass[i]*yh[i];
			sumz = sumz + mass[i]*zh[i];
			sumvx = sumvx + mass[i]*vxh[i];
			sumvy = sumvy + mass[i]*vyh[i];
			sumvz = sumvz + mass[i]*vzh[i];
			
			capx = sumx/eta[i];
			capy = sumy/eta[i];
			capz = sumz/eta[i];
			capvx = sumvx/eta[i];
			capvy = sumvy/eta[i];
			capvz = sumvz/eta[i];
		}
	}

	delete[] eta;
}
