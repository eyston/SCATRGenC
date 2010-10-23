#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

#define iterations 1
#define N 5

void coord_h2j(int nbod, float mass[], float xh[], float yh[], float zh[],
									   float vxh[], float vyh[], float vzh[],
									   float xj[], float yj[], float zj[],
									   float vxj[], float vyj[], float vzj[])
{
	float *eta = new float[nbod];
	float sumx, sumy, sumz, sumvx, sumvy, sumvz;
	float capx, capy, capz, capvx, capvy, capvz;
	
	eta[0] = mass[0];
	
	for(int i=1; i<nbod; i++)
	{
		eta[i] = eta[i-1] + mass[i];
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

int main()
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
	
	/*ifstream input_file_s;
	input_file_s.open("param.in");
	
	input_file_s >> t0 >> tstop >> dt >> tinc;
	cout << t0 << "," << tstop << "," << dt << "," << tinc << endl;
	input_file_s.close();*/
	
	FILE *input_file;
	input_file = fopen("param.in", "r");
	
	fscanf(input_file, "%f %f %f %f\n", &t0, &tstop, &dt, &tinc);
	fscanf(input_file, "%f %f\n", &dtout, &dtdump);
	fscanf(input_file, "%c %c %c %c %c %c\n", &lflag_c[0], &lflag_c[1], &lflag_c[2], &lflag_c[3], &lflag_c[4], &lflag_c[5]);
	fscanf(input_file, "%f %f %f %f %c\n", &rmin, &rmax, &rmaxu, &qmin, &lclose_c);
	fscanf(input_file, "%f\n", &rcrit);
	fscanf(input_file, "%s\n", &output_file_name);
	fscanf(input_file, "%s\n", &frame);
	fscanf(input_file, "%s\n", &fopenstat);

	fclose(input_file);
	
	printf("%f, %f, %f, %f\n", t0, tstop, dt, tinc);
	printf("%f, %f\n", dtout, dtdump);
	printf("%c, %c, %c, %c, %c, %c\n", lflag_c[0], lflag_c[1], lflag_c[2], lflag_c[3], lflag_c[4], lflag_c[5]);
	
	for(int i=0; i<6; i++)
	{
		lflag[i] = lflag_c[i] == 'T';
	}
	
	printf("%d, %d, %d, %d, %d, %d\n", lflag[0], lflag[1], lflag[2], lflag[3], lflag[4], lflag[5]);
	
	lclose = lclose_c == 'T';

	printf("%f, %f, %f, %f, %d\n", rmin, rmax, rmaxu, qmin, lclose);
	printf("%f\n", rcrit);
	printf("%s\n", output_file_name);
	printf("%s\n", frame);
	printf("%s\n", fopenstat);
	
/*	cout << "load n-body values from file (xh, vh)" << endl;
	cout << "convert n-body helio-centric coords to jacobian (xj)" << endl;
	cout << "calculate initial acceleration of n-bodies (ah)" << endl;
	
	for(int i=0; i < iterations; i++)
	{
		cout << " -- apply a heliocentric half-kick to n-bodies (xh)" << endl;
		cout << " -- convert heliocentric velocities to jacobian velocities (vj)" << endl;
		cout << " -- apply a drift in jacobian coordinates (xj)" << endl;
		cout << " -- calculate new heliocentric velocities and positions (xh, vh)" << endl;
		cout << " -- calculate the acceleration of n-bodies in heliocentric coords (ah)" << endl;
		cout << " -- apply a heliocentric half-kick to n-bodies (xh)" << endl;
	} */
	
	return 0;
}