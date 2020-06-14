#ifndef SUBSPACE_H
#define SUBSPACE_H

#include "deformableParticles2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

class subspace{

protected:
	vector<deformableParticles2D*> resident_cells;
	vector<deformableParticles2D*> cashed_cells;

	vector<deformableParticles2D*> cash_out_list;
	stack<deformableParticles2D*> migrate_out_list;

	int NDIM = 2;

	double sigmaXX = 0.0;
	double sigmaXY = 0.0;
	double sigmaYX = 0.0;
	double sigmaYY = 0.0;

	double dt0;
	double PI = 4 * atan(1);

public:
	subspace();
	int look_for_new_box(deformableParticles2D*& cell);
	void cash_out();
	void migrate_out();

	void cash_in(vector<deformableParticles2D*>& cash_list);
	void migrate_in(deformableParticles2D*& migration);

	void calculateForces_insub();
	void activityCOM_brownian_insub(double T, double v0, double Dr, double vtau, double t_scale, int frames, double scaled_v);
	void conserve_momentum()


};







#endif