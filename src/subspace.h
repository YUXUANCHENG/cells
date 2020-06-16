#ifndef SUBSPACE_H
#define SUBSPACE_H

#include "deformableParticles2D.h"
#include "cellPacking2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>

using namespace std;

class cellPacking2D;

class subspace{

protected:

	// resident cells and cashed cells
	vector<deformableParticles2D*> resident_cells;
	vector<deformableParticles2D*> cashed_cells;

	// list indicates near boundary cells that migrate to neighbor boxes
	stack<int> migrate_out_list;
	stack<int> migrate_out_destination;

	// pointer to the whole system (cell_group)
	cellPacking2D * pointer_to_system;

	// which box is this
	int box_id;

	vector<double> L;
	vector<int> N_systems;

	int NDIM = 2;

	double sigmaXX = 0.0;
	double sigmaXY = 0.0;
	double sigmaYX = 0.0;
	double sigmaYY = 0.0;
	double dt0;
	double PI = 4 * atan(1);
	
	// indicate what fraction of the system size will be cashed
	double cashed_fraction;

public:

	void initialize(cellPacking2D * const & pointer, vector<double> const & L, vector<int> const& N_systems,int box_id, double const & dt0) {
		pointer_to_system = pointer;
		this->L = L;
		this->box_id = box_id;
		this->N_systems = N_systems;
		this->dt0 = dt0;
	};

	void cashe_out(int direction);
	void reset_cashe();
	int neighbor_box(int direction, int upper_lower);
	double find_boundary(int direction, int upper_lower);
	void migrate_out();

	void cashe_in(vector<deformableParticles2D*>& cash_list);
	void migrate_in(deformableParticles2D* const & migration);

	void calculateForces_insub();
	void activityCOM_brownian_insub(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void conserve_momentum();


};







#endif