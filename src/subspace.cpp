// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"
#include "random"
#include <cmath>

using namespace std;

// split the packing system into smaller subsystems
void cellPacking2D::split_into_subspace() {
	int box;
	subsystem = new subspace[N_systems[0] * N_systems[1]];

	for (int i = 0; i < N_systems[0] * N_systems[1]; i++) {
		(subsystem[i]).initialize(this, L, N_systems, i, dt0);
	}
	
	// assign cells into subsystems
	for (int ci = 0; ci < NCELLS; ci++) {
		box = look_for_new_box(cell(ci));
		migrate_into(box, &(cell(ci)));
	}
};

// cashe list send to subsystems
void cellPacking2D::cashe_into(int i, vector<deformableParticles2D*>& cash_list) {
	subsystem[i].cashe_in(cash_list);
};

// migrate cells into subsystems
void cellPacking2D::migrate_into(int i, deformableParticles2D* const & migration) {
	subsystem[i].migrate_in(migration);
};

// figure out which box the cells belong to
int cellPacking2D::look_for_new_box(deformableParticles2D & cell) {
	int box_id = 0;
	int x_id = 0;
	int y_id = 0;

	x_id = floor(cell.cpos(0)/(L.at(0)/N_systems[0]));
	y_id = floor(cell.cpos(1) / (L.at(1) / N_systems[1]));
	box_id = y_id * N_systems[0] + x_id;

	return box_id;
}

void cellPacking2D::activityCOM_brownian_subsystem(int N_x, int N_y, double T, double v0, double Dr, double vtau, double t_scale, int frames) {

	N_systems.push_back(N_x);
	N_systems.push_back(N_y);
	split_into_subspace();

#pragma omp parallel for
	for(int i = 0; i < N_systems[0] * N_systems[1]; i++)
		(subsystem[i]).activityCOM_brownian_insub(T, v0, Dr, vtau, t_scale, frames);

}

// cashe cells into cashe list
void subspace::cashe_in(vector<deformableParticles2D*>& cash_list) {
	cashed_cells.insert(cashed_cells.end(),cash_list.begin(),cash_list.end());
};

// add migrated cells
void subspace::migrate_in(deformableParticles2D* const & migration) {
	resident_cells.push_back(migration);
};

// reset cashe system
void subspace::reset_cashe() {
	if (!cashed_cells.empty()) {
		cashed_cells.clear();
	}
}

// send cashe list 
void subspace::cashe_out(int direction) {

	// list indicates near boundary cells that need to be sent to neighbor boxes
	vector<deformableParticles2D*> cash_out_list_lower;
	vector<deformableParticles2D*> cash_out_list_upper;

	// find neighbor boxes
	int lower_index = neighbor_box(direction, -1);
	int upper_index = neighbor_box(direction, +1);

	double spacing = L.at(direction) / N_systems[direction];
	// find boundaries of the current box
	double lower_boundary = find_boundary(direction, -1);
	double upper_boundary = find_boundary(direction, +1);

	// check if resident cells are near boundary
	if (!resident_cells.empty()) {
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			if (resident_cells[ci]->cpos(direction) < lower_boundary + cashed_fraction * spacing &&
					resident_cells[ci]->cpos(direction) > lower_boundary)

				cash_out_list_lower.push_back(resident_cells[ci]);
			else if(resident_cells[ci]->cpos(direction) > upper_boundary - cashed_fraction * spacing &&
					resident_cells[ci]->cpos(direction) < upper_boundary)

				cash_out_list_upper.push_back(resident_cells[ci]);
		}
	}

	// check if cashed cells are near boundary, but only for y direction
	if (!cashed_cells.empty() && direction == 1) {
		for (int ci = 0; ci < cashed_cells.size(); ci++) {
			if (cashed_cells[ci]->cpos(direction) < lower_boundary + cashed_fraction * spacing &&
					cashed_cells[ci]->cpos(direction) > lower_boundary)

				cash_out_list_lower.push_back(cashed_cells[ci]);
			else if (cashed_cells[ci]->cpos(direction) > upper_boundary - cashed_fraction * spacing &&
					cashed_cells[ci]->cpos(direction) < upper_boundary)

				cash_out_list_upper.push_back(cashed_cells[ci]);
		}

	}

	// send to other boxes
	pointer_to_system->cashe_into(lower_index, cash_out_list_lower);
	pointer_to_system->cashe_into(upper_index, cash_out_list_upper);

};

// find neightor box
int subspace::neighbor_box(int direction, int upper_lower) {
	int current[2], neighbor_box_id;

	// current box id in x and y
	current[0] = box_id % N_systems[0];
	current[1] = floor(box_id / N_systems[0]);

	// neighbor box
	current[direction] += upper_lower;
	// periodic boundary
	current[direction] = (current[direction] + N_systems[direction]) % N_systems[direction];

	// neighbor id
	neighbor_box_id = current[0] + current[1] * N_systems[0];

	return neighbor_box_id;
}

double subspace::find_boundary(int direction, int upper_lower) {
	int current[2];
	double boundary;
	double spacing = L.at(direction) / N_systems[direction];

	// current box id in x and y
	current[0] = box_id % N_systems[0];
	current[1] = floor(box_id / N_systems[0]);

	// find boundary
	if (upper_lower == -1)
		boundary = current[direction] * spacing;
	else if(upper_lower == 1)
		boundary = (current[direction] + 1) * spacing;

	return boundary;
}



// migrate cells to other boxes and remove from current box
void subspace::migrate_out() {

	int new_box_index;
	int cell_index;
	deformableParticles2D* target_cell;

	// empty the stack
	while (!migrate_out_list.empty()) {
		migrate_out_list.pop();
	}

	while (!migrate_out_destination.empty()) {
		migrate_out_destination.pop();
	}

	// find the migration list
	for (int ci = 0; ci < resident_cells.size(); ci++) {
		new_box_index = pointer_to_system->look_for_new_box(*resident_cells.at(ci));
		if (new_box_index != box_id) {
			migrate_out_list.push(ci);
			migrate_out_destination.push(new_box_index);
		}
	}

		
	// migrate to other subsystems
	if (!migrate_out_list.empty()){for (int i = 0; i < migrate_out_list.size(); i++) {
			// migrate backwards, otherwise the list order would be changed
			cell_index = migrate_out_list.top();
			target_cell = resident_cells[cell_index];
			// find which subsystem to go
			new_box_index = migrate_out_destination.top();
			// migrate
			pointer_to_system->migrate_into(new_box_index, target_cell);
			// pop from list
			migrate_out_list.pop();
			migrate_out_destination.pop();
			// remove from resident list
			resident_cells.erase(resident_cells.begin() + cell_index);
		}
	}
};



// calculate forces 
void subspace::calculateForces_insub() {
	// local variables
	int ci, cj, ck, vi, d, dd, inContact;

	
	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;


	/*
	// reset contacts before force calculation
	resetContacts();
	Ncc = 0;
	Nvv = 0;
	*/


	// reset forces
	if (!resident_cells.empty()) {
	
		for (ci = 0; ci < resident_cells.size(); ci++) {
			// reset center of mass forces
			for (d = 0; d < NDIM; d++)
				resident_cells[ci]->setCForce(d, 0.0);

			// reset vertex forces and interaction energy
			for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
				// forces
				for (d = 0; d < NDIM; d++)
					resident_cells[ci]->setVForce(vi, d, 0.0);

				// energies
				resident_cells[ci]->setUInt(vi, 0.0);
			}
		}
	}


	// loop over cells and cell pairs, calculate shape and interaction forces
	if (!resident_cells.empty()) {
		for (ci = 0; ci < resident_cells.size(); ci++) {
			// forces between resident cells
			// loop over pairs, add info to contact matrix
			for (cj = ci + 1; cj < resident_cells.size(); cj++) {
				// calculate forces, add to number of vertex-vertex contacts
				inContact = resident_cells[ci]->vertexForce(*resident_cells[cj], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				/*
							if (inContact > 0) {
								// add to cell-cell contacts
								addContact(ci, cj);
								Ncc++;

								// increment vertex-vertex contacts
								Nvv += inContact;
							}
				*/
			}

			if (!cashed_cells.empty()) {
				// forces between resident cell and cashed cell
				for (ck = 0; ck < cashed_cells.size(); ck++) {
					inContact = resident_cells[ci]->vertexForce_cashed(*resident_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				}
			}

			// forces on vertices due to shape
			resident_cells[ci]->shapeForces();
		}
	}


};


// active brownian in subsystems
void subspace::activityCOM_brownian_insub(double T, double v0, double Dr, double vtau, double t_scale, int frames) {


	int ci, vi, d;
	int count = 0;
	double t = 0.0;

	// Cal print frequency
	int print_frequency = floor(T / (dt0 * t_scale * frames));

	// random device class instance, source of 'true' randomness for initializing random seed
	std::random_device rd;

	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen(rd());

	double random_angle;

	// Scale velocity by avg cell radius
	double scaled_v = pointer_to_system->scale_v(v0);

	// Reset velocity
	if (!resident_cells.empty()) {
		for (ci = 0; ci < resident_cells.size(); ci++) {
			for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
				for (d = 0; d < NDIM; d++)
					resident_cells[ci]->setVVel(vi, d, 0.0);
			}
		}
	}
	// Active brownian MD
	for (t = 0.0; t < T; t += dt0 * t_scale) {
		// do verlet update
		if (!resident_cells.empty()) {
			for (ci = 0; ci < resident_cells.size(); ci++) {
				resident_cells[ci]->verletPositionUpdate(dt0 * t_scale);
				resident_cells[ci]->updateCPos();
			}
		}
		// reset contacts before force calculation
		//resetContacts();

		//need to avoid deadlock
#pragma omp barrier
		// To avoid deadlock, migrate sequencially
#pragma omp criticle
		migrate_out();

#pragma omp barrier
		reset_cashe();

		// seems like there is no deadlock when cashe at x direction
#pragma omp barrier
		cashe_out(0);

		// To avoid deadlock, only update odd or even box id
#pragma omp barrier
		int y_id = floor(box_id / N_systems[0]);
		if (y_id % 2 == 0)
			cashe_out(1);
#pragma omp barrier
		if (y_id % 2 == 1)
			cashe_out(1);
#pragma omp barrier



		// calculate forces
		calculateForces_insub();

		// update velocities
		if (!resident_cells.empty()) {
			for (ci = 0; ci < resident_cells.size(); ci++) {
				std::normal_distribution<double> dist(0, 1);

				// get random number with normal distribution using gen as random source
				random_angle = dist(gen);
				resident_cells[ci]->activeVerletVelocityUpdateCOM_brownian(dt0 * t_scale, Dr, random_angle, scaled_v);
			}
		}
		count++;

		// Reset COM velocity
		conserve_momentum();

		/*
		phi = packingFraction();

		if (count % print_frequency == 0) {
			printJammedConfig_yc();
			phiPrintObject << phi << endl;
			printCalA();
			printContact();
			printV();
		}
		*/
	}



};


void subspace::conserve_momentum() {
	double factor = 0.0;
	double system_p[2] = { 0.0, 0.0 };
	double system_mass = 0.0;
	double v_temp = 0.0;
	if (!resident_cells.empty()) {
		// get systems mass
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			system_mass += resident_cells[ci]->getNV() * PI * pow(0.5 * resident_cells[ci]->getdel() * resident_cells[ci]->getl0(), 2);
		}

		//get system momentum
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			for (int d = 0; d < NDIM; d++) {
				system_p[d] += resident_cells[ci]->momentum(d);
			}
		}

		// substract COM velocity
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			for (int vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
				for (int d = 0; d < NDIM; d++) {

					v_temp = resident_cells[ci]->vvel(vi, d) - system_p[d] / system_mass;
					resident_cells[ci]->setVVel(vi, d, v_temp);
				}
			}
		}
	}

}