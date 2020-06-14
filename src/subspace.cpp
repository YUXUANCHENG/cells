// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"
#include "random"


using namespace std;

void cellPacking2D::split_into_subspace() {
	subsystem = new subspace[16];
};


void cellPacking2D::cash_into(int i, vector<deformableParticles2D*>& cash_list) {
	subsystem[i].cash_in(cash_list);
};

void cellPacking2D::migrate_into(int i, deformableParticles2D* migration) {
	subsystem[i].migrate_in(migration);
};


void subspace::cash_in(vector<deformableParticles2D*>& cash_list) {
	cashed_cells.insert(cashed_cells.end(),cash_list.begin(),cash_list.end());
};

void subspace::migrate_in(deformableParticles2D*& migration) {
	resident_cells.push_back(migration);
};

void subspace::cash_out() {


};

void subspace::migrate_out() {


};

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


	// loop over cells and cell pairs, calculate shape and interaction forces
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

		// forces between resident cell and cashed cell
		for (ck = 0; cj < cashed_cells.size(); cj++) {
			inContact = resident_cells[ci]->vertexForce_cashed(*resident_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
		}


		// forces on vertices due to shape
		resident_cells[ci]->shapeForces();
	}


};



void subspace::activityCOM_brownian_insub(double T, double v0, double Dr, double vtau, double t_scale, int frames, double scaled_v) {


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
	//double scaled_v = scale_v(v0);

	// Reset velocity
	for (ci = 0; ci < resident_cells.size(); ci++) {
		for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
			for (d = 0; d < NDIM; d++)
				resident_cells[ci]->setVVel(vi, d, 0.0);
		}
	}

	// Active brownian MD
	for (t = 0.0; t < T; t += dt0 * t_scale) {
		// do verlet update
		for (ci = 0; ci < resident_cells.size(); ci++) {
			resident_cells[ci]->verletPositionUpdate(dt0 * t_scale);
			resident_cells[ci]->updateCPos();
		}

		// reset contacts before force calculation
		//resetContacts();

		// calculate forces
		calculateForces_insub();

		// update velocities
		for (ci = 0; ci < resident_cells.size(); ci++) {
			std::normal_distribution<double> dist(0, 1);

			// get random number with normal distribution using gen as random source
			random_angle = dist(gen);
			resident_cells[ci]->activeVerletVelocityUpdateCOM_brownian(dt0 * t_scale, Dr, random_angle, scaled_v);
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

	for (int ci = 0; ci < resident_cells.size(); ci++) {
		system_mass += resident_cells[ci]->getNV() * PI * pow(0.5 * resident_cells[ci]->getdel() * resident_cells[ci]->getl0(), 2);
	}

	for (int ci = 0; ci < resident_cells.size(); ci++) {
		for (int d = 0; d < NDIM; d++) {
			system_p[d] += resident_cells[ci]->momentum(d);
		}
	}

	for (int ci = 0; ci < resident_cells.size(); ci++) {
		for (int vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
			for (int d = 0; d < NDIM; d++) {

				v_temp = resident_cells[ci]->vvel(vi, d) - system_p[d] / system_mass;
				resident_cells[ci]->setVVel(vi, d, v_temp);
			}
		}
	}


}