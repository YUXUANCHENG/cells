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

void cellPacking2D::initialize_subsystems(int N_x, int N_y) {

	N_systems.push_back(N_x);
	N_systems.push_back(N_y);
	split_into_subspace();

}

void cellPacking2D::paralell_activityCOM_brownian(double T, double v0, double Dr, double vtau, double t_scale, int frames) {

omp_set_num_threads(N_systems.at(0) * N_systems.at(1));
#pragma omp parallel
	{
		int i = omp_get_thread_num();
		(subsystem[i]).activityCOM_brownian_insub(T, v0, Dr, vtau, t_scale, frames);
	}
}


// compress isotropically to fixed packing fraction
void cellPacking2D::paralell_qsIsoCompression(double phiTarget, double deltaPhi, double Ftol) {
	// local variables
	double phi0, phiNew, dphi, Fcheck, Kcheck;
	int NSTEPS, k;

	// get initial packing fraction
	phi = packingFraction();
	phi0 = phi;

	// determine number of steps to target
	NSTEPS = floor(abs((phiTarget - phi0)) / deltaPhi);
	if (NSTEPS == 0)
		NSTEPS = 1;

	// update new dphi to make steps even
	dphi = (phiTarget - phi) / NSTEPS;

	// iterator
	k = 0;

	// loop until phi is the correct value
	while (k < NSTEPS) {
		// update iterator
		k++;

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression with NSTEPS = " << NSTEPS << " and dphi = " << dphi << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* NSTEPS 		= " << NSTEPS << endl;
		cout << "	* dphi 			= " << dphi << endl << endl;
		cout << "	AFTER LAST MINIMIZATION:" << endl;
		cout << "	* phi 			= " << phi << endl;
		//cout << "	* Fcheck 		= " << Fcheck << endl;
		//cout << "	* Kcheck 		= " << Kcheck << endl;
		cout << endl << endl;

		// increase packing fraction to new phi value
		phiNew = phi0 + k * dphi;
		setPackingFraction(phiNew);

		// calculate phi before minimization
		phi = packingFraction();

		// relax shapes (energies calculated in relax function)

		paralell_fireMinimizeF(Ftol, Fcheck, Kcheck);
		if (k % 10 == 0) {
			printJammedConfig_yc();
			printCalA();
			printContact();
		}
	}

	while (phi < phiTarget) {

		phiNew = phi + deltaPhi;
		setPackingFraction(phiNew);

		// calculate phi before minimization
		phi = packingFraction();

		// relax shapes (energies calculated in relax function)
		paralell_fireMinimizeF(Ftol, Fcheck, Kcheck);

	}

}


// FIRE 2.0 force minimization with backstepping
void cellPacking2D::paralell_fireMinimizeF(double Ftol, double& Fcheck, double& Kcheck) {
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 = 0.3;
	const double finc = 1.1;
	const double fdec = 0.5;
	const double falpha = 0.99;
	const double dtmax = 10 * dt0;
	const double dtmin = 1e-8 * dt0;
	const double Trescale = 1e-8 * NCELLS;
	const int NMIN = 20;
	const int NNEGMAX = 2000;
	const int NDELAY = 1000;
	int npPos = 0;
	int npNeg = 0;
	int npPMIN = 0;
	double alpha = 0.0;
	double t = 0.0;
	double P = 0.0;

	// local variables
	int ci, vi, d, k, kmax;
	double vstarnrm, fstarnrm, vtmp, ftmp;
	double K, F, Pcheck;
	double xold, xnew, vold, vnew, fold;

	// variable to test for potential energy minimization
	bool converged = false;

	// reset time step
	dt = dt0;

	// initialize forces
	calculateForces();

	// rescale velocities
	//rescaleVelocities(Trescale);

	// norm of total force vector, kinetic energy
	F = forceRMS();
	K = totalKineticEnergy();
	Pcheck = 0.5 * (sigmaXX + sigmaYY) / (NCELLS * L.at(0) * L.at(1));

	// scale P and K for convergence checking
	Fcheck = F;
	Kcheck = K / NCELLS;

	omp_set_num_threads(N_systems.at(0) * N_systems.at(1));
#pragma omp parallel
	{
		int i = omp_get_thread_num();
		(subsystem[i]).fireMinimizeF_insub(Ftol, Fcheck, Kcheck, P, vstarnrm, fstarnrm, converged);
	}
}



// FIRE 2.0 force minimization with backstepping
void subspace::fireMinimizeF_insub(double Ftol, double& Fcheck, double& Kcheck, double & P, double & vstarnrm, double & fstarnrm, bool & converged) {
	// HARD CODE IN FIRE PARAMETERS
	const double alpha0 = 0.3;
	const double finc = 1.1;
	const double fdec = 0.5;
	const double falpha = 0.99;
	const double dtmax = 10 * dt0;
	const double dtmin = 1e-8 * dt0;
	const int NMIN = 20;
	const int NNEGMAX = 2000;
	const int NDELAY = 1000;
	int npPos = 0;
	int npNeg = 0;
	int npPMIN = 0;
	double alpha = 0.0;
	double t = 0.0;
	double local_P = 0.0;

	// local variables
	int ci, vi, d, k, kmax;
	double local_vstarnrm, local_fstarnrm, vtmp, ftmp;
	double K, F, Pcheck;
	double xold, xnew, vold, vnew, fold;

	// reset time step
	double dt = dt0;

	int NCELLS = pointer_to_system->getNCELLS();

	// calculate cashed fraction
	double spacing = L.at(0) / N_systems[0];
	cashed_fraction = pointer_to_system->scale_v(2) / spacing;

	// iterate until system converged
	kmax = 1e6;
	for (k = 0; k < kmax; k++) {

		// Step 1. calculate P and norms

		// be careful about synchronization

#pragma omp master
		{
			P = 0.0;
			vstarnrm = 0.0;
			fstarnrm = 0.0;
			Fcheck = 0.0;
			Kcheck = 0.0;
		}
#pragma omp barrier
		local_P = 0;
		local_vstarnrm = 0.0;
		local_fstarnrm = 0.0;
		for (ci = 0; ci < resident_cells.size(); ci++) {
			for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
				for (d = 0; d < NDIM; d++) {
					// get tmp variables
					ftmp = resident_cells[ci]->vforce(vi, d);
					vtmp = resident_cells[ci]->vvel(vi, d);

					// calculate based on all vertices on all cells
					local_P += ftmp * vtmp;
					local_vstarnrm += vtmp * vtmp;
					local_fstarnrm += ftmp * ftmp;
				}
			}
		}

		// be careful about synchronization

#pragma omp critical
		{
			P += local_P;
			vstarnrm += local_vstarnrm;
			fstarnrm += local_fstarnrm;
		}
		// only master thread
		// get norms
#pragma omp barrier
#pragma omp master
		{
			vstarnrm = sqrt(vstarnrm);
			fstarnrm = sqrt(fstarnrm);

			// output some information to console
			if (k % pointer_to_system->getNPRINT() == 0) {
				cout << "===================================================" << endl << endl;
				cout << " 	FIRE MINIMIZATION, k = " << k << endl << endl;
				cout << "===================================================" << endl;
				cout << "	* Run data:" << endl;
				cout << "	* Kcheck 	= " << Kcheck << endl;
				cout << "	* Fcheck 	= " << Fcheck << endl;
				//cout << "	* Pcheck 	= " << Pcheck << endl;
				cout << "	* phi 		= " << pointer_to_system->getphi() << endl;
				cout << "	* dt 		= " << dt << endl;
				cout << "	* alpha 	= " << alpha << endl;
				cout << "	* P 		= " << P << endl;
				cout << "	* Pdir 		= " << P / (vstarnrm * fstarnrm) << endl;
				cout << endl << endl;
			}
		}
#pragma omp barrier
		// Step 2. Adjust simulation based on net motion of system
		if (P > 0) {
			// increment pos counter
			npPos++;

			// reset neg counter
			npNeg = 0;

			// alter sim if enough positive steps taken
			if (npPos > NMIN) {
				// change time step
				if (dt * finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else {
			// reset pos counter
			npPos = 0;

			// rest neg counter
			npNeg++;

			// check for stuck sim
			if (npNeg > NNEGMAX)
				break;

			// decrease time step if past initial delay
			if (k > NDELAY) {
				// decrease time step 
				if (dt * fdec > dtmin)
					dt *= fdec;

				// change alpha
				alpha = alpha0;
			}

			// take half step backwards
			for (ci = 0; ci < resident_cells.size(); ci++) {
				for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
					for (d = 0; d < NDIM; d++)
						resident_cells[ci]->setVPos(vi, d, resident_cells[ci]->vpos(vi, d) - 0.5 * dt * resident_cells[ci]->vvel(vi, d));
				}
			}

			// reset velocities to 0
			for (ci = 0; ci < resident_cells.size(); ci++) {
				for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
					for (d = 0; d < NDIM; d++)
						resident_cells[ci]->setVVel(vi, d, 0.0);
				}
			}
		}

		double tmp1, tmp2;
		// update velocities if forces are acting
		if (fstarnrm > 0) {
			for (ci = 0; ci < resident_cells.size(); ci++) {
				for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
					for (d = 0; d < NDIM; d++) {
						vtmp = (1 - alpha) * resident_cells[ci]->vvel(vi, d) + alpha * (resident_cells[ci]->vforce(vi, d) / fstarnrm) * vstarnrm;
						resident_cells[ci]->setVVel(vi, d, vtmp);
					}
				}
			}
		}

		// VV update in FIRE 2.0: position update
		for (ci = 0; ci < resident_cells.size(); ci++) {
			resident_cells[ci]->verletPositionUpdate(dt);
			resident_cells[ci]->updateCPos();
		}

		// perform cashing and migration
#pragma omp barrier
// To avoid deadlock, migrate sequencially
#pragma omp critical
		{
			migrate_out();
		}

#pragma omp barrier

#pragma omp critical
		{
			reset_cashe();
		}
		// seems like there is no deadlock when cashe at x direction
#pragma omp barrier

#pragma omp critical
		{
			cashe_out(0);
		}

		// To avoid deadlock, only update odd or even box id
#pragma omp barrier

#pragma omp critical
		{
			cashe_out(1);
		}
#pragma omp barrier

		// calculate forces
		calculateForces_insub();

		// VV update in FIRE 2.0: Velocity update 2
		for (ci = 0; ci < resident_cells.size(); ci++)
			resident_cells[ci]->verletVelocityUpdate(dt);

		// update t
		t += dt;

		// track energy and forces
		F = forceRMS_insub();
		K = totalKineticEnergy_insub();
		//Pcheck = 0.5 * (sigmaXX + sigmaYY) / (pointer_to_system->getNCELLS() * L.at(0) * L.at(1));

		// scale P and K for convergence checking
		// be careful about synchronization
#pragma omp critical
		{
			Fcheck += F;
			Kcheck += K / NCELLS;
		}
#pragma omp barrier
		// update if Fcheck under tol
		if (abs(Fcheck) < Ftol)
			npPMIN++;
		else
			npPMIN = 0;

		// check that P is not crazy
		if (abs(P) > 800) {
			cout << "	ERROR: P = " << P << ", ending." << endl;
			cout << "	** Kcheck = " << Kcheck << endl;
			cout << "	** Fcheck = " << Fcheck << endl;
			//cout << "	** Pcheck = " << Pcheck << endl;
			cout << "	** k = " << k << ", t = " << t << endl;

			/*
			// print minimized config, energy and contact network
			if (packingPrintObject.is_open()) {
				cout << "	* Printing vetex positions to file" << endl;
				printSystemPositions();
			}

			if (energyPrintObject.is_open()) {
				cout << "	* Printing cell energy to file" << endl;
				printSystemEnergy(k);
			}
			*/
			exit(1);
			
		}

		// check for convergence
		// be careful about synchronization, only master thread update parameter converged

#pragma omp master
		{
			converged = (abs(Fcheck) < Ftol && npPMIN > NMIN);

			if (converged) {
				cout << "	** FIRE has converged!" << endl;
				cout << "	** Fcheck = " << Fcheck << endl;
				cout << "	** Kcheck = " << Kcheck << endl;
				//cout << "	** Pcheck = " << Pcheck << endl;
				cout << "	** k = " << k << ", t = " << t << endl;
				cout << "	** Breaking out of FIRE protocol." << endl;

				/*
				// print minimized config, energy and contact network
				if (packingPrintObject.is_open()) {
				cout << "	* Printing vetex positions to file" << endl;
				pointer_to_system->printSystemPositions();
				}
				
				if (energyPrintObject.is_open()) {
					cout << "	* Printing cell energy to file" << endl;
					printSystemEnergy(k);
				}

				break;
				*/
			}
		}
#pragma omp barrier
		if (converged) {
			break;
		}

	}

	// reset dt to be original value before ending function
	dt = dt0;

	// if no convergence, just stop
	if (k == kmax) {
		cout << "	** ERROR: FIRE not converged in kmax = " << kmax << " force evaluations, ending code" << endl;
		exit(1);
	}
}


double subspace::forceRMS_insub() {
	int ci, vi, d;
	int NVTOTAL = 0;
	double frms = 0.0;

	// loop over forces, calc total force norm
	for (ci = 0; ci < resident_cells.size(); ci++) {
		NVTOTAL += resident_cells[ci]->getNV();
		for (vi = 0; vi < resident_cells[ci]->getNV(); vi++) {
			for (d = 0; d < NDIM; d++)
				frms += pow(resident_cells[ci]->vforce(vi, d), 2);
		}
	}

	// get force scale
	frms = sqrt(frms) / (NDIM * NVTOTAL);

	// return
	return frms;
}


// Calculate total kinetic energy in system
double subspace::totalKineticEnergy_insub() {
	// local variables
	int ci;
	double val = 0.0;

	// loop over cells, add kinetic energy
	for (ci = 0; ci < resident_cells.size(); ci++)
		val += resident_cells[ci]->totalKineticEnergy();

	// return value
	return val;
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
					inContact = resident_cells[ci]->vertexForce_cashed(*cashed_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
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
	double phi = 0.0;

	// Cal print frequency
	int print_frequency = floor(T / (dt0 * t_scale * frames));

	// random device class instance, source of 'true' randomness for initializing random seed
	std::random_device rd;

	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen(rd());

	double random_angle;

	// Scale velocity by avg cell radius
	double scaled_v = pointer_to_system->scale_v(v0);

	// calculate cashed fraction
	double spacing = L.at(0) / N_systems[0];
	cashed_fraction = pointer_to_system->scale_v(2) / spacing;

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
		/*
#pragma omp barrier
		// To avoid deadlock, migrate sequencially
#pragma omp critical
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
	*/

#pragma omp barrier
	// To avoid deadlock, migrate sequencially
#pragma omp critical
		{
			migrate_out();
		}

#pragma omp barrier

#pragma omp critical
		{
		reset_cashe();
		}
		// seems like there is no deadlock when cashe at x direction
#pragma omp barrier

#pragma omp critical
		{
			cashe_out(0);
		}

		// To avoid deadlock, only update odd or even box id
#pragma omp barrier

#pragma omp critical
		{
			cashe_out(1);
		}
#pragma omp barrier



#pragma omp master
		{
		//phi = pointer_to_system->packingFraction();

		if (count % print_frequency == 0) {
			pointer_to_system->printJammedConfig_yc();
			//pointer_to_system->phiPrintObject << phi << endl;
			pointer_to_system->printCalA();
			//pointer_to_system->printContact();
			pointer_to_system->printV();
		}
		}

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