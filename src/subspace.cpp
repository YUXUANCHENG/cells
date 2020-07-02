// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"
#include "random"
#include <cmath>

using namespace std;

// split the packing system into smaller subsystems
void cellPacking2D::split_into_subspace() {
	int box;
	// create N[0] * N[1] subsystems
	if (subsystem == nullptr)
		subsystem = new subspace[N_systems[0] * N_systems[1]];

	// initialize subsystems
	for (int i = 0; i < N_systems[0] * N_systems[1]; i++) {
		(subsystem[i]).initialize(this, L, N_systems, i, dt0);
	}

	// assign cells into subsystems
	for (int ci = 0; ci < NCELLS; ci++) {
		box = look_for_new_box(cell(ci));
		migrate_into(box, &(cell(ci)));
		cell(ci).set_id(ci);
	}
};

// cashe list send to subsystems
void cellPacking2D::cashe_into(int i, vector<deformableParticles2D*>& cash_list) {
	subsystem[i].cashe_in(cash_list);
};

// migrate cells into subsystems
void cellPacking2D::migrate_into(int i, deformableParticles2D* const& migration) {
	subsystem[i].migrate_in(migration);
};

// figure out which box the cells belong to
int cellPacking2D::look_for_new_box(deformableParticles2D& cell) {
	int box_id = 0;
	int x_id = 0;
	int y_id = 0;

	// figure out x and y index
	x_id = floor(cell.cpos(0) / (L.at(0) / N_systems[0]));
	y_id = floor(cell.cpos(1) / (L.at(1) / N_systems[1]));

	//add periodic boundary just in case
	x_id = x_id % N_systems[0];
	y_id = y_id % N_systems[1];

	// convert into box id
	box_id = y_id * N_systems[0] + x_id;

	return box_id;
}

// initialization
void cellPacking2D::initialize_subsystems(int N_x, int N_y) {

	// set how many boxes along each direction
	if (N_systems.size() < 2)
	{
		N_systems.push_back(N_x);
		N_systems.push_back(N_y);
	}
	else
	{
		N_systems.at(0) = N_x;
		N_systems.at(1) = N_y;
	}
	// split
	split_into_subspace();

}


// reset system
void  cellPacking2D::reset_subsystems() {
	for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
		(subsystem[i]).reset();
}

void  cellPacking2D::delete_subsystems() {
	delete[] subsystem;
	subsystem = nullptr;
}
// active brownian simulation
void cellPacking2D::parallel_activityCOM_brownian(double T, double v0, double Dr, double vtau, double t_scale, int frames) {

	omp_set_num_threads(N_systems.at(0) * N_systems.at(1));
#pragma omp parallel
	{
		int i = omp_get_thread_num();
		(subsystem[i]).activityCOM_brownian_insub(T, v0, Dr, vtau, t_scale, frames);
	}
}


double cellPacking2D::max_length() {
	double max_length = 0;
	double length = 0;

	for (int ci = 0; ci < NCELLS; ci++) {
		length = cell(ci).max_length();
		if (length > max_length)
			max_length = length;
	}

	return max_length;
}

// compress isotropically to fixed packing fraction
void cellPacking2D::parallel_qsIsoCompression(double phiTarget, double deltaPhi, double Ftol) {
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

		parallel_fireMinimizeF(Ftol, Fcheck, Kcheck);
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
		parallel_fireMinimizeF(Ftol, Fcheck, Kcheck);

	}

}


// FIRE 2.0 force minimization with backstepping
void cellPacking2D::parallel_fireMinimizeF(double Ftol, double& Fcheck, double& Kcheck) {
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

void cellPacking2D::parallel_findJamming(double dphi0, double Ftol, double Ptol) {
	// local variables
	double Ptest, Ktest, Ftest;
	int NSTEPS, k, kmax, kr, nc, nr, ci, cj;
	cellPacking2D savedState;
	int NDOF;

	// get total number of degrees of freedom
	NDOF = 0;
	for (ci = 0; ci < NCELLS; ci++)
		NDOF += NDIM * cell(ci).getNV();

	// get initial packing fraction
	phi = packingFraction();

	// iterator
	k = 0;
	kmax = 1e5;

	// jamming variables
	bool jammed, overcompressed, undercompressed;
	double rH, rH0, rL, dr0, scaleFactor;

	// compute first dr0 based on current phi (i.e. non root search)
	dr0 = sqrt((phi + dphi0) / phi);

	// save initial state
	// r0 = sqrt(cell(0).geta0());
	saveState(savedState);

	// initialize as unjammed
	jammed = false;

	// phiJ bounds
	rH0 = -1;
	rH = -1;
	rL = -1;

	// initialize velocities
	double Tinit = 1e-6;
	//initializeVelocities(Tinit);

	double P = 0.0, vstarnrm = 0.0, fstarnrm = 0.0;
	bool converged = false;

	// loop until phi is the correct value
	while (!jammed && k < kmax) {
		// update iterator
		k++;

		// relax shapes (energies/forces calculated during FIRE minimization)
		omp_set_num_threads(N_systems.at(0) * N_systems.at(1));
#pragma omp parallel
		{
			int i = omp_get_thread_num();
			(subsystem[i]).fireMinimizeF_insub(Ftol, Ftest, Ktest, P, vstarnrm, fstarnrm, converged);
		}

		// calculate Ptest for comparison
		Ptest = 0.5 * (sigmaXX + sigmaYY) / (NDOF * L.at(0) * L.at(1));

		// remove rattlers
		//kr = 0;
		//nr = removeRattlers(kr);

		// update number of contacts
		//nc = totalNumberOfContacts();
		nc = Ncc;

		// boolean checks
		undercompressed = ((Ptest < 2.0 * Ptol && rH < 0) || (Ptest < Ptol&& rH > 0));
		overcompressed = (Ptest > 2.0 * Ptol && nc > 0);
		jammed = (Ptest < 2.0 * Ptol && Ptest > Ptol && nc > 0 && rH > 0);

		// output to console
		cout << "===================================================" << endl << endl << endl;
		cout << " 	quasistatic isotropic compression to jamming " << endl << endl;
		cout << "===================================================" << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* dphi 			= " << dphi0 << endl;
		cout << "	* phi 			= " << phi << endl;
		cout << "	* rH0 			= " << rH0 << endl;
		cout << "	* rH 			= " << rH << endl;
		cout << "	* rL 			= " << rL << endl;
		cout << "	* Ftest 		= " << Ftest << endl;
		cout << "	* Ktest 		= " << Ktest << endl;
		cout << "	* Ptest 		= " << Ptest << endl;
		cout << "	* # of contacts = " << nc << endl;
		//cout << "	* # of rattlers = " << nr << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl << endl;

		/*
		cout << "	* contact matrix:" << endl;

		// print contact matrix to console
		for (ci=0; ci<NCELLS; ci++){
			for (cj=0; cj<NCELLS; cj++){
				if (cj == ci)
					cout << "0" << "  ";
				else
					cout << contacts(ci,cj) << "  ";
			}
			cout << endl;
		}

		// print final two lines
		cout << endl << endl;
		*/

		// update packing fraction based on jamming check
		if (rL < 0) {
			// if still undercompressed, then grow until overjammed found
			if (undercompressed) {
				// set scale to normal compression
				scaleFactor = dr0;
			}
			// if first overcompressed, return to pre-overcompression state, to midpoint between phi and phiH
			else if (overcompressed) {
				// current = upper bound length scale r
				rH = sqrt(cell(0).geta0());

				// save this length scale
				rH0 = rH;

				// old = old length scale
				rL = rH / scaleFactor;

				// save overcompressed state
				saveState(savedState);

				// compute new scale factor
				scaleFactor = 0.5 * (rH + rL) / rH0;

				// print to console
				cout << "	-- -- overcompressed for first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else {
			// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
			if (undercompressed) {
				// current = new lower bound length scale r
				rL = sqrt(cell(0).geta0());

				// load state
				loadState(savedState);

				// compute new scale factor
				scaleFactor = 0.5 * (rH + rL) / rH0;

				// print to console
				cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;

			}
			else if (overcompressed) {
				// current = new upper bound length scale r
				rH = sqrt(cell(0).geta0());

				// load state
				loadState(savedState);

				// compute new scale factor
				scaleFactor = 0.5 * (rH + rL) / rH0;

				// print to console
				cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
			}
			else if (jammed) {
				cout << "	** At k = 0, jamming found!" << endl;
				cout << "	** phiJ = " << phi << endl;
				cout << "	** F = " << Ftest << endl;
				cout << "	** P = " << Ptest << endl;
				cout << "	** K = " << Ktest << endl;
				cout << "	** nc = " << nc << endl;
				cout << " WRITING JAMMED CONFIG TO .jam FILE" << endl;
				cout << " ENDING COMPRESSION SIMULATION" << endl;
				printJammedConfig_yc();
				phiPrintObject << phi << endl;
				printCalA();
				printContact();
				break;
			}
		}

		//saveState(savedState);

		// grow or shrink particles by scale factor
		scaleLengths(scaleFactor);

		// update new phi (only update here, do NOT calculate relaxed phi value)
		phi = packingFraction();

		if (k % 1 == 0) {
			printJammedConfig_yc();
			printCalA();
			printContact();
		}
	}

	if (k == kmax) {
		cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
		exit(1);
	}
}


// FIRE 2.0 force minimization with backstepping
void subspace::fireMinimizeF_insub(double Ftol, double& Fcheck, double& Kcheck, double& P, double& vstarnrm, double& fstarnrm, bool& converged) {
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
	double K, F;
	double xold, xnew, vold, vnew, fold;

	// reset time step
	double dt = dt0;

	int NCELLS = pointer_to_system->getNCELLS();
	int NVTOTAL = pointer_to_system->getNVTOTAL();

	// system stress
	double& sigmaXX_t = pointer_to_system->getSigmaXX();
	double& sigmaYY_t = pointer_to_system->getSigmaYY();
	double Pcheck = 0.5 * (sigmaXX_t + sigmaYY_t) / (NCELLS * L.at(0) * L.at(1));

	// system contact number
	double& Ncc_t = pointer_to_system->getNcc();
	double& Nvv_t = pointer_to_system->getNvv();

	// max length scale in the system
	double length = pointer_to_system->max_length();

	// calculate cashed fraction
	for (d = 0; d < NDIM; d++) {
		double spacing = L.at(d) / N_systems[d];
		//cashed_fraction.at(d) = pointer_to_system->scale_v(cashed_length) / spacing;
		cashed_fraction.at(d) = cashed_length * length / spacing;
		if (cashed_fraction.at(d) > 0.99) {
			cout << " Too much boxes for two little cells " << endl;
			cashed_fraction.at(d) = 0.99;
		}
	}

	calculateForces_insub();
#pragma omp barrier
#pragma omp master
	{
		for (d = 0; d < NDIM; d++)
			cout << "cashed fraction at " << d << " = " << cashed_fraction.at(d) << endl;
		const double Trescale = 1e-8 * NCELLS;
		//pointer_to_system->rescaleVelocities(Trescale);
		Fcheck = pointer_to_system->forceRMS();
		Kcheck = pointer_to_system->totalKineticEnergy();
	}
#pragma omp barrier

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
			sigmaXX_t = 0.0;
			sigmaYY_t = 0.0;
			Ncc_t = 0.0;
			Nvv_t = 0.0;
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
			if (k % pointer_to_system->getNPRINT() == 1)
				print_information();
		}
		// only master thread
		// get norms
#pragma omp barrier
#pragma omp master
		{
			vstarnrm = sqrt(vstarnrm);
			fstarnrm = sqrt(fstarnrm);

			// output some information to console
			if (k % pointer_to_system->getNPRINT() == 1) {
				cout << "===================================================" << endl << endl;
				cout << " 	FIRE MINIMIZATION, k = " << k << endl << endl;
				cout << "===================================================" << endl;
				cout << "	* Run data:" << endl;
				cout << "	* Kcheck 	= " << Kcheck << endl;
				cout << "	* Fcheck 	= " << Fcheck << endl;
				cout << "	* Pcheck 	= " << Pcheck << endl;
				cout << "	* phi 		= " << pointer_to_system->getphi() << endl;
				cout << "	* dt 		= " << dt << endl;
				cout << "	* alpha 	= " << alpha << endl;
				cout << "	* P 		= " << P << endl;
				cout << "	* Pdir 		= " << P / (vstarnrm * fstarnrm) << endl;
				cout << endl << endl;
			}

			Fcheck = 0.0;
			Kcheck = 0.0;
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

#pragma omp barrier
		if (k % update_freqency == 0) {
// perform cashing and migration
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
			// sent cashe to neighbors in x direction
#pragma omp barrier
#pragma omp critical
			{
				cashe_out(0);
			}
			// sent cashe to neighbors in y direction
#pragma omp barrier
#pragma omp critical
			{
				cashe_out(1);
			}
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


		// scale P and K for convergence checking
		// be careful about synchronization
#pragma omp critical
		{
			Fcheck += F;
			Kcheck += K / NCELLS;
			sigmaXX_t += sigmaXX;
			sigmaYY_t += sigmaYY;
			Ncc_t += Ncc;
			Nvv_t += Nvv;
		}
#pragma omp barrier
		// calculate force RMS
#pragma omp master
		{
			Fcheck = sqrt(Fcheck) / (NDIM * NVTOTAL);
			//Fcheck = pointer_to_system->forceRMS();
		}
#pragma omp barrier

		Pcheck = 0.5 * (sigmaXX_t + sigmaYY_t) / (NCELLS * L.at(0) * L.at(1));

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
			cout << "	** Pcheck = " << Pcheck << endl;
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
				cout << "	** Pcheck = " << Pcheck << endl;
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

void subspace::print_information() {

	cout << "This is box " << box_id << endl;
	cout << "resident: " << endl;
	for (int ci = 0; ci < resident_cells.size(); ci++)
		cout << resident_cells[ci]->get_id() << " ";
	cout << endl;
	cout << "cashed: " << endl;
	for (int ci = 0; ci < cashed_cells.size(); ci++)
		cout << cashed_cells[ci]->get_id() << " ";
	cout << endl;

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
	//frms = sqrt(frms) / (NDIM * NVTOTAL);

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
	cashed_cells.insert(cashed_cells.end(), cash_list.begin(), cash_list.end());
	//for (int i = 0; i < cash_list.size(); i++)
	//	cashed_cells.push_back(cash_list.at(i));
};

// add migrated cells
void subspace::migrate_in(deformableParticles2D* const& migration) {
	resident_cells.push_back(migration);
};

// reset cashe system
void subspace::reset_cashe() {
	if (!cashed_cells.empty()) {
		cashed_cells.clear();
	}
}

// reset the whole system
void subspace::reset() {

	reset_cashe();
	if (!resident_cells.empty()) {
		resident_cells.clear();
	}
}





// send cashe list 
void subspace::cashe_out(int direction) {

	// find neighbor boxes
	int lower_index = neighbor_box(direction, -1);
	int upper_index = neighbor_box(direction, +1);

	double spacing = L.at(direction) / N_systems[direction];
	// find boundaries of the current box
	double lower_boundary = find_boundary(direction, -1);
	double upper_boundary = find_boundary(direction, +1);

	// reset list
	cash_out_list_lower.clear();
	cash_out_list_upper.clear();

	// check if resident cells are near boundary
	if (!resident_cells.empty()) {
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			if (resident_cells[ci]->cpos(direction) < lower_boundary + cashed_fraction.at(direction) * spacing)

				cash_out_list_lower.push_back(resident_cells[ci]);
			if (resident_cells[ci]->cpos(direction) > upper_boundary - cashed_fraction.at(direction) * spacing)

				cash_out_list_upper.push_back(resident_cells[ci]);
		}
	}

	// check if cashed cells are near boundary, but only for y direction
	if (!cashed_cells.empty() && direction == 1) {
		for (int ci = 0; ci < cashed_cells.size(); ci++) {
			if (cashed_cells[ci]->cpos(direction) < lower_boundary + cashed_fraction.at(direction) * spacing &&
				cashed_cells[ci]->cpos(direction) >= lower_boundary)

				cash_out_list_lower.push_back(cashed_cells[ci]);
			if (cashed_cells[ci]->cpos(direction) > upper_boundary - cashed_fraction.at(direction) * spacing &&
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
	else if (upper_lower == 1)
		boundary = (current[direction] + 1) * spacing;

	return boundary;
}



// migrate cells to other boxes and remove from current box
void subspace::migrate_out() {

	int new_box_index;
	int cell_index;
	deformableParticles2D* target_cell;

	// stack indicates near boundary cells that migrate to neighbor boxes
	stack<int> migrate_out_list;
	stack<int> migrate_out_destination;

	/*
	// empty the stack
	while (!migrate_out_list.empty()) {
		migrate_out_list.pop();
	}

	while (!migrate_out_destination.empty()) {
		migrate_out_destination.pop();
	}
	*/

	// find the migration list
	for (int ci = 0; ci < resident_cells.size(); ci++) {
		new_box_index = pointer_to_system->look_for_new_box(*resident_cells.at(ci));
		if (new_box_index != box_id) {
			migrate_out_list.push(ci);
			migrate_out_destination.push(new_box_index);
		}
	}


	// migrate to other subsystems
	if (!migrate_out_list.empty()) {
		for (int i = 0; i < migrate_out_list.size(); i++) {
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



	// reset contacts before force calculation
	//resetContacts();
	Ncc = 0;
	Nvv = 0;



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
				if (inContact > 0) {
					// add to cell-cell contacts
					//addContact(ci, cj);
					Ncc++;

					// increment vertex-vertex contacts
					Nvv += inContact;
				}

			}

			if (!cashed_cells.empty()) {
				// forces between resident cell and cashed cell
				for (ck = 0; ck < cashed_cells.size(); ck++) {
					// notice that stress between resident and cashed cells are double counted
					inContact = resident_cells[ci]->vertexForce_cashed(*cashed_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
					if (inContact > 0) {
						// add to cell-cell contacts
						//addContact(ci, cj);
						// notice this is double counted
						Ncc += 1 / 2;

						// increment vertex-vertex contacts
						Nvv += inContact / 2;
					}
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

	// max length scale in the system
	double length = pointer_to_system->max_length();

	// calculate cashed fraction
	for (d = 0; d < NDIM; d++) {
		double spacing = L.at(d) / N_systems[d];
		//cashed_fraction.at(d) = pointer_to_system->scale_v(cashed_length) / spacing;
		cashed_fraction.at(d) = cashed_length * length / spacing;
		if (cashed_fraction.at(d) > 0.99) {
			cout << " Too much boxes for two little cells " << endl;
			cashed_fraction.at(d) = 0.99;
		}
	}

#pragma omp master
	{
		for (d = 0; d < NDIM; d++)
			cout << "cashed fraction at " << d << " = " << cashed_fraction.at(d) << endl;
	}

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
		if (count % update_freqency == 0) {
			// To avoid deadlock, migrate sequencially
#pragma omp critical
			{
				migrate_out();
			}

#pragma omp barrier
			// reset cashe
#pragma omp critical
			{
				reset_cashe();
			}
#pragma omp barrier
			// sent cashe to neighbors in x direction
#pragma omp critical
			{
				cashe_out(0);
			}
#pragma omp barrier
			// sent cashe to neighbors in y direction
#pragma omp critical
			{
				cashe_out(1);
			}
		}
#pragma omp barrier

		// use master thread to print
#pragma omp master
		{
			//phi = pointer_to_system->packingFraction();

			if (count % print_frequency == 0) {
				pointer_to_system->printJammedConfig_yc();
				//pointer_to_system->phiPrintObject << phi << endl;
				pointer_to_system->printCalA();
				//pointer_to_system->printContact();
				pointer_to_system->printV();
				cout << "t = " << t << endl;
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
		// Might cause some problem if only conserve momentum in subsystems
		//conserve_momentum();
#pragma omp barrier
#pragma omp master
		{
			pointer_to_system->conserve_momentum();
		}
#pragma omp barrier
	}



};

// conserve COM momentum
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

double subspace::max_length() {
	double max_length = 0;
	double length = 0;

	for (int ci = 0; ci < resident_cells.size(); ci++) {
		length = resident_cells[ci]->max_length();
		if (length > max_length)
			max_length = length;
	}

	return max_length;
}


