/*

	Example file to generate a gel of cells
	from an initial bidisperse sphere packing,
	and to decompress at a fixed rate
	rather than QS

*/

// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>

// use std name space
using namespace std;

// define PI
const double PI = 4.0*atan(1);

// length paramaters
const int NT 					= 5e4;
const int NPRINT				= 50;

// simulation constants
const double sizeDispersion 	= 0.1;			// size dispersion (std dev of cell sizes)
const double timeStepMag 		= 0.005;			// time step in MD units (zeta * lenscale / forcescale)

// disk constants
const double phiDisk	 		= 0.3;			// initial packing fraction of disks (sets boundary)

// compression constants
const double phiTarget			= 0.95;			// cell packing fraction (regardless of final pressure)
const double deltaPhi			= 0.001;		// compression step size

// gelation constants
const double phiGel 			= 0.3;			// final packing fraction
const double gelRate 			= 1e-4;			// rate of size decrease (i.e. area loss relative to initial box area)
const double varPerimRate 		= 0.01;			// rate of relaxation to deformed perimeter
const double aGelation			= 0.05;			// attraction parameter during gelation sim

// force parameters
const double kl 			= 1.0;				// perimeter force constant
const double ka 			= 10.0;				// area force constant
const double gam 			= 0.0;				// surface tension force constant
const double kb 			= 2.0;				// bending energy constant
const double kint 			= 1.0;				// interaction energy constant
const double del 			= 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
const double aInitial 		= 0.0;				// attraction parameter to start
const double da 			= 0.001;			// attraction increment

// deformability
const double calA0 			= 1.2;				// ratio of preferred perimeter^2 to preferred area

// main function
int main()
{
	// local variables

	// output files
	string posFile = "pos.test";
	string enFile = "en.test";

	// system details
	int NCELLS 		= 4;
	int NV			= 50;
	int seed 		= 1;
	double Ltmp 	= 1.0;
	double plThresh = 1e-2;
	double dl0		= 1e-3;
	double calA0max = 1.8;

	// instantiate object
	cout << "	** Instantiating object for initial disk packing to be turned into a cell packing" << endl;
	cout << "	** NCELLS = " << NCELLS << endl;
	cellPacking2D packingObject(NCELLS,NT,NPRINT,Ltmp,seed); 	// NOTE: NEED TO MAKE NEW CONSTRUCTOR, EVERYTHING ELSE DONE IN initializeGel AND regularPolygon FUNCTIONS

	// set initial conditions as if disks in box with given packing fraction (sets boundary size)
	cout << "	** Initializing gel at phiDisk = " << phiDisk << " using SP model" << endl;
	packingObject.initializeGel(NV, phiDisk, sizeDispersion, del, ka);

	// set deformability, force values
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,aInitial);

	// update time scale
	packingObject.setdt(timeStepMag);

	// open position output file
	packingObject.openPackingObject(posFile);
	packingObject.openEnergyObject(enFile);

	// set preferred angle
	// double theta0 = 2.0*PI/NV;
	// packingObject.cell(0).setc0Angle(theta0);

	// run test NVE
	packingObject.initializeVelocities(0.1);
	packingObject.cellNVE();

	/*
	// compress to set packing fraction using FIRE, pressure relaxation
	cout << "	** QS compresison protocol to phiTarget = " << phiTarget << endl;
	packingObject.qsIsoCompression(phiTarget,deltaPhi);

	// set new attraction
	packingObject.gelForceVals(calA0,kl,ka,gam,kb,kint,del,aGelation);

	// -- decrease phi as if boundary was growing: phi(t) = phi(0)/(1 + a*t)
	cout << "	** Running qs gel extension simulation with deltaPhi = " << deltaPhi << ", calA0max = " << calA0max << endl;
	packingObject.qsIsoGelRatchet(phiGel,deltaPhi,plThresh,dl0,calA0max,timeStepMag);
	*/

	return 0;
}
























