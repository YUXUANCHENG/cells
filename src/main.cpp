#include "cell_jamming.h"
#include "CLI.h"



// use std name space
using namespace std;



int main(int argc, char const *argv[]){

	jamming main_function;

	//main_function.unjam();
	//main_function.soft_particle_limit();
	//main_function.confluency();
	//main_function.confluency(argv);
	//main_function.unjam_N();
	//main_function.active_brownian();
	//main_function.test();
	//main_function.test_ground_shape();
	//main_function.const_ground_shape();
	//main_function.const_ground_shape_arg(argv);
	//main_function.cell_NVE_arg(argv);
	//main_function.sp_NVE_arg(argv);
	//main_function.random_NVE_arg(argv);
	//main_function.Hopper(argv);
	//main_function.Hopper_width(argv);
	//main_function.sp_NVE_probe_arg(argv);
	//main_function.cell_NVE_probe_arg(argv);
	//main_function.bumpy_NVE_arg(argv);
	//main_function.mesaure_ground_shape();
	//main_function.gravity();

	//Bumpy_CLI cli;
	DPM_CLI cli;
	//cli.findJamming(argv);
	cli.NVE(argv);

	system("pause");
	return 0;

}


