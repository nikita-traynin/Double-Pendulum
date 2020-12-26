#include "helper_functions.h"

void error_graph(char *argv[])
{
	/*
	argv[1] (float) and argv[2] (integer) are the temporal step size of our approximation in scientific notation. 
	argv[3] (integer) is the max_time of the simulation
	argv[4] (integer) is the resolution of the graph (Anything more than 75 takes a while)
	*/
	
	//sets the step size for the temporal discretization
    const float h_coeff = atof(argv[1]);
    const int h_exp = atoi(argv[2]);
    const double h = h_coeff * pow(10.0, (double)(h_exp));

    //sets the max time the simulation can be run 
    const double max_time = atoi(argv[3]);

	//set initial angular velocities (usually 0, but experiment!)
    const double theta_v0 = 0.0;
    const double phi_v0 = 0.0;

	//no maximum error - we store all the error values

	
    //set the domain, usually (0, pi) X (0, pi)
	const double min_angle_a = 0;
    const double max_angle_a = pi;
	const double min_angle_c = 0;
    const double max_angle_c = pi;
	
	//grid resolution and grid intervals. (usually a square graph with m=n, but can experiment with rectangular)
    const int m = atoi(argv[4]);
    const int n = m;
	const double angle_step_a = (max_angle_a - min_angle_a) / m;
    const double angle_step_c = (max_angle_c - min_angle_c) / n;

	//declare the axes of the grid
    double as[m];
    double cs[n];

	//the four values to solve for: first angle, first angular velocity, second angle, second angular velocity
	//also time variable 
    double a,b,c,d; 
	double t;
	
	//stats of the double pendulum
    double energy_estimate;
	
	//allocate space to store the energy relative errors
	double* energy_relative_errors = (double*)malloc(sizeof(double)*m*n);

	//initialize grid axes
    for (int i = 0; i < m; i++)
        as[i] = min_angle_a + i*angle_step_a;

    for (int i = 0; i < n; i++)
        cs[i] = min_angle_c + i*angle_step_c;
	
	
	//loop through all initial angle values
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++) 
        {
			//Set the initial positions and angular velocities 
            t = 0.0;
            a = as[i];
            b = theta_v0;
            c = cs[j];
            d = phi_v0;
			
			//calculate the energy we start with 
            double initial_energy = calculate_energy(a,b,c,d);

			//begin the simulation
            while (1 == 1)
            {
                calculate_next_values(&a, &b, &c, &d, &t, h);
				
				//end the simulation
                if (t > max_time)
                {
                    energy_estimate = calculate_energy(a, b, c, d);
                    *(energy_relative_errors + j*n + i) = fabs((energy_estimate - initial_energy)/initial_energy);
                    break;
                }

            }
        }
    }
	draw_graph(energy_relative_errors, .05, h, m, n);
}