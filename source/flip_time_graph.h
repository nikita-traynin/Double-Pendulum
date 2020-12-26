#include "helper_functions.h"

void flip_time_graph(char *argv[])
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
	
	//set the maximum energy error we tolerate 
    const double energy_error = .01; //in Newton-meters 

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
    double energy_estimate, energy_diff;

	//allocate space for storing final flipping times 
	double* time_taken = (double*)malloc(sizeof(double)*m*n);
	
	//initialize grid axes
    for(int i = 0; i < m; i++)
        as[i] = min_angle_a + i*angle_step_a;
    
    for(int i = 0; i < n; i++)
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
			
			//check if there is enough energy for a flip; if not, the skip to the next iteration
			if (3*cos(a) + cos(c) > 2) {
				*(time_taken + j*n + i) = max_time + h;
				continue;
			}
			
			//begin the simulation 
            while (1)
            {
				//calculate the next time and spatial values 
                calculate_next_values(&a, &b, &c, &d, &t, h);
				
				//calculate the energy of the system at the moment 
                energy_estimate = calculate_energy(a,b,c,d); 
                energy_diff = fabs((energy_estimate - initial_energy)/initial_energy);
				
				//check if our error is tolerable 
                if(energy_diff < energy_error)
                {
					//check if there's been a flip 
                    if( a > (as[i] + 2.0*pi) ||  a < (as[i] - 2.0*pi) || c > (cs[j] + 2.0*pi) || c < (cs[j] - 2.0*pi) )
                        {*(time_taken + j*n + i) = t;  break;}

					//check if we've surpassed the maximum simulation time 
                    else if (t > max_time)
                        {*(time_taken + j*n + i) = t;  break;}
                }
				
				//error too big => store special code 
                else
                {
                    *(time_taken + j*n + i) = -1.0; //out of range for color scheme
                    break;
                }

            }
        }
    }
	printf("DONE SIMULATING");
    draw_graph(time_taken, max_time, h, m, n);

    return;
}

