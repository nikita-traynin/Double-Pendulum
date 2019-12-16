#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define pi 3.14159265359


void graph_first_inversion_time(char *argv[])
{
	//first arg is size, second is max_time
	//third and 4th args are for the maximum error
	//third arg is the value, and fourth arg is the coefficient of 
		//the power of ten
	
    const int RED = 0;
    const int GREEN = 1;
    const int BLUE = 2;
    const double SEGS = 4.0;
    const double FIR = .2;
    const double SEC = .3;
    const double THIR = .5;
    const double FOUR = 3.0;


    int i,j;
	
	//changes step size 
    const double t_0 = 0.0;
    const int h_coeff = atoi(argv[3]);
    const int h_exp = atoi(argv[4]);
    const double h = h_coeff * pow(10.0, (double)(h_exp));

	//changes maximum time the motion is simulated for
    const double max_time = atoi(argv[2]);
    const double segment = max_time/SEGS;
    const double first_segment = (FIR)*segment;
    const double second_segment = (FIR + SEC)*segment;
    const double third_segment = (FIR + SEC + THIR)*segment;
    const double fourth_segment = (FIR + SEC + THIR + FOUR)*segment + h;
    const double theta_v0 = 0.0;
    const double phi_v0 = 0.0;

	//physical values - mass, length, and g constant
    const double m1 = 1.0;
    const double m2 = 1.0;
    const double L1 = 1.0;
    const double L2 = 1.0;
    const double g = 1.0;

    const double energy_error = .01;

    //zooming on the section of the domain we are interested in 
    const double a_left_1 = 0;
    const double a_right_1 = 1;
    const double c_bottom_1 = 0;
    const double c_top_1 = 1;

    const double max_angle_a = (a_right_1)*pi;
    const double min_angle_a = (a_left_1)*pi;
    const double max_angle_c = (c_top_1)*pi;
    const double min_angle_c = (c_bottom_1)*pi;

    const double angle_difference_a = max_angle_a - min_angle_a;
    const double angle_difference_c = max_angle_c - min_angle_c;

    printf("%f, %f", angle_difference_a, angle_difference_c);

    const int m = atoi(argv[1]);
    const int n = m;

	printf("\n%d, %d, %d, %d", atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));

    const double angle_step_a = angle_difference_a / m;
    const double angle_step_c = angle_difference_c / n;

    double as[m];
    double cs[n];


    int* output;
    output = (int*)malloc(sizeof(int)*m*n*3);

    double a,b,c,d; //the four values to solve for: first angle, first angular velocity, second angle, second angular velocity
    double k1a, k2a, k3a, k4a, k1b, k2b, k3b, k4b, k1c, k2c, k3c, k4c,k1d,k2d, k3d, k4d;
    double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, avg_a, avg_b, avg_c, avg_d; //vars needed for runge kutta
    double time_taken, t, energy_est, energy_diff;


    for(i = 0; i < m; i++)
    {
        as[i] = min_angle_a + i*angle_step_a;
    }

    for(i = 0; i < n; i++)
    {
        cs[i] = min_angle_c + i*angle_step_c;

    }


    for (i = 0; i < m; i++) //for angles min_angle to max_angle (theta)
    {
        for (j = 0; j < n; j++) //for angles min_angle to max_angle (phi)
        {
            t = t_0;
            a = as[i];
            b = theta_v0;
            c = cs[j];
            d = phi_v0;
            double initial_energy = -(m1+m2)*g*(L1*cos(a)) - m2*g*L2*cos(c);

            while (1 == 1)
            {

                k1a = b;
                k1b = (  ( m2*L2*d*d*sin(a-c) ) + ( m1+m2 )*(g*sin(a)) + (m2*L1*b*b)*sin(a-c)*cos(a-c) - m2*g*sin(c)*cos(a-c)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a-c)),2) ));
                k1c = d;
                k1d = ( (cos(a-c)) * (  ((m1+m2)*g*sin(a)) + (m2*L2*d*d*sin(a-c))  )  +  (m1+m2)*( -g*sin(c) + L1*b*b*sin(a-c)) )  /  (L2*( m1 + m2 - m2*pow((cos(a-c)),2) ));

                a1 = a + k1a*(h/2);
                b1 = b + k1b*(h/2);
                c1 = c + k1c*(h/2);
                d1 = d + k1d*(h/2);

                k2a = b1;
                k2b = (  ( m2*L2*d1*d1*sin(a1-c1) ) + ( m1+m2 )*(g*sin(a1)) + (m2*L1*b1*b1)*sin(a1-c1)*cos(a1-c1) - m2*g*sin(c1)*cos(a1-c1)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a1-c1)),2) ));
                k2c = d1;
                k2d = ( (cos(a1-c1)) * (  ((m1+m2)*g*sin(a1)) + (m2*L2*d1*d1*sin(a1-c1))  )  +  (m1+m2)*( -g*sin(c1) + L1*b1*b1*sin(a1-c1)) )  /  (L2*( m1 + m2 - m2*pow((cos(a1-c1)),2) ));

                a2 = a + k2a*(h/2);
                b2 = b + k2b*(h/2);
                c2 = c + k2c*(h/2);
                d2 = d + k2d*(h/2);

                k3a = b2;
                k3b = (  ( m2*L2*d2*d2*sin(a2-c2) ) + ( m1+m2 )*(g*sin(a2)) + (m2*L1*b2*b2)*sin(a2-c2)*cos(a2-c2) - m2*g*sin(c2)*cos(a2-c2)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a2-c2)),2) ));
                k3c = d2;
                k3d = ( (cos(a2-c2)) * (  ((m1+m2)*g*sin(a2)) + (m2*L2*d2*d2*sin(a2-c2))  )  +  (m1+m2)*( -g*sin(c2) + L1*b2*b2*sin(a2-c2)) )  /  (L2*( m1 + m2 - m2*pow((cos(a2-c2)),2) ));

                a3 = a + k3a*h;
                b3 = b + k3b*h;
                c3 = c + k3c*h;
                d3 = d + k3d*h;

                k4a = b3;
                k4b = (  ( m2*L2*d3*d3*sin(a3-c3) ) + ( m1+m2 )*(g*sin(a3)) + (m2*L1*b3*b3)*sin(a3-c3)*cos(a3-c3) - m2*g*sin(c3)*cos(a3-c3)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a3-c3)),2) ));
                k4c = d3;
                k4d = ( (cos(a3-c3)) * (  ((m1+m2)*g*sin(a3)) + (m2*L2*d3*d3*sin(a3-c3))  )  +  (m1+m2)*( -g*sin(c3) + L1*b3*b3*sin(a3-c3)) )  /  (L2*( m1 + m2 - m2*pow((cos(a3-c3)),2) ));

                avg_a = (1.0/6.0)*(k1a + 2*k2a + 2*k3a + k4a);
                avg_b = (1.0/6.0)*(k1b + 2*k2b + 2*k3b + k4b);
                avg_c = (1.0/6.0)*(k1c + 2*k2c + 2*k3c + k4c);
                avg_d = (1.0/6.0)*(k1d + 2*k2d + 2*k3d + k4d);

                a += (avg_a)*h;
                b += (avg_b)*h;
                c += (avg_c)*h;
                d += (avg_d)*h;
                t += h;

                energy_est = .5*(m1+m2)*L1*L1*b*b + .5*m2*(L2*L2*d*d + 2*L1*L2*b*d*cos(a-c)) - (m1+m2)*g*(L1*cos(a)) - m2*g*L2*cos(c);
                energy_diff = fabs((energy_est - initial_energy)/initial_energy);

                if(energy_diff < energy_error)
                {
                    if( a > (as[i] + 2.0*pi) ||  a < (as[i] - 2.0*pi) || c > (cs[j] + 2.0*pi) || c < (cs[j] - 2.0*pi) )
                        {time_taken = t;  /*printf("%f, %f: %f for %f\n", as[i], cs[j], energy_diff, time_taken);*/ break;}

                    else if (t > max_time)
                        {time_taken = t;  break;}
                }
                else
                {
                    time_taken = fourth_segment + 1; //out of range for color scheme
                    //printf("%f, %f: %f\n", as[i], cs[j], energy_diff);
                    break;
                }

            }

            //CYANS-BLUES
            if(time_taken < first_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = (time_taken) * (255.0/first_segment);
                *(output + 3*j*n + i*3 + BLUE) = 255;
            }
            //GREENS-CYANS
            else if(time_taken >= first_segment && time_taken < second_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = (second_segment - time_taken) * (255.0/(second_segment-first_segment));
            }
            //YELLOWS-GREENS
            else if(time_taken >= second_segment && time_taken < third_segment)
            {
                *(output + 3*j*n + i*3 + RED) = (time_taken - second_segment) * (255.0/(third_segment-second_segment));
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //REDS-YELLOWS
            else if(time_taken >= third_segment && time_taken <= fourth_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 255;
                *(output + 3*j*n + i*3 + GREEN) = (fourth_segment - time_taken) * (255.0/(fourth_segment-third_segment));
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //BLACK: LARGE ERROR IN CALCULATION
            else
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = 0;
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
        }
    }

    FILE * fp;

    char* buffer;
    buffer = (char*)malloc(10*sizeof(char));

    char filename[50];
    strcpy(filename, "image");

    strcat(filename, itoa(m, buffer, 10));

    strcat(filename, "x");

    strcat(filename, itoa(n, buffer,10));

    strcat(filename, "_");

    strcat(filename, itoa(round(max_time), buffer, 10));

    strcat(filename, "_");

    strcat(filename, itoa(h_coeff, buffer, 10));

    strcat(filename, "e");

    strcat(filename, itoa(h_exp, buffer,10));

    strcat(filename, "_");

    fp = fopen(filename, "w");

    //PPM HEADER with width,height,and gray-scale (ASCII PIXELMAP)
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", m, n);
    fprintf(fp, "%d\n", 255);

    for(j = n-1; j >= 0; j--)
    {
        for(i = 0; i < m; i++)
        {
            fprintf(fp, "%d %d %d", *(output + 3*n*j + 3*i + RED), *(output + 3*n*j + 3*i + GREEN), *(output + 3*n*j + 3*i + BLUE));
            fprintf(fp, "   ");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    return;
}

void error_graph()
{
    const int RED = 0;
    const int GREEN = 1;
    const int BLUE = 2;
    const double SEGS = 4.0;
    const double FIR = .2;
    const double SEC = .3;
    const double THIR = .5;
    const double FOUR = 3.0;


    int i,j;

    const double t_0 = 0.0;
    const int h_coeff = 3;
    const int h_exp = -2;
    const double h = h_coeff * pow(10.0, (double)(h_exp));

    const double max_time = 300.0;

    const double max_error = .7;
    const double segment = max_error/SEGS;
    const double first_segment = (FIR)*segment;
    const double second_segment = (FIR + SEC)*segment;
    const double third_segment = (FIR + SEC + THIR)*segment;
    const double fourth_segment = (FIR + SEC + THIR + FOUR)*segment + h;

    const double theta_v0 = 0.0;
    const double phi_v0 = 0.0;

    const double m1 = 1.0;
    const double m2 = 1.0;
    const double L1 = 1.0;
    const double L2 = 1.0;
    const double g = 1.0;

    const double max_angle_a = pi;
    const double min_angle_a = 0;
    const double max_angle_c = pi;
    const double min_angle_c = 0;

    const double angle_difference_a = max_angle_a - min_angle_a;
    const double angle_difference_c = max_angle_c - min_angle_c;

    printf("%f, %f", angle_difference_a, angle_difference_c);

    const int m = 1000;
    const int n = 1000;

    const double angle_step_a = angle_difference_a / m;
    const double angle_step_c = angle_difference_c / n;

    double as[m];
    double cs[n];


    int* output;
    output = (int*)malloc(sizeof(int)*m*n*3);

    double a,b,c,d; //the four values to solve for: first angle, first angular velocity, second angle, second angular velocity
    double k1a, k2a, k3a, k4a, k1b, k2b, k3b, k4b, k1c, k2c, k3c, k4c,k1d,k2d, k3d, k4d;
    double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, avg_a, avg_b, avg_c, avg_d; //vars needed for runge kutta
    double t, energy_est, energy_diff;


    for(i = 0; i < m; i++)
    {
        as[i] = min_angle_a + i*angle_step_a;
    }

    for(i = 0; i < n; i++)
    {
        cs[i] = min_angle_c + i*angle_step_c;

    }


    for (i = 0; i < m; i++) //for angles min_angle to max_angle (theta)
    {
        for (j = 0; j < n; j++) //for angles min_angle to max_angle (phi)
        {
            t = t_0;
            a = as[i];
            b = theta_v0;
            c = cs[j];
            d = phi_v0;
            double initial_energy = -(m1+m2)*g*(L1*cos(a)) - m2*g*L2*cos(c);

            while (1 == 1)
            {

                k1a = b;
                k1b = (  ( m2*L2*d*d*sin(a-c) ) + ( m1+m2 )*(g*sin(a)) + (m2*L1*b*b)*sin(a-c)*cos(a-c) - m2*g*sin(c)*cos(a-c)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a-c)),2) ));
                k1c = d;
                k1d = ( (cos(a-c)) * (  ((m1+m2)*g*sin(a)) + (m2*L2*d*d*sin(a-c))  )  +  (m1+m2)*( -g*sin(c) + L1*b*b*sin(a-c)) )  /  (L2*( m1 + m2 - m2*pow((cos(a-c)),2) ));

                a1 = a + k1a*(h/2);
                b1 = b + k1b*(h/2);
                c1 = c + k1c*(h/2);
                d1 = d + k1d*(h/2);

                k2a = b1;
                k2b = (  ( m2*L2*d1*d1*sin(a1-c1) ) + ( m1+m2 )*(g*sin(a1)) + (m2*L1*b1*b1)*sin(a1-c1)*cos(a1-c1) - m2*g*sin(c1)*cos(a1-c1)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a1-c1)),2) ));
                k2c = d1;
                k2d = ( (cos(a1-c1)) * (  ((m1+m2)*g*sin(a1)) + (m2*L2*d1*d1*sin(a1-c1))  )  +  (m1+m2)*( -g*sin(c1) + L1*b1*b1*sin(a1-c1)) )  /  (L2*( m1 + m2 - m2*pow((cos(a1-c1)),2) ));

                a2 = a + k2a*(h/2);
                b2 = b + k2b*(h/2);
                c2 = c + k2c*(h/2);
                d2 = d + k2d*(h/2);

                k3a = b2;
                k3b = (  ( m2*L2*d2*d2*sin(a2-c2) ) + ( m1+m2 )*(g*sin(a2)) + (m2*L1*b2*b2)*sin(a2-c2)*cos(a2-c2) - m2*g*sin(c2)*cos(a2-c2)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a2-c2)),2) ));
                k3c = d2;
                k3d = ( (cos(a2-c2)) * (  ((m1+m2)*g*sin(a2)) + (m2*L2*d2*d2*sin(a2-c2))  )  +  (m1+m2)*( -g*sin(c2) + L1*b2*b2*sin(a2-c2)) )  /  (L2*( m1 + m2 - m2*pow((cos(a2-c2)),2) ));

                a3 = a + k3a*h;
                b3 = b + k3b*h;
                c3 = c + k3c*h;
                d3 = d + k3d*h;

                k4a = b3;
                k4b = (  ( m2*L2*d3*d3*sin(a3-c3) ) + ( m1+m2 )*(g*sin(a3)) + (m2*L1*b3*b3)*sin(a3-c3)*cos(a3-c3) - m2*g*sin(c3)*cos(a3-c3)  )  /  (L1*(-m1 - m2 + m2*pow((cos(a3-c3)),2) ));
                k4c = d3;
                k4d = ( (cos(a3-c3)) * (  ((m1+m2)*g*sin(a3)) + (m2*L2*d3*d3*sin(a3-c3))  )  +  (m1+m2)*( -g*sin(c3) + L1*b3*b3*sin(a3-c3)) )  /  (L2*( m1 + m2 - m2*pow((cos(a3-c3)),2) ));

                avg_a = (1.0/6.0)*(k1a + 2*k2a + 2*k3a + k4a);
                avg_b = (1.0/6.0)*(k1b + 2*k2b + 2*k3b + k4b);
                avg_c = (1.0/6.0)*(k1c + 2*k2c + 2*k3c + k4c);
                avg_d = (1.0/6.0)*(k1d + 2*k2d + 2*k3d + k4d);

                a += (avg_a)*h;
                b += (avg_b)*h;
                c += (avg_c)*h;
                d += (avg_d)*h;
                t += h;

                if (t > max_time)
                {
                    energy_est = .5*(m1+m2)*L1*L1*b*b + .5*m2*(L2*L2*d*d + 2*L1*L2*b*d*cos(a-c)) - (m1+m2)*g*(L1*cos(a)) - m2*g*L2*cos(c);
                    energy_diff = fabs((energy_est - initial_energy)/initial_energy);
                    break;
                }

            }

            //CYANS-BLUES (SMALL ERROR)
            if(energy_diff < first_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = (energy_diff) * (255.0/first_segment);
                *(output + 3*j*n + i*3 + BLUE) = 255;
            }
            //GREENS-CYANS
            else if(energy_diff >= first_segment && energy_diff < second_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = (second_segment - energy_diff) * (255.0/(second_segment-first_segment));
            }
            //YELLOWS-GREENS
            else if(energy_diff >= second_segment && energy_diff < third_segment)
            {
                *(output + 3*j*n + i*3 + RED) = (energy_diff - second_segment) * (255.0/(third_segment-second_segment));
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //REDS-YELLOWS (LARGE ERROR)
            else if(energy_diff >= third_segment && energy_diff <= fourth_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 255;
                *(output + 3*j*n + i*3 + GREEN) = (fourth_segment - energy_diff) * (255.0/(fourth_segment-third_segment));
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //BLACK: LARGE ERROR IN CALCULATION
            else
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = 0;
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
        }
    }

    FILE * fp;

    char* buffer;
    buffer = (char*)malloc(10*sizeof(char));

    char filename[50];
    strcpy(filename, "error_image");

    strcat(filename, itoa(m, buffer, 10));

    strcat(filename, "x");

    strcat(filename, itoa(n, buffer,10));

    strcat(filename, "_");

    strcat(filename, itoa(round(max_time), buffer, 10));

    strcat(filename, "_");

    strcat(filename, itoa(h_coeff, buffer, 10));

    strcat(filename, "e");

    strcat(filename, itoa(h_exp, buffer,10));

    strcat(filename, "_");

    fp = fopen(filename, "w");

    //PPM HEADER with width,height,and gray-scale (ASCII PIXELMAP)
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", m, n);
    fprintf(fp, "%d\n", 255);

    for(j = n-1; j >= 0; j--)
    {
        for(i = 0; i < m; i++)
        {
            fprintf(fp, "%d %d %d", *(output + 3*n*j + 3*i + RED), *(output + 3*n*j + 3*i + GREEN), *(output + 3*n*j + 3*i + BLUE));
            fprintf(fp, "   ");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}
