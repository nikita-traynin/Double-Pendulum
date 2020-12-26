#ifndef HELPER_H
#define HELPER_H

#include <stdlib.h> //for malloc 
#include <stdio.h> //for FILE, fprintf, printf, etc.
#include <math.h> //for cos and sin

#define pi 3.14159265359


//physical values - mass, length, and g constant (kilograms, meters, m^3*kg^-1*s^-2, respectively)
const double m1 = 1.0;
const double m2 = 1.0;
const double L1 = 1.0;
const double L2 = 1.0;
const double g = 1.0;

//stored helper methods 
void calculate_next_values(double* a, double* b, double* c, double* d, double* t, double h);
void draw_graph(double* time_taken, double max_time, double h, int m, int n);
double calculate_energy(double a, double b, double c, double d);
double first_equation(double a, double b, double c, double d);
double second_equation(double a, double b, double c, double d);

void calculate_next_values(double* a, double* b, double* c, double* d, double* t, double h)
{
	//the four stage runge kutta approximate slopes for each four variables
    double k1a, k2a, k3a, k4a, k1b, k2b, k3b, k4b, k1c, k2c, k3c, k4c,k1d,k2d, k3d, k4d;
	
	//the intermediate values 
    double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3; 
	
	k1a = *b;
	k1b = first_equation(*a, *b, *c, *d);
	k1c = *d;
	k1d = second_equation(*a, *b, *c, *d);

	a1 = *a + k1a*(h/2);
	b1 = *b + k1b*(h/2);
	c1 = *c + k1c*(h/2);
	d1 = *d + k1d*(h/2);

	k2a = b1;
	k2b = first_equation(a1, b1, c1, d1);
	k2c = d1;
	k2d = second_equation(a1, b1, c1, d1);

	a2 = *a + k2a*(h/2);
	b2 = *b + k2b*(h/2);
	c2 = *c + k2c*(h/2);
	d2 = *d + k2d*(h/2);

	k3a = b2;
	k3b = first_equation(a2, b2, c2, d2);
	k3c = d2;
	k3d = second_equation(a2, b2, c2, d2);

	a3 = *a + k3a*h;
	b3 = *b + k3b*h;
	c3 = *c + k3c*h;
	d3 = *d + k3d*h;

	k4a = b3;
	k4b = first_equation(a3, b3, c3, d3);
	k4c = d3;
	k4d = second_equation(a3, b3, c3, d3);
	
	*a += (1.0/6.0)*(k1a + 2*k2a + 2*k3a + k4a)*h;
	*b += (1.0/6.0)*(k1b + 2*k2b + 2*k3b + k4b)*h;
	*c += (1.0/6.0)*(k1c + 2*k2c + 2*k3c + k4c)*h;
	*d += (1.0/6.0)*(k1d + 2*k2d + 2*k3d + k4d)*h;
	*t += h;
}


void draw_graph(double* time_taken, double max_time, double h, int m, int n) 
{
	//variables for coloring the graph later (4 color segments) 
    enum colors{RED, GREEN, BLUE};
    const double SEGS = 4.0;
    const double FIR = .2;
    const double SEC = .3;
    const double THIR = .5;
    const double FOUR = 3.0;
	
	//the ranges of the color palette
    const double segment = max_time/SEGS;
    const double first_segment = (FIR)*segment;
    const double second_segment = (FIR + SEC)*segment;
    const double third_segment = (FIR + SEC + THIR)*segment;
	
	//allocate space for output
	int * output = (int*)malloc(sizeof(int)*3*m*n);
	
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			double time = *(time_taken + j*n + i);
			//CYANS-BLUES (VERY QUICKLY)
			if(time < first_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = (time) * (255.0/first_segment);
                *(output + 3*j*n + i*3 + BLUE) = 255;
            }
            //GREENS-CYANS (QUICKLY)
            else if(time >= first_segment && time < second_segment)
            {
                *(output + 3*j*n + i*3 + RED) = 0;
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = (second_segment - time) * (255.0/(second_segment-first_segment));
            }
            //YELLOWS-GREENS (SLOWLY)
            else if(time >= second_segment && time < third_segment)
            {
                *(output + 3*j*n + i*3 + RED) = (time - second_segment) * (255.0/(third_segment-second_segment));
                *(output + 3*j*n + i*3 + GREEN) = 255;
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //REDS-YELLOWS (VERY SLOWLY OR NEVER)
            else if(time >= third_segment && time <= max_time + h)
            {
                *(output + 3*j*n + i*3 + RED) = 255;
                *(output + 3*j*n + i*3 + GREEN) = (max_time - time) * (255.0/(max_time-third_segment));
                *(output + 3*j*n + i*3 + BLUE) = 0;
            }
            //BLACK: LARGE ERROR IN CALCULATION
            else if (time == -1.0)
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

    char* filename = "test";
    /* strcpy(filename, "image");
    strcat(filename, itoa(m, buffer, 10));
    strcat(filename, "x");
    strcat(filename, itoa(n, buffer,10));
    strcat(filename, "_");
    strcat(filename, itoa(round(max_time), buffer, 10));
    strcat(filename, "_");
    strcat(filename, itoa(h_coeff, buffer, 10));
    strcat(filename, "e");
    strcat(filename, itoa(h_exp, buffer,10));
    strcat(filename, "_"); */

    fp = fopen(filename, "w");
	
	if(!fp) {
		printf("File couldn't open!");
		return;
	}
    //PPM HEADER with width,height,and gray-scale (ASCII PIXELMAP)
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", m, n);
    fprintf(fp, "%d\n", 255);

    for(int j = n-1; j >= 0; j--)
    {
        for(int i = 0; i < m; i++)
        {
            fprintf(fp, "%d %d %d", *(output + 3*n*j + 3*i + RED), *(output + 3*n*j + 3*i + GREEN), *(output + 3*n*j + 3*i + BLUE));
            fprintf(fp, "   ");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}


double calculate_energy(double a, double b, double c, double d)
{
	return .5*(m1+m2)*L1*L1*b*b + .5*m2*(L2*L2*d*d + 2*L1*L2*b*d*cos(a-c)) - (m1+m2)*g*(L1*cos(a)) - m2*g*L2*cos(c);
}
double first_equation(double a, double b, double c, double d) 
{
	double sine_val = sin(a-c);
	double cos_val = cos(a-c);
	return (  ( m2*L2*d*d*sine_val ) + ( m1+m2 )*(g*sin(a)) + (m2*L1*b*b)*sine_val*cos_val - m2*g*sin(c)*cos_val  )  /  (L1*(-m1 - m2 + m2*cos_val*cos_val));
}

double second_equation(double a, double b, double c, double d)
{
	double sine_val = sin(a-c);
	double cos_val = cos(a-c);
	return  (	cos_val * (((m1+m2)*g*sin(a)) + (m2*L2*d*d*sine_val))  +  (m1+m2)*( -g*sin(c) + L1*b*b*sine_val) )  /  (L2*( m1 + m2 - m2*cos_val*cos_val));
}

#endif
