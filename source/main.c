#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define pi 3.14159265359

int main(int argc, char *argv[])
{
    int user_response;
	int waiting_variable;

    printf("\nInput 0 for inversion time graph, or 1 for error graph.");
    scanf("%d", &user_response);

    if(user_response == 0)
    {
        graph_first_inversion_time(argv);
    }
    else if(user_response == 1)
    {
        error_graph(argv);
    }
    else
    {
        printf("\nInvalid response. Exiting.\n");
    }


    return 0;
}