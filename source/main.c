#include "flip_time_graph.h"
#include "error_graph.h"

int main(int argc, char *argv[])
{
    int user_response;

    printf("\nInput 0 for inversion time graph, or 1 for error graph.");
    scanf("%d", &user_response);

    if(user_response == 0)
    {
        flip_time_graph(argv);
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