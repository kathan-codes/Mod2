#include <stdio.h>



int main()
{
    int n, gridSize;
    
    // Number of charges
    printf("Enter number of charges: ");
    scanf("%d", &n);

    double q[n];
    double x[n], y[n];


    // Charge values and coordinates
    for (int i = 0; i < n; i++) 
    {
        printf("\nCharge %d:\n", i + 1);

        printf("Enter charge value: ");
        scanf("%lf", &q[i]);

        printf("Enter x coordinate: ");
        scanf("%lf", &x[i]);

        printf("Enter y coordinate: ");
        scanf("%lf", &y[i]);
    }

    // Grid size
    printf("\nEnter grid size (1 to %d): ", MAX_GRID_SIZE);
    scanf("%d", &gridSize);

    if (gridSize < 1 || gridSize > MAX_GRID_SIZE)
     {
        printf("Invalid grid size.\n");
        return 1;
    }

    return 0;
}