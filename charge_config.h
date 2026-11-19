#ifndef CHARGE_CONFIG_H
#define CHARGE_CONFIG_H

#define MAX_CHARGES 100
#define MAX_GRID_SIZE 100

/*
 * Shared data passed from Module 1 to Modules 2, 3, 4, 5, etc.
 *
 * The original Module 1 variable names are retained:
 * n, q, x, y, and gridSize.
 */
typedef struct
{
    int n;
    double q[MAX_CHARGES];
    double x[MAX_CHARGES];
    double y[MAX_CHARGES];
    int gridSize;
} ChargeConfig;

/* Module 1: read and validate the charge configuration. */
int readChargeConfiguration(ChargeConfig *config);

/* Module 1: display the configuration for checking by the user. */
void displayChargeConfiguration(const ChargeConfig *config);

#endif
