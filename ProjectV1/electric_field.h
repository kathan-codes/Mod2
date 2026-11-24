#ifndef ELECTRIC_FIELD_H
#define ELECTRIC_FIELD_H

#include "charge_config.h"

/* Physical constants for dry air (20 deg C, 101.325 kPa) */
#define EPS0 8.8541878188e-12      /* Vacuum permittivity in F/m */
#define EPS_R_AIR 1.0005364        /* Relative permittivity of dry air */
#define PI 3.14159265358979323846  /* Pi constant */

/*
 * Structure to store electric field results:
 * Ex, Ey: components in N/C
 * magnitude: field magnitude in N/C
 * angle: counterclockwise angle [0, 360) deg from +X axis (-1 if undefined)
 * defined: 1 if field is valid, 0 if undefined (e.g. coincident charge)
 */
typedef struct
{
    double Ex;
    double Ey;
    double magnitude;
    double angle;
    int defined;
} ElectricFieldResult;

/* Globally accessible result variables */
extern ElectricFieldResult chargeField[MAX_CHARGES];
extern ElectricFieldResult resultantField;
extern ElectricFieldResult externalField;
extern ElectricFieldResult netField;

/* Current reference-point coordinates in metres */
extern double referenceX;
extern double referenceY;

/* Public function prototypes */
double calculate_coulomb_constant(void);
double calculate_angle(double ex, double ey);
void update_net_field(void);
void calculate_electric_field(const ChargeConfig *config, double rx, double ry);
void change_reference_point(double new_x, double new_y);
void set_external_field(double ex, double ey);
void clear_external_field(void);
double calculate_shortest_distance(int charge_index, double desired_magnitude);
void display_results(void);
void run_electric_field_calculator(const ChargeConfig *config);

#endif /* ELECTRIC_FIELD_H */
