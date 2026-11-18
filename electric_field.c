/*
 * Module 2: Electric Field Calculator
 * Course: IC103 Programming Project
 * File: electric_field.c
 *
 * Description:
 * Rebuilt single-file implementation of Module 2.
 * Calculates individual electric fields, resultant field, external field
 * superposition, reference point recalculations, and shortest distance
 * for target field magnitude. All shared arrays and functions are public (non-static)
 * so Modules 3, 4 and later can access them directly.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Physical constants for dry air (20 deg C, 101.325 kPa) */
#define EPS0 8.8541878188e-12      /* Vacuum permittivity in F/m */
#define EPS_R 1.0005364            /* Relative permittivity of dry air */
#define PI 3.14159265358979323846  /* Pi */
#define MAX_CHARGES 100            /* Maximum supported charges */

/* Globally accessible arrays for individual charge results (0-indexed: Charge 1 at index 0) */
double Ex[MAX_CHARGES];
double Ey[MAX_CHARGES];
double E_magnitude[MAX_CHARGES];
double E_angle[MAX_CHARGES];
int E_defined[MAX_CHARGES];

/* Stored charge configuration from Module 1 */
int num_charges = 0;
double charges_q[MAX_CHARGES];
double charges_x[MAX_CHARGES];
double charges_y[MAX_CHARGES];

/* Current reference point coordinates (in metres) */
double ref_x = 0.0;
double ref_y = 0.0;

/* Globally accessible variables for final charge-only resultant field */
double final_Ex = 0.0;
double final_Ey = 0.0;
double final_magnitude = 0.0;
double final_angle = 0.0;
int final_defined = 0;

/* Globally accessible variables for external electric field */
double ext_Ex = 0.0;
double ext_Ey = 0.0;
int ext_active = 0;

/* Globally accessible variables for final net electric field (charges + external) */
double net_Ex = 0.0;
double net_Ey = 0.0;
double net_magnitude = 0.0;
double net_angle = 0.0;
int net_defined = 0;

/* Function prototypes (all public / non-static) */
double calculate_coulomb_constant(void);
double calculate_angle(double ex, double ey, int *is_defined);
void update_net_field(void);
void calculate_electric_field(int n, const double q[], const double x[], const double y[], double rx, double ry);
void change_reference_point(double new_x, double new_y);
void set_external_field(double ex, double ey);
void set_external_field_polar(double magnitude, double angle_deg);
void clear_external_field(void);
double calculate_shortest_distance(int charge_index, double desired_magnitude);
double calculate_shortest_distance_by_id(int charge_id, double desired_magnitude);
void display_results(void);
void run_electric_field_calculator(int n, const double q[], const double x[], const double y[]);


/*
 * Calculates Coulomb's constant in dry air:
 * k = 1 / (4 * PI * EPS0 * EPS_R)
 */
double calculate_coulomb_constant(void)
{
    return 1.0 / (4.0 * PI * EPS0 * EPS_R);
}

/*
 * Calculates counterclockwise angle in degrees [0, 360) from +X axis.
 * Sets *is_defined to 0 if magnitude is zero (undefined direction).
 */
double calculate_angle(double ex, double ey, int *is_defined)
{
    double mag = hypot(ex, ey);
    if (mag <= 1e-15)
    {
        *is_defined = 0;
        return 0.0;
    }

    *is_defined = 1;
    double angle = atan2(ey, ex) * (180.0 / PI);
    if (angle < 0.0)
    {
        angle += 360.0;
    }
    if (angle >= 360.0)
    {
        angle -= 360.0;
    }
    return angle;
}

/*
 * Updates the net electric field by superimposing charge resultant and external field.
 */
void update_net_field(void)
{
    if (final_defined)
    {
        net_Ex = final_Ex + ext_Ex;
        net_Ey = final_Ey + ext_Ey;
        net_magnitude = hypot(net_Ex, net_Ey);
        net_angle = calculate_angle(net_Ex, net_Ey, &net_defined);
    }
    else if (ext_active)
    {
        net_Ex = ext_Ex;
        net_Ey = ext_Ey;
        net_magnitude = hypot(net_Ex, net_Ey);
        net_angle = calculate_angle(net_Ex, net_Ey, &net_defined);
    }
    else
    {
        net_Ex = 0.0;
        net_Ey = 0.0;
        net_magnitude = 0.0;
        net_angle = 0.0;
        net_defined = 0;
    }
}

/*
 * Calculates individual electric fields and resultant for n charges at (rx, ry).
 * Updates Ex[], Ey[], E_magnitude[], E_angle[], E_defined[],
 * and updates final resultant and net fields.
 */
void calculate_electric_field(int n, const double q[], const double x[], const double y[], double rx, double ry)
{
    if (n <= 0 || q == NULL || x == NULL || y == NULL)
    {
        num_charges = 0;
        final_defined = 0;
        net_defined = 0;
        return;
    }

    num_charges = (n > MAX_CHARGES) ? MAX_CHARGES : n;
    ref_x = rx;
    ref_y = ry;

    for (int i = 0; i < num_charges; i++)
    {
        charges_q[i] = q[i];
        charges_x[i] = x[i];
        charges_y[i] = y[i];
    }

    double k = calculate_coulomb_constant();
    double total_ex = 0.0;
    double total_ey = 0.0;
    int valid_charge_count = 0;

    for (int i = 0; i < num_charges; i++)
    {
        /* Coincident charge: field is undefined at the point charge location */
        if (charges_x[i] == ref_x && charges_y[i] == ref_y)
        {
            Ex[i] = 0.0;
            Ey[i] = 0.0;
            E_magnitude[i] = 0.0;
            E_angle[i] = 0.0;
            E_defined[i] = 0;
            continue;
        }

        double dx = ref_x - charges_x[i];
        double dy = ref_y - charges_y[i];
        double r = hypot(dx, dy);

        /* Convert microcoulombs to coulombs */
        double q_coulombs = charges_q[i] * 1e-6;

        /* E = k * q / r^2 in direction of r vector (dx/r, dy/r) => k * q * dx / r^3 */
        double r3 = r * r * r;
        Ex[i] = (k * q_coulombs * dx) / r3;
        Ey[i] = (k * q_coulombs * dy) / r3;
        E_magnitude[i] = hypot(Ex[i], Ey[i]);
        E_angle[i] = calculate_angle(Ex[i], Ey[i], &E_defined[i]);

        total_ex += Ex[i];
        total_ey += Ey[i];
        valid_charge_count++;
    }

    /* Final resultant field due to charges */
    final_Ex = total_ex;
    final_Ey = total_ey;
    if (valid_charge_count == 0)
    {
        final_magnitude = 0.0;
        final_angle = 0.0;
        final_defined = 0;
    }
    else
    {
        final_magnitude = hypot(final_Ex, final_Ey);
        final_angle = calculate_angle(final_Ex, final_Ey, &final_defined);
    }

    /* Update net field with external field */
    update_net_field();
}

/*
 * Public function: Changes reference point and recalculates fields
 * using stored charges. Preserves existing external field.
 */
void change_reference_point(double new_x, double new_y)
{
    calculate_electric_field(num_charges, charges_q, charges_x, charges_y, new_x, new_y);
}

/*
 * Public function: Sets or changes external electric field by component values.
 */
void set_external_field(double ex, double ey)
{
    ext_Ex = ex;
    ext_Ey = ey;
    ext_active = 1;
    update_net_field();
}

/*
 * Public function: Sets or changes external electric field by magnitude (N/C)
 * and direction angle (degrees CCW from +X axis).
 */
void set_external_field_polar(double magnitude, double angle_deg)
{
    if (magnitude <= 0.0)
    {
        clear_external_field();
        return;
    }

    double rad = angle_deg * (PI / 180.0);
    ext_Ex = magnitude * cos(rad);
    ext_Ey = magnitude * sin(rad);
    ext_active = 1;
    update_net_field();
}

/*
 * Public function: Clears the external electric field.
 */
void clear_external_field(void)
{
    ext_Ex = 0.0;
    ext_Ey = 0.0;
    ext_active = 0;
    update_net_field();
}

/*
 * Public function: Calculates shortest radial distance (in metres) from selected
 * charge (0-based index) where electric field magnitude equals desired_magnitude.
 * Formula: E = k * |q| / r^2  =>  r = sqrt(k * |q| / E).
 * Excludes external field and other charges. Returns -1.0 if invalid.
 */
double calculate_shortest_distance(int charge_index, double desired_magnitude)
{
    if (charge_index < 0 || charge_index >= num_charges || desired_magnitude <= 0.0)
    {
        return -1.0;
    }

    double q_coulombs = fabs(charges_q[charge_index]) * 1e-6;
    if (q_coulombs <= 0.0)
    {
        return -1.0;
    }

    double k = calculate_coulomb_constant();
    return sqrt((k * q_coulombs) / desired_magnitude);
}

/*
 * Convenience function: Same as calculate_shortest_distance using 1-based charge ID (1, 2, ...).
 */
double calculate_shortest_distance_by_id(int charge_id, double desired_magnitude)
{
    return calculate_shortest_distance(charge_id - 1, desired_magnitude);
}

/*
 * Public function: Prints summary of calculated electric fields to stdout.
 */
void display_results(void)
{
    printf("\n===============================================================\n");
    printf("                  ELECTRIC FIELD CALCULATOR (MODULE 2)         \n");
    printf("===============================================================\n");
    printf("Reference Point: (%.4f, %.4f) m\n\n", ref_x, ref_y);
    printf("--- Individual Charge Contributions ---\n");

    for (int i = 0; i < num_charges; i++)
    {
        if (E_defined[i])
        {
            printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n",
                   i + 1, charges_q[i], charges_x[i], charges_y[i]);
            printf("  Ex = %+.4e N/C\n", Ex[i]);
            printf("  Ey = %+.4e N/C\n", Ey[i]);
            printf("  Magnitude = %.4e N/C\n", E_magnitude[i]);
            printf("  Angle = %.2f degrees\n\n", E_angle[i]);
        }
        else
        {
            printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n",
                   i + 1, charges_q[i], charges_x[i], charges_y[i]);
            printf("  Field Undefined (coincides with reference point)\n\n");
        }
    }

    printf("--- Charge-only Resultant Field ---\n");
    if (final_defined)
    {
        printf("  Ex = %+.4e N/C\n", final_Ex);
        printf("  Ey = %+.4e N/C\n", final_Ey);
        printf("  Magnitude = %.4e N/C\n", final_magnitude);
        printf("  Angle = %.2f degrees\n", final_angle);
    }
    else
    {
        printf("  Field is Undefined (all charges coincident or undefined)\n");
    }

    if (ext_active)
    {
        printf("\n--- External Electric Field ---\n");
        printf("  Ex = %+.4e N/C\n", ext_Ex);
        printf("  Ey = %+.4e N/C\n", ext_Ey);
        printf("  Magnitude = %.4e N/C\n", hypot(ext_Ex, ext_Ey));

        printf("\n--- Final Net Field (Charges + External) ---\n");
        if (net_defined)
        {
            printf("  Ex = %+.4e N/C\n", net_Ex);
            printf("  Ey = %+.4e N/C\n", net_Ey);
            printf("  Magnitude = %.4e N/C\n", net_magnitude);
            printf("  Angle = %.2f degrees\n", net_angle);
        }
        else
        {
            printf("  Net Field is Undefined\n");
        }
    }
    printf("===============================================================\n");
}

/*
 * Public integration entry point if called by Module 1.
 * Prompts user for initial reference point and calculates fields.
 */
void run_electric_field_calculator(int n, const double q[], const double x[], const double y[])
{
    double rx = 0.0;
    double ry = 0.0;
    printf("\n[Module 2] Enter reference point X coordinate (m): ");
    if (scanf("%lf", &rx) != 1) rx = 0.0;
    printf("[Module 2] Enter reference point Y coordinate (m): ");
    if (scanf("%lf", &ry) != 1) ry = 0.0;

    calculate_electric_field(n, q, x, y, rx, ry);
    display_results();
}
