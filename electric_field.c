/*
 * Module 2: Electric Field Calculator
 * Course: IC103 Programming Project
 * File: electric_field.c
 *
 * Description:
 * Single-file implementation of Module 2 integrated with the new ChargeConfig
 * interface from Module 1 (charge_config.h).
 *
 * Calculates individual electric fields, resultant field, external field
 * superposition, reference point recalculations, and shortest distance for
 * target field magnitude. All shared structures, variables and functions are
 * public (non-static) so other modules can access them directly.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "charge_config.h"

/* Physical constants for dry air (20 deg C, 101.325 kPa) */
#define EPS0 8.8541878188e-12      /* Vacuum permittivity in F/m */
#define EPS_R 1.0005364            /* Relative permittivity of dry air */
#define PI 3.14159265358979323846  /* Pi constant */

/*
 * Structure to store electric field results:
 * Ex, Ey: field components in N/C
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

/* Globally accessible array for individual charge results (0-indexed: Charge 1 at index 0) */
ElectricFieldResult chargeField[MAX_CHARGES];

/* Globally accessible variables for resultant, external, and net fields */
ElectricFieldResult resultantField;
ElectricFieldResult externalField;
ElectricFieldResult netField;

/* Current reference point coordinates (in metres) */
double ref_x = 0.0;
double ref_y = 0.0;

/* Stored charge configuration received from Module 1 */
ChargeConfig currentConfig;

/* Function prototypes (all public / non-static) */
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
 * Returns -1.0 if magnitude is zero (direction is undefined).
 */
double calculate_angle(double ex, double ey)
{
    double mag = hypot(ex, ey);
    if (mag <= 1e-15)
    {
        return -1.0; /* Undefined direction for zero field vector */
    }

    double angle = atan2(ey, ex) * (180.0 / PI);
    if (angle < 0.0)
    {
        angle += 360.0;
    }
    if (angle >= 360.0)
    {
        angle -= 360.0;
    }
    if (fabs(angle) < 1e-12)
    {
        angle = 0.0;
    }
    return angle;
}

/*
 * Updates the net electric field by superimposing charge resultant and external field.
 */
void update_net_field(void)
{
    if (resultantField.defined)
    {
        netField.Ex = resultantField.Ex + externalField.Ex;
        netField.Ey = resultantField.Ey + externalField.Ey;
        netField.magnitude = hypot(netField.Ex, netField.Ey);
        netField.angle = calculate_angle(netField.Ex, netField.Ey);
        netField.defined = 1;
    }
    else if (externalField.defined)
    {
        netField.Ex = externalField.Ex;
        netField.Ey = externalField.Ey;
        netField.magnitude = externalField.magnitude;
        netField.angle = externalField.angle;
        netField.defined = 1;
    }
    else
    {
        netField.Ex = 0.0;
        netField.Ey = 0.0;
        netField.magnitude = 0.0;
        netField.angle = -1.0;
        netField.defined = 0;
    }
}

/*
 * Calculates individual electric fields and charge-only resultant for charges at (rx, ry).
 * Updates chargeField[], resultantField, and netField.
 */
void calculate_electric_field(const ChargeConfig *config, double rx, double ry)
{
    if (config == NULL || config->n <= 0)
    {
        resultantField.defined = 0;
        netField.defined = 0;
        return;
    }

    currentConfig = *config;
    ref_x = rx;
    ref_y = ry;

    double k = calculate_coulomb_constant();
    double total_ex = 0.0;
    double total_ey = 0.0;
    int valid_count = 0;

    for (int i = 0; i < currentConfig.n; i++)
    {
        /* Check if charge coincides with reference point (distance is zero) */
        if (currentConfig.x[i] == ref_x && currentConfig.y[i] == ref_y)
        {
            chargeField[i].Ex = 0.0;
            chargeField[i].Ey = 0.0;
            chargeField[i].magnitude = 0.0;
            chargeField[i].angle = -1.0;
            chargeField[i].defined = 0;
            continue;
        }

        double dx = ref_x - currentConfig.x[i];
        double dy = ref_y - currentConfig.y[i];
        double r = hypot(dx, dy);

        /* Charge q is directly in Coulombs (C) */
        double q_coulombs = currentConfig.q[i];

        /* E = k * q / r^2 along r vector (dx/r, dy/r) => k * q * dx / r^3 */
        double r3 = r * r * r;
        chargeField[i].Ex = (k * q_coulombs * dx) / r3;
        chargeField[i].Ey = (k * q_coulombs * dy) / r3;
        chargeField[i].magnitude = hypot(chargeField[i].Ex, chargeField[i].Ey);
        chargeField[i].angle = calculate_angle(chargeField[i].Ex, chargeField[i].Ey);
        chargeField[i].defined = 1;

        total_ex += chargeField[i].Ex;
        total_ey += chargeField[i].Ey;
        valid_count++;
    }

    /* Charge-only resultant field */
    resultantField.Ex = total_ex;
    resultantField.Ey = total_ey;

    if (valid_count == 0)
    {
        resultantField.magnitude = 0.0;
        resultantField.angle = -1.0;
        resultantField.defined = 0;
    }
    else
    {
        resultantField.magnitude = hypot(total_ex, total_ey);
        resultantField.angle = calculate_angle(total_ex, total_ey);
        resultantField.defined = 1;
    }

    /* Superimpose with external field to update netField */
    update_net_field();
}

/*
 * Public function: Changes reference point and recalculates all fields.
 * Preserves the existing external field.
 */
void change_reference_point(double new_x, double new_y)
{
    calculate_electric_field(&currentConfig, new_x, new_y);
}

/*
 * Public function: Sets or replaces the external electric field (Ex, Ey in N/C).
 */
void set_external_field(double ex, double ey)
{
    externalField.Ex = ex;
    externalField.Ey = ey;
    externalField.magnitude = hypot(ex, ey);
    externalField.angle = calculate_angle(ex, ey);
    externalField.defined = (externalField.magnitude > 1e-15) ? 1 : 0;
    update_net_field();
}

/*
 * Public function: Clears the external electric field and resets its components to zero.
 */
void clear_external_field(void)
{
    externalField.Ex = 0.0;
    externalField.Ey = 0.0;
    externalField.magnitude = 0.0;
    externalField.angle = -1.0;
    externalField.defined = 0;
    update_net_field();
}

/*
 * Public function: Calculates shortest radial distance (in metres) from selected
 * charge (0-based index) where its electric field magnitude equals desired_magnitude.
 * Formula: E = k * |q| / r^2  =>  r = sqrt(k * |q| / E).
 * Excludes other charges and external field. Returns -1.0 if invalid.
 */
double calculate_shortest_distance(int charge_index, double desired_magnitude)
{
    if (charge_index < 0 || charge_index >= currentConfig.n || desired_magnitude <= 0.0)
    {
        return -1.0;
    }

    double q_coulombs = fabs(currentConfig.q[charge_index]);
    if (q_coulombs <= 0.0)
    {
        return -1.0;
    }

    double k = calculate_coulomb_constant();
    return sqrt((k * q_coulombs) / desired_magnitude);
}

/*
 * Public function: Displays formatted electric field results to the console.
 */
void display_results(void)
{
    printf("\n===============================================================\n");
    printf("                  ELECTRIC FIELD RESULTS                       \n");
    printf("===============================================================\n");
    printf("Reference Point: (%.4f, %.4f) m\n\n", ref_x, ref_y);
    printf("--- Individual Charge Contributions ---\n");

    for (int i = 0; i < currentConfig.n; i++)
    {
        printf("Charge %d (q = %.6g C at (%.4f, %.4f) m):\n",
               i + 1, currentConfig.q[i], currentConfig.x[i], currentConfig.y[i]);

        if (chargeField[i].defined)
        {
            printf("  Ex = %+.4e N/C\n", chargeField[i].Ex);
            printf("  Ey = %+.4e N/C\n", chargeField[i].Ey);
            printf("  Magnitude = %.4e N/C\n", chargeField[i].magnitude);
            if (chargeField[i].angle >= 0.0)
            {
                printf("  Angle = %.2f degrees\n\n", chargeField[i].angle);
            }
            else
            {
                printf("  Angle = Undefined (zero magnitude)\n\n");
            }
        }
        else
        {
            printf("  Field is Undefined (charge coincides with reference point)\n\n");
        }
    }

    printf("--- Charge-only Resultant Field ---\n");
    if (resultantField.defined)
    {
        printf("  Ex = %+.4e N/C\n", resultantField.Ex);
        printf("  Ey = %+.4e N/C\n", resultantField.Ey);
        printf("  Magnitude = %.4e N/C\n", resultantField.magnitude);
        if (resultantField.angle >= 0.0)
        {
            printf("  Angle = %.2f degrees\n", resultantField.angle);
        }
        else
        {
            printf("  Angle = Undefined (zero magnitude)\n");
        }
    }
    else
    {
        printf("  Field is Undefined (all charges coincident with reference point)\n");
    }

    if (externalField.defined)
    {
        printf("\n--- External Electric Field ---\n");
        printf("  Ex = %+.4e N/C\n", externalField.Ex);
        printf("  Ey = %+.4e N/C\n", externalField.Ey);
        printf("  Magnitude = %.4e N/C\n", externalField.magnitude);
        if (externalField.angle >= 0.0)
        {
            printf("  Angle = %.2f degrees\n", externalField.angle);
        }
        else
        {
            printf("  Angle = Undefined\n");
        }

        printf("\n--- Final Net Field (Charges + External) ---\n");
        if (netField.defined)
        {
            printf("  Ex = %+.4e N/C\n", netField.Ex);
            printf("  Ey = %+.4e N/C\n", netField.Ey);
            printf("  Magnitude = %.4e N/C\n", netField.magnitude);
            if (netField.angle >= 0.0)
            {
                printf("  Angle = %.2f degrees\n", netField.angle);
            }
            else
            {
                printf("  Angle = Undefined\n");
            }
        }
        else
        {
            printf("  Net Field is Undefined\n");
        }
    }
    printf("===============================================================\n");
}

/*
 * Main public entry point for Module 2:
 * Receives charge configuration directly from Module 1.
 * Prompts for initial reference point and enters the repeating interactive menu.
 */
void run_electric_field_calculator(const ChargeConfig *config)
{
    if (config == NULL || config->n <= 0)
    {
        printf("Error: Invalid charge configuration received.\n");
        return;
    }

    /* Initialize external field to zero */
    clear_external_field();

    /* Flush any leftover characters in stdin */
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF)
    {
        /* discard */
    }

    double rx = 0.0;
    double ry = 0.0;
    printf("\n[Module 2] Enter initial reference point X coordinate (m): ");
    if (scanf("%lf", &rx) != 1) rx = 0.0;
    printf("[Module 2] Enter initial reference point Y coordinate (m): ");
    if (scanf("%lf", &ry) != 1) ry = 0.0;

    /* Initial calculation and display */
    calculate_electric_field(config, rx, ry);
    display_results();

    /* Repeating interactive menu */
    while (1)
    {
        printf("\n=================================================\n");
        printf("             ELECTRIC FIELD CALCULATOR\n");
        printf("=================================================\n\n");
        printf("Current Reference Point: (%.4f, %.4f)\n\n", ref_x, ref_y);
        printf("1. Display Electric Field Results\n");
        printf("2. Change Reference Point\n");
        printf("3. Set / Change External Electric Field\n");
        printf("4. Clear External Electric Field\n");
        printf("5. Calculate Shortest Distance\n");
        printf("6. Exit Module 2\n\n");
        printf("Enter your choice: ");

        int choice;
        if (scanf("%d", &choice) != 1)
        {
            /* Handle non-numeric input or EOF to avoid infinite loops */
            while ((ch = getchar()) != '\n' && ch != EOF)
            {
                /* discard */
            }
            if (ch == EOF)
            {
                printf("\nInput ended (EOF). Exiting Module 2.\n");
                break;
            }
            printf("Invalid input. Please enter a number between 1 and 6.\n");
            continue;
        }

        if (choice == 1)
        {
            display_results();
        }
        else if (choice == 2)
        {
            double new_x = 0.0, new_y = 0.0;
            printf("\nEnter new reference point X coordinate (m): ");
            if (scanf("%lf", &new_x) != 1) new_x = 0.0;
            printf("Enter new reference point Y coordinate (m): ");
            if (scanf("%lf", &new_y) != 1) new_y = 0.0;

            change_reference_point(new_x, new_y);
            display_results();
        }
        else if (choice == 3)
        {
            double ex = 0.0, ey = 0.0;
            printf("\nEnter external electric field Ex (N/C): ");
            if (scanf("%lf", &ex) != 1) ex = 0.0;
            printf("Enter external electric field Ey (N/C): ");
            if (scanf("%lf", &ey) != 1) ey = 0.0;

            set_external_field(ex, ey);
            printf("\nExternal electric field updated.\n");
            display_results();
        }
        else if (choice == 4)
        {
            clear_external_field();
            printf("\nExternal electric field cleared.\n");
            display_results();
        }
        else if (choice == 5)
        {
            printf("\nAvailable charges:\n");
            for (int i = 0; i < currentConfig.n; i++)
            {
                printf("  Charge %d: q = %.6g C at (%.4f, %.4f) m\n",
                       i + 1, currentConfig.q[i], currentConfig.x[i], currentConfig.y[i]);
            }

            int charge_id;
            printf("Select charge number (1 to %d): ", currentConfig.n);
            if (scanf("%d", &charge_id) != 1 || charge_id < 1 || charge_id > currentConfig.n)
            {
                printf("Invalid charge selection.\n");
                continue;
            }

            double desired_mag;
            printf("Enter desired electric field magnitude (N/C): ");
            if (scanf("%lf", &desired_mag) != 1 || desired_mag <= 0.0)
            {
                printf("Desired field magnitude must be positive.\n");
                continue;
            }

            double dist = calculate_shortest_distance(charge_id - 1, desired_mag);
            if (dist < 0.0)
            {
                printf("Cannot calculate distance: charge value is zero.\n");
            }
            else
            {
                printf("\n--- Shortest Distance Result ---\n");
                printf("Selected Charge: Charge %d (q = %.6g C at (%.4f, %.4f) m)\n",
                       charge_id, currentConfig.q[charge_id - 1],
                       currentConfig.x[charge_id - 1], currentConfig.y[charge_id - 1]);
                printf("Desired Field Magnitude: %.4e N/C\n", desired_mag);
                printf("Shortest Distance from Charge: %.4f m\n", dist);
            }
        }
        else if (choice == 6)
        {
            printf("\nExiting Electric Field Calculator (Module 2). Returning control to caller.\n");
            break;
        }
        else
        {
            printf("Invalid choice. Please select an option between 1 and 6.\n");
        }
    }
}
