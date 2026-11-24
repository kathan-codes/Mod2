/*
 * Module 2: Electric Field Calculator
 * Course: IC103 Programming Project
 * File: electric_field.c
 *
 * Description:
 * Implementation of Module 2 Electric Field Calculator.
 * Calculates individual electric fields, charge-only resultant,
 * uniform background field superposition, reference point updates,
 * and shortest distance to a target field circle (Output 4).
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "electric_field.h"

/* Globally accessible result variables defined here */
ElectricFieldResult chargeField[MAX_CHARGES];
ElectricFieldResult resultantField;
ElectricFieldResult externalField;
ElectricFieldResult netField;

/* Current reference-point coordinates in metres */
double referenceX = 0.0;
double referenceY = 0.0;

/* Number of charges coincident with current reference point (skipped from sum) */
int skippedChargesCount = 0;

/* Stored charge configuration received from Module 1 */
ChargeConfig currentConfig;

static void reset_field_result(ElectricFieldResult *field)
{
    if (field == NULL)
    {
        return;
    }

    field->Ex = 0.0;
    field->Ey = 0.0;
    field->magnitude = 0.0;
    field->angle = -1.0;
    field->defined = 0;
}

static void reset_calculation_results(void)
{
    for (int i = 0; i < MAX_CHARGES; i++)
    {
        reset_field_result(&chargeField[i]);
    }
    reset_field_result(&resultantField);
    reset_field_result(&netField);
}

static int valid_charge_config(const ChargeConfig *config)
{
    if (config == NULL || config->n < 1 || config->n > MAX_CHARGES ||
        config->gridSize < 1 || config->gridSize > MAX_GRID_SIZE)
    {
        return 0;
    }

    for (int i = 0; i < config->n; i++)
    {
        if (!isfinite(config->q[i]) || !isfinite(config->x[i]) ||
            !isfinite(config->y[i]))
        {
            return 0;
        }
    }
    return 1;
}

/* Read one finite number, retrying bad text and stopping cleanly at EOF. */
static int read_finite_double(const char *prompt, double *value)
{
    int ch;

    for (;;)
    {
        printf("%s", prompt);
        if (scanf("%lf", value) == 1)
        {
            if (isfinite(*value))
            {
                return 1;
            }
            printf("Value must be finite. Please try again.\n");
            continue;
        }

        while ((ch = getchar()) != '\n' && ch != EOF)
        {
            /* discard invalid input */
        }
        if (ch == EOF)
        {
            return 0;
        }
        printf("Invalid number. Please enter a finite value.\n");
    }
}

/*
 * Calculates Coulomb's constant in dry air:
 * k = 1 / (4 * PI * EPS0 * EPS_R_AIR)
 */
double calculate_coulomb_constant(void)
{
    return 1.0 / (4.0 * PI * EPS0 * EPS_R_AIR);
}

/*
 * Calculates counterclockwise angle in degrees [0, 360) from +X axis.
 * Returns -1.0 if magnitude is zero (direction is undefined).
 */
double calculate_angle(double ex, double ey)
{
    if (!isfinite(ex) || !isfinite(ey))
    {
        return -1.0;
    }

    double mag = hypot(ex, ey);
    if (!isfinite(mag) || mag <= 1e-15)
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
 * If charge resultant is undefined/zero but background field is active, net field is background field.
 */
void update_net_field(void)
{
    if (resultantField.defined)
    {
        netField.Ex = resultantField.Ex + externalField.Ex;
        netField.Ey = resultantField.Ey + externalField.Ey;
        netField.magnitude = hypot(netField.Ex, netField.Ey);
        if (!isfinite(netField.Ex) || !isfinite(netField.Ey) ||
            !isfinite(netField.magnitude))
        {
            reset_field_result(&netField);
            return;
        }

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
        reset_field_result(&netField);
    }
}

/*
 * Calculates individual electric fields and charge-only resultant for charges at (rx, ry).
 * Handles coincident charges safely by skipping them and summing remaining charges.
 * Updates chargeField[], resultantField, and netField.
 */
void calculate_electric_field(const ChargeConfig *config, double rx, double ry)
{
    if (!valid_charge_config(config))
    {
        currentConfig.n = 0;
        currentConfig.gridSize = 0;
        referenceX = 0.0;
        referenceY = 0.0;
        skippedChargesCount = 0;
        reset_calculation_results();
        return;
    }

    currentConfig = *config;
    if (!isfinite(rx) || !isfinite(ry))
    {
        referenceX = 0.0;
        referenceY = 0.0;
        skippedChargesCount = 0;
        reset_calculation_results();
        return;
    }

    referenceX = rx;
    referenceY = ry;
    skippedChargesCount = 0;

    double k = calculate_coulomb_constant();
    double total_ex = 0.0;
    double total_ey = 0.0;
    int valid_count = 0;

    for (int i = 0; i < currentConfig.n; i++)
    {
        double dx = referenceX - currentConfig.x[i];
        double dy = referenceY - currentConfig.y[i];
        double r = hypot(dx, dy);

        /* Check if charge coincides with reference point (distance < 1e-9 m) */
        if (!isfinite(r) || r < 1e-9)
        {
            reset_field_result(&chargeField[i]);
            skippedChargesCount++;
            continue;
        }

        double q_coulombs = currentConfig.q[i];

        /* Compute k*q/r^2 by successive divisions to avoid forming r^3 */
        double field_scale = (k * (q_coulombs / r)) / r;
        double unit_x = dx / r;
        double unit_y = dy / r;
        chargeField[i].Ex = field_scale * unit_x;
        chargeField[i].Ey = field_scale * unit_y;
        chargeField[i].magnitude = hypot(chargeField[i].Ex, chargeField[i].Ey);

        if (!isfinite(field_scale) || !isfinite(chargeField[i].Ex) ||
            !isfinite(chargeField[i].Ey) || !isfinite(chargeField[i].magnitude))
        {
            reset_field_result(&chargeField[i]);
            continue;
        }

        chargeField[i].angle = calculate_angle(chargeField[i].Ex, chargeField[i].Ey);
        chargeField[i].defined = 1;

        total_ex += chargeField[i].Ex;
        total_ey += chargeField[i].Ey;
        valid_count++;
    }

    /* Charge-only resultant field */
    if (valid_count == 0)
    {
        resultantField.Ex = 0.0;
        resultantField.Ey = 0.0;
        resultantField.magnitude = 0.0;
        resultantField.angle = -1.0;
        resultantField.defined = 0;
    }
    else
    {
        resultantField.Ex = total_ex;
        resultantField.Ey = total_ey;
        resultantField.magnitude = hypot(total_ex, total_ey);
        if (!isfinite(resultantField.magnitude))
        {
            reset_field_result(&resultantField);
        }
        else
        {
            resultantField.angle = calculate_angle(total_ex, total_ey);
            resultantField.defined = 1;
        }
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
 * Public function: Sets or replaces the external electric field using components (Ex, Ey in N/C).
 */
void set_external_field(double ex, double ey)
{
    if (!isfinite(ex) || !isfinite(ey))
    {
        clear_external_field();
        return;
    }

    double magnitude = hypot(ex, ey);
    if (!isfinite(magnitude) || magnitude <= 1e-15)
    {
        clear_external_field();
        return;
    }

    externalField.Ex = ex;
    externalField.Ey = ey;
    externalField.magnitude = magnitude;
    externalField.angle = calculate_angle(ex, ey);
    externalField.defined = 1;
    update_net_field();
}

/*
 * Public function: Sets or replaces the external electric field using magnitude (N/C)
 * and direction (degrees, measured counterclockwise from the positive X-axis).
 * Entering magnitude <= 0 clears the external field.
 */
void set_external_field_polar(double magnitude, double angle_deg)
{
    if (!isfinite(magnitude) || magnitude <= 0.0 || !isfinite(angle_deg))
    {
        clear_external_field();
        return;
    }

    double rad = angle_deg * (PI / 180.0);
    double ex = magnitude * cos(rad);
    double ey = magnitude * sin(rad);
    set_external_field(ex, ey);
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
 * Public function: Calculates shortest distance from the current reference point
 * to a location where the selected charge's electric field magnitude equals desired_magnitude.
 * Formula:
 *   r_target = sqrt(k * |q| / E_target)
 *   d = sqrt((x_r - x_q)^2 + (y_r - y_q)^2)
 *   d_shortest = |d - r_target|
 * Excludes other charges and external field. Returns -1.0 if invalid.
 */
double calculate_shortest_distance(int charge_index, double desired_magnitude)
{
    if (charge_index < 0 || charge_index >= currentConfig.n ||
        !isfinite(desired_magnitude) || desired_magnitude <= 0.0)
    {
        return -1.0;
    }

    double q_coulombs = fabs(currentConfig.q[charge_index]);
    if (!isfinite(q_coulombs) || q_coulombs <= 0.0 ||
        !isfinite(referenceX) || !isfinite(referenceY))
    {
        return -1.0;
    }

    double k = calculate_coulomb_constant();
    double r_target = (sqrt(k) * sqrt(q_coulombs)) / sqrt(desired_magnitude);
    double dx = referenceX - currentConfig.x[charge_index];
    double dy = referenceY - currentConfig.y[charge_index];
    double d = hypot(dx, dy);

    if (!isfinite(r_target) || !isfinite(d))
    {
        return -1.0;
    }

    double distance = fabs(d - r_target);
    if (distance < 1e-12)
    {
        distance = 0.0;
    }
    return isfinite(distance) ? distance : -1.0;
}

/*
 * Public function: Displays formatted electric field results to the console.
 */
void display_results(void)
{
    printf("\n===============================================================\n");
    printf("                  ELECTRIC FIELD RESULTS                       \n");
    printf("===============================================================\n");
    printf("Reference Point: (%.4f, %.4f) m\n\n", referenceX, referenceY);
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
            printf("  Field is Undefined (charge coincides with reference point; excluded from sum)\n\n");
        }
    }

    if (skippedChargesCount > 0)
    {
        printf("Note: %d charge(s) coincident with the reference point were excluded from the resultant calculation.\n\n",
               skippedChargesCount);
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
        printf("  Charge-only Resultant: Undefined (all charges coincident with reference point)\n");
        printf("  Displayed zero sum is over 0 included charges; physical field is undefined at this point.\n");
    }

    if (externalField.defined)
    {
        printf("\n--- External Background Electric Field ---\n");
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

        printf("\n--- Final Net Field (Charges + Background) ---\n");
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
    if (!valid_charge_config(config))
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
    if (!read_finite_double("\n[Module 2] Enter initial reference point X coordinate (m): ", &rx) ||
        !read_finite_double("[Module 2] Enter initial reference point Y coordinate (m): ", &ry))
    {
        printf("\nInput ended. Exiting Module 2.\n");
        return;
    }

    /* Initial calculation and display */
    calculate_electric_field(config, rx, ry);
    display_results();

    /* Repeating interactive menu */
    while (1)
    {
        printf("\n===== ELECTRIC FIELD CALCULATOR =====\n\n");
        printf("1. Display Electric Field Results\n");
        printf("2. Change Reference Point\n");
        printf("3. Set / Change Background Electric Field (Enter magnitude 0 to clear)\n");
        printf("4. Clear Background Electric Field\n");
        printf("5. Find Shortest Distance for a Target Electric Field (Output 4)\n");
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
            printf("Invalid choice. Please enter a number between 1 and 6.\n");
            continue;
        }

        if (choice == 1)
        {
            display_results();
        }
        else if (choice == 2)
        {
            double new_x = 0.0, new_y = 0.0;
            if (!read_finite_double("\nEnter new reference point X coordinate (m): ", &new_x) ||
                !read_finite_double("Enter new reference point Y coordinate (m): ", &new_y))
            {
                printf("\nInput ended. Exiting Module 2.\n");
                return;
            }

            change_reference_point(new_x, new_y);
            display_results();
        }
        else if (choice == 3)
        {
            double mag = 0.0;
            if (!read_finite_double("\nEnter background field magnitude in N/C (0 to clear): ", &mag))
            {
                printf("\nInput ended. Exiting Module 2.\n");
                return;
            }

            if (mag < 0.0)
            {
                printf("Background field magnitude cannot be negative.\n");
                continue;
            }

            if (mag == 0.0)
            {
                clear_external_field();
                printf("\nBackground electric field cleared.\n");
                display_results();
                continue;
            }

            double angle_deg = 0.0;
            if (!read_finite_double("Enter background field direction in degrees (CCW from +X axis): ", &angle_deg))
            {
                printf("\nInput ended. Exiting Module 2.\n");
                return;
            }

            set_external_field_polar(mag, angle_deg);
            printf("\nBackground electric field updated successfully.\n");
            display_results();
        }
        else if (choice == 4)
        {
            clear_external_field();
            printf("\nBackground electric field cleared.\n");
            display_results();
        }
        else if (choice == 5)
        {
            /* Pre-check: Verify at least one nonzero charge exists */
            int nonzero_count = 0;
            for (int i = 0; i < currentConfig.n; i++)
            {
                if (fabs(currentConfig.q[i]) > 1e-15)
                {
                    nonzero_count++;
                }
            }

            if (nonzero_count == 0)
            {
                printf("\nOutput 4 cannot be performed: no nonzero charges are available.\n");
                continue;
            }

            printf("\nAvailable charges:\n");
            for (int i = 0; i < currentConfig.n; i++)
            {
                printf("  Charge %d: q = %.6g C at (%.4f, %.4f) m%s\n",
                       i + 1, currentConfig.q[i], currentConfig.x[i], currentConfig.y[i],
                       (fabs(currentConfig.q[i]) <= 1e-15) ? " [zero charge - ineligible]" : "");
            }

            int charge_id;
            printf("Select charge number (1 to %d): ", currentConfig.n);
            if (scanf("%d", &charge_id) != 1 || charge_id < 1 || charge_id > currentConfig.n)
            {
                printf("Invalid charge selection.\n");
                continue;
            }

            if (fabs(currentConfig.q[charge_id - 1]) <= 1e-15)
            {
                printf("Selected charge has zero value; cannot produce a target electric field.\n");
                continue;
            }

            double desired_mag;
            for (;;)
            {
                if (!read_finite_double("Enter desired electric field magnitude (N/C): ", &desired_mag))
                {
                    printf("\nInput ended. Exiting Module 2.\n");
                    return;
                }

                if (desired_mag < 0.0)
                {
                    printf("Invalid: Target field magnitude cannot be negative. Please try again.\n");
                    continue;
                }

                if (desired_mag == 0.0)
                {
                    printf("For a point charge, no finite matching location exists where the electric field magnitude is zero (field approaches zero only as distance approaches infinity).\n");
                    break;
                }
                break;
            }

            if (desired_mag == 0.0)
            {
                continue;
            }

            double dist = calculate_shortest_distance(charge_id - 1, desired_mag);
            if (dist < 0.0)
            {
                printf("Cannot calculate distance: computation error.\n");
            }
            else
            {
                printf("\n--- Shortest Distance Result (Output 4) ---\n");
                printf("Selected Charge: Charge %d (q = %.6g C at (%.4f, %.4f) m)\n",
                       charge_id, currentConfig.q[charge_id - 1],
                       currentConfig.x[charge_id - 1], currentConfig.y[charge_id - 1]);
                printf("Current Reference Point: (%.4f, %.4f) m\n", referenceX, referenceY);
                printf("Desired Field Magnitude: %.4e N/C\n", desired_mag);
                printf("Shortest Distance to Target-Field Location: %.4f m\n", dist);

                if (dist == 0.0)
                {
                    printf("Note: The current reference point lies exactly on the target-field circle.\n");
                }
                else
                {
                    double d = hypot(referenceX - currentConfig.x[charge_id - 1],
                                     referenceY - currentConfig.y[charge_id - 1]);
                    if (d < 1e-9)
                    {
                        printf("Note: The reference point coincides with the charge location; shortest distance equals target radius.\n");
                    }
                }
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
