/*
 * Module 2: Electric Field Calculator
 * Course: IC103 Programming Project
 * Medium: Dry air at 20 deg C, 101.325 kPa
 *
 * This module receives charge configuration data from Module 1
 * and calculates individual and resultant electric fields,
 * allows reference point updates, uniform background field superposition,
 * and nearest distance determination for a target field (Output 4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Physical constants for dry air approximation */
#define EPS0 8.8541878188e-12      /* Vacuum permittivity in F/m */
#define EPS_R 1.0005364            /* Relative permittivity of dry air */
#define PI 3.14159265358979323846  /* Pi constant */

/*
 * Calculate Coulomb's constant in air:
 * k = 1 / (4 * pi * eps0 * eps_r)
 */
static double get_coulomb_constant(void)
{
    return 1.0 / (4.0 * PI * EPS0 * EPS_R);
}

/*
 * Reads a line from stdin and safely parses a finite double.
 * Returns 1 on success, 0 on EOF.
 */
static int read_double(const char *prompt, double *val)
{
    char buffer[128];
    while (1)
    {
        printf("%s", prompt);
        fflush(stdout);
        if (!fgets(buffer, sizeof(buffer), stdin))
        {
            return 0; /* EOF */
        }
        if (sscanf(buffer, "%lf", val) == 1 && isfinite(*val))
        {
            return 1;
        }
        printf("Invalid input. Please enter a valid finite number.\n");
    }
}

/*
 * Reads a line from stdin and safely parses an integer.
 * Returns 1 on success, 0 on EOF.
 */
static int read_int(const char *prompt, int *val)
{
    char buffer[128];
    while (1)
    {
        printf("%s", prompt);
        fflush(stdout);
        if (!fgets(buffer, sizeof(buffer), stdin))
        {
            return 0; /* EOF */
        }
        if (sscanf(buffer, "%d", val) == 1)
        {
            return 1;
        }
        printf("Invalid input. Please enter a valid integer.\n");
    }
}

/*
 * Calculates direction angle in degrees [0, 360) counterclockwise from +X axis.
 * Sets *is_defined to 1 if direction is defined, or 0 if magnitude is zero.
 */
static double calculate_angle(double ex, double ey, double mag, int *is_defined)
{
    /* If field magnitude is zero or within zero tolerance, angle is undefined */
    if (mag <= 1e-12)
    {
        *is_defined = 0;
        return 0.0;
    }

    *is_defined = 1;

    /* Special case: Ex == 0 */
    if (fabs(ex) <= 1e-12)
    {
        if (ey > 0.0)
            return 90.0;
        else if (ey < 0.0)
            return 270.0;
        else
        {
            *is_defined = 0;
            return 0.0;
        }
    }

    /* General case using atan2 */
    double angle = atan2(ey, ex) * (180.0 / PI);
    if (angle < 0.0)
    {
        angle += 360.0;
    }
    if (angle >= 360.0 || fabs(angle - 360.0) < 1e-9 || fabs(angle) < 1e-9)
    {
        angle = 0.0;
    }
    return angle;
}

/*
 * Evaluates electric field contributions at (ref_x, ref_y)
 * and displays individual breakdown, charge resultant,
 * and (if active) uniform background field and net resultant.
 */
static void evaluate_and_display_fields(int n, const double q[], const double x[], const double y[],
                                        double ref_x, double ref_y,
                                        int bg_active, double bg_ex, double bg_ey, double bg_mag, double bg_angle)
{
    const double k = get_coulomb_constant();
    double total_ex = 0.0;
    double total_ey = 0.0;
    int coincident_count = 0;
    int included_count = 0;

    printf("\n===============================================================\n");
    printf("                  ELECTRIC FIELD CALCULATOR\n");
    printf("===============================================================\n");
    printf("Reference Point: (%.4f, %.4f) m\n\n", ref_x, ref_y);
    printf("--- Individual Charge Contributions ---\n");

    for (int i = 0; i < n; i++)
    {
        double dx = ref_x - x[i];
        double dy = ref_y - y[i];
        double r = sqrt(dx * dx + dy * dy);

        /* Handle coincident charge (charge at the reference point) */
        if (r <= 1e-12)
        {
            printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n", i + 1, q[i], x[i], y[i]);
            printf("  Skipped because it coincides with reference point; field is undefined here.\n\n");
            coincident_count++;
            continue;
        }

        /* Convert charge from microcoulombs to Coulombs */
        double q_coulombs = q[i] * 1e-6;

        /* E = k * q / r^2 in direction of r vector (dx/r, dy/r) */
        double r3 = r * r * r;
        double ex_i = (k * q_coulombs * dx) / r3;
        double ey_i = (k * q_coulombs * dy) / r3;
        if (fabs(ex_i) <= 1e-12) ex_i = 0.0;
        if (fabs(ey_i) <= 1e-12) ey_i = 0.0;
        double mag_i = sqrt(ex_i * ex_i + ey_i * ey_i);

        int angle_defined = 0;
        double angle_i = calculate_angle(ex_i, ey_i, mag_i, &angle_defined);

        printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n", i + 1, q[i], x[i], y[i]);
        printf("  Ex = %+.4e N/C\n", ex_i);
        printf("  Ey = %+.4e N/C\n", ey_i);
        if (angle_defined)
        {
            printf("  Magnitude = %.4e N/C\n", mag_i);
            printf("  Angle = %.2f degrees\n\n", angle_i);
        }
        else
        {
            printf("  Magnitude = 0 N/C\n");
            printf("  Angle = Undefined\n\n");
        }

        total_ex += ex_i;
        total_ey += ey_i;
        included_count++;
    }

    if (coincident_count > 0)
    {
        printf("Total skipped coincident charge(s): %d\n\n", coincident_count);
    }

    /* Resultant field due to charges */
    printf("--- Charge-only Resultant Field ---\n");
    if (included_count == 0)
    {
        printf("  Note: All charges coincide with the reference point.\n");
        printf("  Physical electric field is undefined at this point.\n");
        printf("  Ex = 0.0000e+00 N/C\n");
        printf("  Ey = 0.0000e+00 N/C\n");
        printf("  Magnitude = 0 N/C\n");
        printf("  Angle = Undefined\n");
    }
    else
    {
        if (fabs(total_ex) <= 1e-12) total_ex = 0.0;
        if (fabs(total_ey) <= 1e-12) total_ey = 0.0;
        double total_mag = sqrt(total_ex * total_ex + total_ey * total_ey);
        int total_angle_defined = 0;
        double total_angle = calculate_angle(total_ex, total_ey, total_mag, &total_angle_defined);

        printf("  Ex = %+.4e N/C\n", total_ex);
        printf("  Ey = %+.4e N/C\n", total_ey);
        if (total_angle_defined)
        {
            printf("  Magnitude = %.4e N/C\n", total_mag);
            printf("  Angle = %.2f degrees\n", total_angle);
        }
        else
        {
            printf("  Magnitude = 0 N/C\n");
            printf("  Angle = Undefined\n");
        }
    }

    /* Uniform background field & Net resultant (if background is active) */
    if (bg_active)
    {
        printf("\n--- Uniform Background Electric Field ---\n");
        printf("  Ex = %+.4e N/C\n", bg_ex);
        printf("  Ey = %+.4e N/C\n", bg_ey);
        printf("  Magnitude = %.4e N/C\n", bg_mag);
        printf("  Angle = %.2f degrees\n", bg_angle);

        double net_ex = total_ex + bg_ex;
        double net_ey = total_ey + bg_ey;
        if (fabs(net_ex) <= 1e-12) net_ex = 0.0;
        if (fabs(net_ey) <= 1e-12) net_ey = 0.0;
        double net_mag = sqrt(net_ex * net_ex + net_ey * net_ey);
        int net_angle_defined = 0;
        double net_angle = calculate_angle(net_ex, net_ey, net_mag, &net_angle_defined);

        printf("\n--- Net Resultant Field (Charges + Background) ---\n");
        printf("  Ex = %+.4e N/C\n", net_ex);
        printf("  Ey = %+.4e N/C\n", net_ey);
        if (net_angle_defined)
        {
            printf("  Magnitude = %.4e N/C\n", net_mag);
            printf("  Angle = %.2f degrees\n", net_angle);
        }
        else
        {
            printf("  Magnitude = 0 N/C\n");
            printf("  Angle = Undefined\n");
        }
    }
    printf("===============================================================\n");
}

/*
 * Handles Option 2: Set / Change / Clear Uniform Background Electric Field
 */
static void handle_background_field(int *bg_active, double *bg_mag, double *bg_angle,
                                    double *bg_ex, double *bg_ey)
{
    double mag;
    while (1)
    {
        if (!read_double("\nEnter uniform background field magnitude (N/C, enter 0 to clear): ", &mag))
        {
            printf("\nInput ended (EOF). Operation canceled.\n");
            return;
        }
        if (mag < 0.0)
        {
            printf("Error: Background field magnitude cannot be negative. Please enter a non-negative value.\n");
            continue;
        }
        break;
    }

    if (mag == 0.0)
    {
        *bg_active = 0;
        *bg_mag = 0.0;
        *bg_angle = 0.0;
        *bg_ex = 0.0;
        *bg_ey = 0.0;
        printf("\nUniform background electric field cleared.\n");
    }
    else
    {
        double angle;
        if (!read_double("Enter background field direction (degrees counterclockwise from +X axis): ", &angle))
        {
            printf("\nInput ended (EOF). Operation canceled.\n");
            return;
        }

        /* Normalize angle to [0, 360) */
        double norm_angle = fmod(angle, 360.0);
        if (norm_angle < 0.0)
        {
            norm_angle += 360.0;
        }
        if (norm_angle >= 360.0 || fabs(norm_angle - 360.0) < 1e-9)
        {
            norm_angle = 0.0;
        }

        double rad = norm_angle * (PI / 180.0);
        *bg_active = 1;
        *bg_mag = mag;
        *bg_angle = norm_angle;
        *bg_ex = mag * cos(rad);
        *bg_ey = mag * sin(rad);
        printf("\nUniform background electric field set successfully.\n");
    }
}

/*
 * Handles Option 3: Find Shortest Distance for a Target Electric Field (Output 4)
 */
static void handle_output_4(int n, const double q[], const double x[], const double y[],
                            double ref_x, double ref_y)
{
    /* Find nonzero charges */
    int nonzero_indices[256];
    int nonzero_count = 0;

    for (int i = 0; i < n && nonzero_count < 256; i++)
    {
        if (fabs(q[i]) > 1e-15)
        {
            nonzero_indices[nonzero_count++] = i;
        }
    }

    if (nonzero_count == 0)
    {
        printf("\nOutput 4 cannot be performed: No nonzero charges are available.\n");
        return;
    }

    int selected_idx = -1;
    if (nonzero_count == 1)
    {
        selected_idx = nonzero_indices[0];
        printf("\nUsing single eligible charge: Charge %d (%+.4f microC at (%.4f, %.4f) m)\n",
               selected_idx + 1, q[selected_idx], x[selected_idx], y[selected_idx]);
    }
    else
    {
        printf("\nEligible nonzero charges for Output 4:\n");
        for (int k = 0; k < nonzero_count; k++)
        {
            int idx = nonzero_indices[k];
            printf("  Charge %d: %+.4f microC at (%.4f, %.4f) m\n",
                   idx + 1, q[idx], x[idx], y[idx]);
        }
        while (1)
        {
            int chosen_id;
            if (!read_int("Select charge by ID: ", &chosen_id))
            {
                printf("\nInput ended (EOF). Operation canceled.\n");
                return;
            }
            for (int k = 0; k < nonzero_count; k++)
            {
                if (nonzero_indices[k] + 1 == chosen_id)
                {
                    selected_idx = nonzero_indices[k];
                    break;
                }
            }
            if (selected_idx != -1)
            {
                break;
            }
            printf("Error: Invalid selection. Please choose an ID from the listed eligible charges.\n");
        }
    }

    /* Target electric field magnitude */
    double e_target;
    while (1)
    {
        if (!read_double("Enter target electric field magnitude (N/C): ", &e_target))
        {
            printf("\nInput ended (EOF). Operation canceled.\n");
            return;
        }
        if (e_target < 0.0)
        {
            printf("Error: Target field magnitude cannot be negative. Please enter a positive value.\n");
            continue;
        }
        if (e_target == 0.0)
        {
            printf("Note: For a nonzero point charge, an electric field magnitude of zero has no finite matching location (field approaches zero only as distance approaches infinity).\n");
            return;
        }
        break;
    }

    const double k = get_coulomb_constant();
    double q_coulombs = fabs(q[selected_idx]) * 1e-6;

    /* Radius where E = k * |q| / r^2 == e_target => r = sqrt(k * |q| / e_target) */
    double target_radius = sqrt((k * q_coulombs) / e_target);

    /* Distance from reference point to selected charge */
    double dx = ref_x - x[selected_idx];
    double dy = ref_y - y[selected_idx];
    double d = sqrt(dx * dx + dy * dy);

    /* Shortest distance to target field circle */
    double shortest_distance = fabs(d - target_radius);

    printf("\n--- Output 4 Result ---\n");
    printf("Selected Charge: Charge %d (%+.4f microC at (%.4f, %.4f) m)\n",
           selected_idx + 1, q[selected_idx], x[selected_idx], y[selected_idx]);
    printf("Target Field Magnitude: %.4e N/C\n", e_target);
    printf("Target Field Circle Radius: %.5f m\n", target_radius);
    printf("Distance from Reference Point to Charge: %.5f m\n", d);

    if (fabs(d - target_radius) < 1e-4)
    {
        printf("Shortest Distance to Target Field: 0.00000 m (reference point lies directly on the target circle)\n");
    }
    else if (d <= 1e-12)
    {
        printf("Shortest Distance to Target Field: %.5f m (reference point is at the charge location; distance is defined, no unique direction)\n",
               shortest_distance);
    }
    else
    {
        printf("Shortest Distance to Target Field: %.5f m\n", shortest_distance);
    }
}

/*
 * Main integration entry point for Module 2.
 * Called directly from Module 1 (propro.c) with collected charge data.
 */
void run_electric_field_calculator(int n, const double q[], const double x[], const double y[])
{
    if (n <= 0 || q == NULL || x == NULL || y == NULL)
    {
        printf("Error: Invalid charge data received from Module 1.\n");
        return;
    }

    /* Consume any newline leftover in stdin from Module 1's scanf */
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF)
    {
        /* discard */
    }

    printf("\n===============================================================\n");
    printf("           WELCOME TO MODULE 2: ELECTRIC FIELD CALCULATOR       \n");
    printf("===============================================================\n");
    printf("Successfully received %d charge(s) from Module 1.\n", n);

    /* Module 1 did not collect reference point, so ask user for initial reference point */
    double ref_x = 0.0, ref_y = 0.0;
    printf("Please enter the initial reference point coordinates:\n");
    if (!read_double("  Enter reference point X coordinate (m): ", &ref_x) ||
        !read_double("  Enter reference point Y coordinate (m): ", &ref_y))
    {
        printf("\nInput ended (EOF). Exiting Module 2.\n");
        return;
    }

    /* Uniform background field state */
    int bg_active = 0;
    double bg_mag = 0.0;
    double bg_angle = 0.0;
    double bg_ex = 0.0;
    double bg_ey = 0.0;

    /* Initial calculation and display */
    evaluate_and_display_fields(n, q, x, y, ref_x, ref_y,
                                bg_active, bg_ex, bg_ey, bg_mag, bg_angle);

    /* Interactive Menu Loop */
    while (1)
    {
        printf("\n==================== MENU ====================\n");
        printf("1. Change Reference Point\n");
        printf("2. Set / Change / Clear Background Electric Field\n");
        printf("3. Find Shortest Distance for a Target Electric Field (Output 4)\n");
        printf("4. Exit\n");
        printf("==============================================\n");

        int choice;
        if (!read_int("Enter your choice (1-4): ", &choice))
        {
            printf("\nInput ended (EOF). Exiting.\n");
            break;
        }

        if (choice == 1)
        {
            double new_x, new_y;
            if (!read_double("\nEnter new reference point X coordinate (m): ", &new_x) ||
                !read_double("Enter new reference point Y coordinate (m): ", &new_y))
            {
                printf("\nInput ended (EOF). Exiting.\n");
                break;
            }
            ref_x = new_x;
            ref_y = new_y;
            evaluate_and_display_fields(n, q, x, y, ref_x, ref_y,
                                        bg_active, bg_ex, bg_ey, bg_mag, bg_angle);
        }
        else if (choice == 2)
        {
            handle_background_field(&bg_active, &bg_mag, &bg_angle, &bg_ex, &bg_ey);
            evaluate_and_display_fields(n, q, x, y, ref_x, ref_y,
                                        bg_active, bg_ex, bg_ey, bg_mag, bg_angle);
        }
        else if (choice == 3)
        {
            handle_output_4(n, q, x, y, ref_x, ref_y);
        }
        else if (choice == 4 || choice == 0)
        {
            printf("\nExiting Electric Field Calculator. Thank you!\n");
            break;
        }
        else
        {
            printf("Error: Invalid choice. Please enter 1, 2, 3, or 4.\n");
        }
    }
}
