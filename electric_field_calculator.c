/*
 * Module 2: Electric Field Calculator
 * Course: IC103 Programming Project
 * Medium: Dry air at 20 deg C, 101.325 kPa
 *
 * Description:
 * This module calculates individual and resultant electric fields
 * due to multiple point charges provided by Module 1 (propro.c).
 * It supports reference point updates, uniform background field superposition,
 * and finding the shortest distance to a target field circle (Output 4).
 */

#include <stdio.h>
#include <math.h>

/* Physical constants for dry-air approximation */
#define EPS0 8.8541878188e-12      /* Permittivity of free space in F/m */
#define EPS_R 1.0005364            /* Relative permittivity of dry air */
#define PI 3.14159265358979323846  /* Pi constant */

/* Function Declarations */
static double calculate_coulomb_k(void);
static double calculate_angle_deg(double ex, double ey);
static int read_double(const char *prompt, double *val);
static int read_int(const char *prompt, int *val);
static void display_fields(int n, const double q[], const double x[], const double y[],
                           double ref_x, double ref_y,
                           int bg_active, double bg_ex, double bg_ey, double bg_mag, double bg_angle);
static void handle_background_field(int *bg_active, double *bg_mag, double *bg_angle,
                                    double *bg_ex, double *bg_ey);
static void handle_output_4(int n, const double q[], const double x[], const double y[],
                            double ref_x, double ref_y);

/* Public Entry Point called by Module 1 */
void run_electric_field_calculator(int n, const double q[], const double x[], const double y[]);

/*
 * Calculates Coulomb's constant k in dry air:
 * k = 1 / (4 * PI * EPS0 * EPS_R)
 */
static double calculate_coulomb_k(void)
{
    return 1.0 / (4.0 * PI * EPS0 * EPS_R);
}

/*
 * Calculates the direction angle in degrees [0, 360) counterclockwise from +X axis.
 * Returns -1.0 if the field magnitude is zero (undefined direction).
 */
static double calculate_angle_deg(double ex, double ey)
{
    double mag = hypot(ex, ey);
    if (mag <= 1e-15)
    {
        return -1.0; /* Magnitude is zero; direction is undefined */
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
    return angle;
}

/*
 * Safely reads a floating-point number from standard input.
 * Rejects input with unwanted trailing text.
 * Returns 1 on success, or 0 on EOF.
 */
static int read_double(const char *prompt, double *val)
{
    char buffer[128];
    char extra;

    while (1)
    {
        printf("%s", prompt);
        fflush(stdout);

        if (fgets(buffer, sizeof(buffer), stdin) == NULL)
        {
            return 0; /* EOF */
        }

        /* If the user entered an empty line, prompt again */
        if (sscanf(buffer, " %c", &extra) != 1)
        {
            continue;
        }

        /* Validate finite number without trailing extra text */
        if (sscanf(buffer, "%lf %c", val, &extra) == 1 && isfinite(*val))
        {
            return 1;
        }

        printf("Invalid input. Please enter a valid number.\n");
    }
}

/*
 * Safely reads an integer from standard input.
 * Rejects input with unwanted trailing text.
 * Returns 1 on success, or 0 on EOF.
 */
static int read_int(const char *prompt, int *val)
{
    char buffer[128];
    char extra;

    while (1)
    {
        printf("%s", prompt);
        fflush(stdout);

        if (fgets(buffer, sizeof(buffer), stdin) == NULL)
        {
            return 0; /* EOF */
        }

        /* If the user entered an empty line, prompt again */
        if (sscanf(buffer, " %c", &extra) != 1)
        {
            continue;
        }

        /* Validate integer without trailing extra text */
        if (sscanf(buffer, "%d %c", val, &extra) == 1)
        {
            return 1;
        }

        printf("Invalid input. Please enter a valid integer.\n");
    }
}

/*
 * Evaluates and displays individual electric field contributions,
 * resultant charge field, and optional background / net fields.
 */
static void display_fields(int n, const double q[], const double x[], const double y[],
                           double ref_x, double ref_y,
                           int bg_active, double bg_ex, double bg_ey, double bg_mag, double bg_angle)
{
    double k = calculate_coulomb_k();
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
        /* Exact coordinate equality check for coincident charge */
        if (x[i] == ref_x && y[i] == ref_y)
        {
            printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n", i + 1, q[i], x[i], y[i]);
            printf("  Skipped: Charge coincides with the reference point.\n");
            printf("  Electric field of a point charge is undefined at its own location.\n\n");
            coincident_count++;
            continue;
        }

        double dx = ref_x - x[i];
        double dy = ref_y - y[i];
        double r = hypot(dx, dy);

        /* Convert charge from microcoulombs (10^-6 C) to Coulombs */
        double q_coulombs = q[i] * 1e-6;

        /* E = k * q / r^2 in direction of r vector: Ex = k * q * dx / r^3, Ey = k * q * dy / r^3 */
        double r3 = r * r * r;
        double ex_i = (k * q_coulombs * dx) / r3;
        double ey_i = (k * q_coulombs * dy) / r3;
        double mag_i = hypot(ex_i, ey_i);
        double angle_i = calculate_angle_deg(ex_i, ey_i);

        printf("Charge %d (%+.4f microC at (%.4f, %.4f) m):\n", i + 1, q[i], x[i], y[i]);
        printf("  Ex = %+.4e N/C\n", ex_i);
        printf("  Ey = %+.4e N/C\n", ey_i);
        printf("  Magnitude = %.4e N/C\n", mag_i);
        if (angle_i >= 0.0)
        {
            printf("  Angle = %.2f degrees\n\n", angle_i);
        }
        else
        {
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

    /* Charge-only resultant field */
    printf("--- Charge-only Resultant Field ---\n");
    if (included_count == 0)
    {
        printf("  All charges coincide with the reference point.\n");
        printf("  Physical electric field is undefined at this point.\n");
    }
    else
    {
        double total_mag = hypot(total_ex, total_ey);
        double total_angle = calculate_angle_deg(total_ex, total_ey);

        printf("  Ex = %+.4e N/C\n", total_ex);
        printf("  Ey = %+.4e N/C\n", total_ey);
        printf("  Magnitude = %.4e N/C\n", total_mag);
        if (total_angle >= 0.0)
        {
            printf("  Angle = %.2f degrees\n", total_angle);
        }
        else
        {
            printf("  Angle = Undefined\n");
        }
    }

    /* Background field and Net resultant field */
    if (bg_active)
    {
        printf("\n--- Uniform Background Electric Field ---\n");
        printf("  Ex = %+.4e N/C\n", bg_ex);
        printf("  Ey = %+.4e N/C\n", bg_ey);
        printf("  Magnitude = %.4e N/C\n", bg_mag);
        printf("  Angle = %.2f degrees\n", bg_angle);

        if (included_count > 0)
        {
            double net_ex = total_ex + bg_ex;
            double net_ey = total_ey + bg_ey;
            double net_mag = hypot(net_ex, net_ey);
            double net_angle = calculate_angle_deg(net_ex, net_ey);

            printf("\n--- Net Resultant Field (Charges + Background) ---\n");
            printf("  Ex = %+.4e N/C\n", net_ex);
            printf("  Ey = %+.4e N/C\n", net_ey);
            printf("  Magnitude = %.4e N/C\n", net_mag);
            if (net_angle >= 0.0)
            {
                printf("  Angle = %.2f degrees\n", net_angle);
            }
            else
            {
                printf("  Angle = Undefined\n");
            }
        }
        else
        {
            printf("\n--- Net Resultant Field (Charges + Background) ---\n");
            printf("  Cannot calculate net field because charge contribution is undefined.\n");
        }
    }
    printf("===============================================================\n");
}

/*
 * Option 2: Set, change, or clear uniform background electric field.
 */
static void handle_background_field(int *bg_active, double *bg_mag, double *bg_angle,
                                    double *bg_ex, double *bg_ey)
{
    double mag;
    while (1)
    {
        if (!read_double("\nEnter uniform background field magnitude (N/C, enter 0 to clear): ", &mag))
        {
            return;
        }
        if (mag < 0.0)
        {
            printf("Error: Background field magnitude cannot be negative.\n");
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
        printf("Uniform background electric field cleared.\n");
    }
    else
    {
        double angle;
        if (!read_double("Enter background field direction (degrees counterclockwise from +X axis): ", &angle))
        {
            return;
        }

        /* Normalize angle to [0, 360) */
        double norm_angle = fmod(angle, 360.0);
        if (norm_angle < 0.0)
        {
            norm_angle += 360.0;
        }
        if (norm_angle >= 360.0)
        {
            norm_angle = 0.0;
        }

        double rad = norm_angle * (PI / 180.0);
        *bg_active = 1;
        *bg_mag = mag;
        *bg_angle = norm_angle;
        *bg_ex = mag * cos(rad);
        *bg_ey = mag * sin(rad);
        printf("Uniform background electric field set successfully.\n");
    }
}

/*
 * Option 3: Find shortest distance from current reference point to
 * where a selected charge produces target electric field magnitude (Output 4).
 */
static void handle_output_4(int n, const double q[], const double x[], const double y[],
                            double ref_x, double ref_y)
{
    /* Count nonzero charges */
    int nonzero_count = 0;
    int single_nonzero_idx = -1;

    for (int i = 0; i < n; i++)
    {
        if (q[i] != 0.0)
        {
            nonzero_count++;
            single_nonzero_idx = i;
        }
    }

    if (nonzero_count == 0)
    {
        printf("\nOutput 4 is unavailable: No nonzero charges are present.\n");
        return;
    }

    int selected_idx = -1;
    if (nonzero_count == 1)
    {
        selected_idx = single_nonzero_idx;
        printf("\nUsing the only nonzero charge: Charge %d (%+.4f microC at (%.4f, %.4f) m)\n",
               selected_idx + 1, q[selected_idx], x[selected_idx], y[selected_idx]);
    }
    else
    {
        printf("\nEligible nonzero charges:\n");
        for (int i = 0; i < n; i++)
        {
            if (q[i] != 0.0)
            {
                printf("  Charge %d: %+.4f microC at (%.4f, %.4f) m\n",
                       i + 1, q[i], x[i], y[i]);
            }
        }

        while (1)
        {
            int chosen_id;
            if (!read_int("Select charge by ID: ", &chosen_id))
            {
                return;
            }
            if (chosen_id >= 1 && chosen_id <= n && q[chosen_id - 1] != 0.0)
            {
                selected_idx = chosen_id - 1;
                break;
            }
            printf("Invalid selection. Please choose an ID from the listed eligible charges.\n");
        }
    }

    double e_target;
    while (1)
    {
        if (!read_double("Enter target electric field magnitude (N/C): ", &e_target))
        {
            return;
        }
        if (e_target < 0.0)
        {
            printf("Error: Target field magnitude cannot be negative. Please enter a positive value.\n");
            continue;
        }
        if (e_target == 0.0)
        {
            printf("A zero target electric field has no finite matching location for a point charge.\n");
            return;
        }
        break;
    }

    double k = calculate_coulomb_k();
    double q_coulombs = fabs(q[selected_idx]) * 1e-6;
    double target_radius = sqrt((k * q_coulombs) / e_target);

    double dx = ref_x - x[selected_idx];
    double dy = ref_y - y[selected_idx];
    double d = hypot(dx, dy);

    double shortest_distance = fabs(d - target_radius);

    printf("\n--- Output 4 Result ---\n");
    printf("Selected Charge: Charge %d (%+.4f microC at (%.4f, %.4f) m)\n",
           selected_idx + 1, q[selected_idx], x[selected_idx], y[selected_idx]);
    printf("Target Field Magnitude: %.4e N/C\n", e_target);
    printf("Target Field Circle Radius: %.4f m\n", target_radius);
    printf("Distance from Reference Point to Charge: %.4f m\n", d);
    printf("Shortest Distance to Target Field: %.4f m\n", shortest_distance);
}

/*
 * Public entry point for Module 2: Electric Field Calculator.
 * Receives charge data directly from Module 1.
 */
void run_electric_field_calculator(int n, const double q[], const double x[], const double y[])
{
    if (n <= 0 || q == NULL || x == NULL || y == NULL)
    {
        printf("Error: Invalid charge data received from Module 1.\n");
        return;
    }

    /* Clear any leftover newline from Module 1's input */
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF)
    {
        /* discard */
    }

    printf("\n===============================================================\n");
    printf("           WELCOME TO MODULE 2: ELECTRIC FIELD CALCULATOR       \n");
    printf("===============================================================\n");
    printf("Successfully received %d charge(s) from Module 1.\n", n);

    double ref_x = 0.0;
    double ref_y = 0.0;
    printf("\nPlease enter the initial reference point coordinates:\n");
    if (!read_double("  Enter reference point X coordinate (m): ", &ref_x) ||
        !read_double("  Enter reference point Y coordinate (m): ", &ref_y))
    {
        return;
    }

    /* Background field state */
    int bg_active = 0;
    double bg_mag = 0.0;
    double bg_angle = 0.0;
    double bg_ex = 0.0;
    double bg_ey = 0.0;

    /* Initial calculation and display */
    display_fields(n, q, x, y, ref_x, ref_y,
                   bg_active, bg_ex, bg_ey, bg_mag, bg_angle);

    /* Follow-up menu loop */
    while (1)
    {
        printf("\n==================== MENU ====================\n");
        printf("1. Change reference point\n");
        printf("2. Set / change / clear uniform background electric field\n");
        printf("3. Find shortest distance for a target electric field (Output 4)\n");
        printf("0. Exit\n");
        printf("==============================================\n");

        int choice;
        if (!read_int("Enter your choice: ", &choice))
        {
            break;
        }

        if (choice == 1)
        {
            double new_x, new_y;
            printf("\nEnter new reference point coordinates:\n");
            if (!read_double("  Enter reference point X coordinate (m): ", &new_x) ||
                !read_double("  Enter reference point Y coordinate (m): ", &new_y))
            {
                break;
            }
            ref_x = new_x;
            ref_y = new_y;
            display_fields(n, q, x, y, ref_x, ref_y,
                           bg_active, bg_ex, bg_ey, bg_mag, bg_angle);
        }
        else if (choice == 2)
        {
            handle_background_field(&bg_active, &bg_mag, &bg_angle, &bg_ex, &bg_ey);
            display_fields(n, q, x, y, ref_x, ref_y,
                           bg_active, bg_ex, bg_ey, bg_mag, bg_angle);
        }
        else if (choice == 3)
        {
            handle_output_4(n, q, x, y, ref_x, ref_y);
        }
        else if (choice == 0)
        {
            printf("\nExiting Electric Field Calculator. Thank you!\n");
            break;
        }
        else
        {
            printf("Invalid choice. Please enter 0, 1, 2, or 3.\n");
        }
    }
}
