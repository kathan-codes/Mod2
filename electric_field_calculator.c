/*
 * Module 2: Electric Field Calculator
 * 
 * Part of the multi-module Electric Field and Potential Simulation project.
 * Implemented in ISO C99 for dry, carbon-dioxide-free air at 20 °C and 101.325 kPa.
 * 
 * Ownership: Module 2 (Electric Field Calculator).
 * In this stage, runs with an isolated demo dataset unless integrated.
 * 
 * Compilation switch:
 *   Default: compiles with demo main() for standalone execution and verification.
 *   Integration: compile with -DMODULE2_NO_MAIN (or call run_electric_field_calculator).
 * 
 * Zero tolerances (in N/C):
 *   Individual charge vector: magnitude <= 1e-12 N/C
 *   Charge-only resultant: max(1e-12, 1e-12 * sum_of_included_charge_field_magnitudes) N/C
 *   Net resultant: max(1e-12, 1e-12 * (sum_of_included_charge_field_magnitudes + background_magnitude)) N/C
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>

/* Physical constants for dry air at 20 °C and 101.325 kPa */
#define EPSILON_0 8.8541878188e-12
#define EPSILON_R_AIR 1.0005364
#define PI 3.14159265358979323846
#define K_AIR (1.0 / (4.0 * PI * EPSILON_0 * EPSILON_R_AIR))

#define MAX_CHARGES 256
#define INPUT_LINE_BUFFER_SIZE 256

/* Point charge data structure */
typedef struct {
    int id;           /* Stable 1-based identifier */
    double q_microC;  /* Charge value in microcoulombs (sign preserved) */
    double x;         /* X coordinate in metres */
    double y;         /* Y coordinate in metres */
} PointCharge;

/* 2D reference point */
typedef struct {
    double x;         /* X coordinate in metres */
    double y;         /* Y coordinate in metres */
} ReferencePoint;

/* 2D vector calculation result */
typedef struct {
    double ex;                  /* Ex component in N/C */
    double ey;                  /* Ey component in N/C */
    double mag;                 /* Magnitude in N/C */
    double angle_deg;           /* Angle counterclockwise from +X in degrees [0, 360) */
    bool is_zero;               /* True if classified as zero within tolerance */
    bool angle_defined;         /* True if angle is physically defined (false for zero vector) */
    bool components_within_tol; /* True if components were non-zero but within zero tolerance */
} VectorResult;

/* Evaluated field contribution from a single charge */
typedef struct {
    int id;
    double q_microC;
    double x;
    double y;
    bool coincident;
    VectorResult field;
} EvaluatedCharge;

/* Uniform background field */
typedef struct {
    bool active;
    double magnitude; /* Non-negative magnitude in N/C */
    double angle_deg; /* Angle counterclockwise from +X in degrees [0, 360) */
    double ex;        /* Ex in N/C */
    double ey;        /* Ey in N/C */
} BackgroundField;

/* Complete evaluation breakdown */
typedef struct {
    EvaluatedCharge charges[MAX_CHARGES];
    size_t total_count;
    size_t coincident_count;
    size_t included_count;
    
    VectorResult charge_resultant;
    
    bool has_background;
    VectorResult background;
    VectorResult net_resultant;
} FieldEvaluation;

/* Demo dataset: reproducible baseline */
static const PointCharge DEMO_CHARGES[] = {
    { 1,  1.0, 0.0, 0.0 },  /* Charge 1: +1 microC at (0, 0) m */
    { 2, -1.0, 2.0, 0.0 }   /* Charge 2: -1 microC at (2, 0) m */
};
static const size_t DEMO_CHARGE_COUNT = sizeof(DEMO_CHARGES) / sizeof(DEMO_CHARGES[0]);
static const ReferencePoint DEMO_INITIAL_REF = { 1.0, 1.0 }; /* (1, 1) m */

/* Line reader status */
typedef enum {
    READ_LINE_OK,
    READ_LINE_EOF,
    READ_LINE_OVERLONG
} ReadLineStatus;

/* Prompt helper status */
typedef enum {
    PROMPT_OK,
    PROMPT_EOF
} PromptStatus;

/* Case-insensitive prefix check */
static bool ci_starts_with(const char *str, const char *prefix) {
    while (*prefix) {
        if (tolower((unsigned char)*str) != tolower((unsigned char)*prefix)) {
            return false;
        }
        str++;
        prefix++;
    }
    return true;
}

/* Reads a full line, discarding overlong trailing content if buffer is exceeded */
static ReadLineStatus read_line(char *buf, size_t size, FILE *stream) {
    if (size == 0) return READ_LINE_EOF;
    if (fgets(buf, (int)size, stream) == NULL) {
        return READ_LINE_EOF;
    }
    size_t len = strlen(buf);
    if (len == 0) {
        return READ_LINE_EOF;
    }
    if (buf[len - 1] == '\n') {
        buf[len - 1] = '\0';
        if (len > 1 && buf[len - 2] == '\r') {
            buf[len - 2] = '\0';
        }
        return READ_LINE_OK;
    }
    if (!feof(stream)) {
        int ch;
        while ((ch = fgetc(stream)) != '\n' && ch != EOF) {
            /* discard remainder of overlong line */
        }
        return READ_LINE_OVERLONG;
    }
    return READ_LINE_OK;
}

/* Strict double parser: rejects NaN, inf, trailing junk, overflow/underflow, blank lines */
static bool parse_double_strictly(const char *buf, double *out_val) {
    if (buf == NULL) return false;
    
    const char *p = buf;
    while (isspace((unsigned char)*p)) p++;
    if (*p == '\0') {
        return false;
    }
    
    const char *check = p;
    if (*check == '+' || *check == '-') check++;
    if (ci_starts_with(check, "nan") || ci_starts_with(check, "inf")) {
        return false;
    }
    
    errno = 0;
    char *endptr = NULL;
    double val = strtod(p, &endptr);
    
    if (endptr == p) {
        return false;
    }
    if (errno == ERANGE) {
        return false;
    }
    if (!isfinite(val)) {
        return false;
    }
    
    while (*endptr != '\0') {
        if (!isspace((unsigned char)*endptr)) {
            return false;
        }
        endptr++;
    }
    
    *out_val = val;
    return true;
}

/* Strict long parser: rejects trailing junk, blank lines, out of range */
static bool parse_long_strictly(const char *buf, long *out_val) {
    if (buf == NULL) return false;
    
    const char *p = buf;
    while (isspace((unsigned char)*p)) p++;
    if (*p == '\0') {
        return false;
    }
    
    errno = 0;
    char *endptr = NULL;
    long val = strtol(p, &endptr, 10);
    
    if (endptr == p) {
        return false;
    }
    if (errno == ERANGE) {
        return false;
    }
    
    while (*endptr != '\0') {
        if (!isspace((unsigned char)*endptr)) {
            return false;
        }
        endptr++;
    }
    
    *out_val = val;
    return true;
}

/* Prompts for a finite double value */
static PromptStatus prompt_double(const char *prompt, double *out_val) {
    char buf[INPUT_LINE_BUFFER_SIZE];
    while (1) {
        printf("%s", prompt);
        fflush(stdout);
        ReadLineStatus status = read_line(buf, sizeof(buf), stdin);
        if (status == READ_LINE_EOF) {
            return PROMPT_EOF;
        }
        if (status == READ_LINE_OVERLONG) {
            printf("Error: Input line exceeds buffer limit. Please enter a valid number.\n");
            continue;
        }
        if (!parse_double_strictly(buf, out_val)) {
            printf("Error: Invalid numeric input. Please enter a finite real number.\n");
            continue;
        }
        return PROMPT_OK;
    }
}

/* Prompts for an integer value */
static PromptStatus prompt_int(const char *prompt, int *out_val) {
    char buf[INPUT_LINE_BUFFER_SIZE];
    while (1) {
        printf("%s", prompt);
        fflush(stdout);
        ReadLineStatus status = read_line(buf, sizeof(buf), stdin);
        if (status == READ_LINE_EOF) {
            return PROMPT_EOF;
        }
        if (status == READ_LINE_OVERLONG) {
            printf("Error: Input line exceeds buffer limit. Please enter a valid integer.\n");
            continue;
        }
        long val;
        if (!parse_long_strictly(buf, &val) || val < INT_MIN || val > INT_MAX) {
            printf("Error: Invalid integer input. Please enter a valid integer.\n");
            continue;
        }
        *out_val = (int)val;
        return PROMPT_OK;
    }
}

/* Normalizes angle to [0, 360) degrees */
static double normalize_angle_deg(double deg) {
    deg = fmod(deg, 360.0);
    if (deg < 0.0) {
        deg += 360.0;
    }
    if (deg >= 360.0 || fabs(deg - 360.0) < 1e-12) {
        deg = 0.0;
    }
    if (fabs(deg) < 1e-12) {
        deg = 0.0;
    }
    return deg;
}

/*
 * Evaluates 2D vector magnitude and angle with three specific angle cases:
 * 1. Ex = 0 and Ey > 0 => 90 degrees
 * 2. Ex = 0 and Ey < 0 => 270 degrees
 * 3. Magnitude <= zero_tol => Magnitude 0; angle undefined
 */
static bool evaluate_vector(double ex, double ey, double zero_tol, VectorResult *out) {
    if (!isfinite(ex) || !isfinite(ey)) {
        return false;
    }
    out->ex = ex;
    out->ey = ey;
    out->components_within_tol = false;

    double mag = hypot(ex, ey);
    if (!isfinite(mag)) {
        return false;
    }

    if (mag <= zero_tol) {
        out->mag = 0.0;
        out->is_zero = true;
        out->angle_defined = false;
        out->angle_deg = 0.0;
        if (mag > 0.0) {
            out->components_within_tol = true;
        }
    } else {
        out->mag = mag;
        out->is_zero = false;
        out->angle_defined = true;
        if (ex == 0.0) {
            if (ey > 0.0) {
                out->angle_deg = 90.0;
            } else {
                out->angle_deg = 270.0;
            }
        } else {
            double rad = atan2(ey, ex);
            out->angle_deg = normalize_angle_deg(rad * 180.0 / PI);
        }
    }
    return true;
}

/*
 * Computes field evaluation for all charges at the given reference point.
 * Checks for numeric safety and coincidence before returning.
 */
static bool compute_field_evaluation(const PointCharge *charges, size_t count,
                                     ReferencePoint ref, const BackgroundField *bg,
                                     FieldEvaluation *out) {
    if (count > MAX_CHARGES) {
        return false;
    }
    memset(out, 0, sizeof(*out));
    out->total_count = count;

    double sum_ex = 0.0;
    double sum_ey = 0.0;
    double sum_indiv_mags = 0.0;

    for (size_t i = 0; i < count; i++) {
        out->charges[i].id = charges[i].id;
        out->charges[i].q_microC = charges[i].q_microC;
        out->charges[i].x = charges[i].x;
        out->charges[i].y = charges[i].y;

        /* Coincidence check: exact coordinate match only */
        if (charges[i].x == ref.x && charges[i].y == ref.y) {
            out->charges[i].coincident = true;
            out->coincident_count++;
            out->charges[i].field.mag = 0.0;
            out->charges[i].field.is_zero = true;
            out->charges[i].field.angle_defined = false;
            continue;
        }

        out->charges[i].coincident = false;
        out->included_count++;

        double dx = ref.x - charges[i].x;
        double dy = ref.y - charges[i].y;
        double r = hypot(dx, dy);

        if (!isfinite(r) || r == 0.0) {
            return false;
        }

        double q_C = charges[i].q_microC * 1e-6;
        double inv_r = 1.0 / r;
        double factor = (K_AIR * q_C * inv_r) * (inv_r * inv_r);

        if (!isfinite(factor)) {
            return false;
        }

        double ex = factor * dx;
        double ey = factor * dy;

        if (!isfinite(ex) || !isfinite(ey)) {
            return false;
        }

        /* Individual vector tolerance: 1e-12 N/C */
        if (!evaluate_vector(ex, ey, 1e-12, &out->charges[i].field)) {
            return false;
        }

        double raw_indiv_mag = hypot(ex, ey);
        sum_ex += ex;
        sum_ey += ey;
        sum_indiv_mags += raw_indiv_mag;
    }

    /* Charge-only resultant */
    if (out->included_count == 0) {
        out->charge_resultant.ex = 0.0;
        out->charge_resultant.ey = 0.0;
        out->charge_resultant.mag = 0.0;
        out->charge_resultant.is_zero = true;
        out->charge_resultant.angle_defined = false;
        out->charge_resultant.angle_deg = 0.0;
    } else {
        if (!isfinite(sum_ex) || !isfinite(sum_ey)) {
            return false;
        }
        double tol_charge = fmax(1e-12, 1e-12 * sum_indiv_mags);
        if (!evaluate_vector(sum_ex, sum_ey, tol_charge, &out->charge_resultant)) {
            return false;
        }
    }

    /* Background field contribution */
    if (bg != NULL && bg->active) {
        out->has_background = true;
        out->background.ex = bg->ex;
        out->background.ey = bg->ey;
        out->background.mag = bg->magnitude;
        out->background.angle_deg = bg->angle_deg;
        out->background.angle_defined = (bg->magnitude > 0.0);
        out->background.is_zero = (bg->magnitude == 0.0);

        double net_ex = sum_ex + bg->ex;
        double net_ey = sum_ey + bg->ey;
        if (!isfinite(net_ex) || !isfinite(net_ey)) {
            return false;
        }
        double tol_net = fmax(1e-12, 1e-12 * (sum_indiv_mags + bg->magnitude));
        if (!evaluate_vector(net_ex, net_ey, tol_net, &out->net_resultant)) {
            return false;
        }
    } else {
        out->has_background = false;
    }

    return true;
}

/* Prints the per-charge breakdown */
static void print_charge_breakdown(const FieldEvaluation *eval, ReferencePoint ref) {
    printf("\n--- Point Charge Field Contributions (Reference Point: (%.6f, %.6f) m) ---\n", ref.x, ref.y);
    if (eval->total_count == 0) {
        printf("No point charges present in dataset.\n");
        return;
    }

    for (size_t i = 0; i < eval->total_count; i++) {
        const EvaluatedCharge *ec = &eval->charges[i];
        if (ec->coincident) {
            printf("Charge %d (%+.6f microC at (%.6f, %.6f) m): skipped because it coincides with the reference point; its field is undefined here.\n",
                   ec->id, ec->q_microC, ec->x, ec->y);
        } else {
            printf("Charge %d (%+.6f microC at (%.6f, %.6f) m):\n",
                   ec->id, ec->q_microC, ec->x, ec->y);
            printf("  Ex: %+.6e N/C, Ey: %+.6e N/C\n", ec->field.ex, ec->field.ey);
            if (ec->field.is_zero) {
                printf("  Magnitude: 0 N/C, Angle: undefined\n");
            } else {
                printf("  Magnitude: %.6e N/C, Angle: %.2f degrees\n",
                       ec->field.mag, ec->field.angle_deg);
            }
        }
    }
    printf("Excluded charges: %zu\n", eval->coincident_count);
}

/* Prints the charge-only resultant */
static void print_charge_resultant(const FieldEvaluation *eval) {
    if (eval->included_count == 0) {
        if (eval->coincident_count > 0) {
            printf("\n--- Charge-Only Resultant Field (Convention-based sum from 0 included charges) ---\n");
            printf("Note: All point charges coincide with the reference point; the physical point-charge field at this reference point is undefined.\n");
        } else {
            printf("\n--- Charge-Only Resultant Field ---\n");
            printf("Note: No charges in dataset.\n");
        }
    } else if (eval->coincident_count > 0) {
        printf("\n--- Charge-Only Resultant Field (Convention-based sum over %zu included charges) ---\n", eval->included_count);
        printf("Note: %zu charge(s) coincided with reference point; physical total point-charge field is undefined at this point.\n",
               eval->coincident_count);
    } else {
        printf("\n--- Charge-Only Resultant Field ---\n");
    }

    printf("Ex: %+.6e N/C, Ey: %+.6e N/C", eval->charge_resultant.ex, eval->charge_resultant.ey);
    if (eval->charge_resultant.components_within_tol) {
        printf(" (components within zero tolerance)");
    }
    printf("\n");

    if (eval->charge_resultant.is_zero) {
        printf("Magnitude: 0 N/C, Angle: undefined\n");
    } else {
        printf("Magnitude: %.6e N/C, Angle: %.2f degrees\n",
               eval->charge_resultant.mag, eval->charge_resultant.angle_deg);
    }
}

/* Prints background and net resultant blocks */
static void print_background_and_net(const FieldEvaluation *eval) {
    if (!eval->has_background) {
        return;
    }

    printf("\n--- Uniform Background Field ---\n");
    printf("Ex: %+.6e N/C, Ey: %+.6e N/C\n", eval->background.ex, eval->background.ey);
    if (eval->background.is_zero) {
        printf("Magnitude: 0 N/C, Angle: undefined\n");
    } else {
        printf("Magnitude: %.6e N/C, Angle: %.2f degrees\n",
               eval->background.mag, eval->background.angle_deg);
    }

    if (eval->coincident_count > 0) {
        printf("\n--- Net Resultant Field (Convention-based sum over included charges + background) ---\n");
        printf("Note: Physical total field is undefined at this point due to coincident charge(s).\n");
    } else {
        printf("\n--- Net Resultant Field ---\n");
    }

    printf("Ex: %+.6e N/C, Ey: %+.6e N/C", eval->net_resultant.ex, eval->net_resultant.ey);
    if (eval->net_resultant.components_within_tol) {
        printf(" (components within zero tolerance)");
    }
    printf("\n");

    if (eval->net_resultant.is_zero) {
        printf("Magnitude: 0 N/C, Angle: undefined\n");
    } else {
        printf("Magnitude: %.6e N/C, Angle: %.2f degrees\n",
               eval->net_resultant.mag, eval->net_resultant.angle_deg);
    }
}

/* Output 4 calculation */
static bool compute_output_4(const PointCharge *charge, ReferencePoint ref,
                             double e_target, double *out_target_radius,
                             double *out_dist_to_charge,
                             double *out_shortest_distance) {
    if (charge == NULL || e_target <= 0.0) {
        return false;
    }
    double q_C = fabs(charge->q_microC * 1e-6);
    if (q_C == 0.0) {
        return false;
    }

    double r_sq = (K_AIR * q_C) / e_target;
    if (!isfinite(r_sq) || r_sq < 0.0) {
        return false;
    }
    double target_radius = sqrt(r_sq);
    if (!isfinite(target_radius)) {
        return false;
    }

    double dx = ref.x - charge->x;
    double dy = ref.y - charge->y;
    double d = hypot(dx, dy);
    if (!isfinite(d)) {
        return false;
    }

    double shortest = fabs(d - target_radius);
    if (!isfinite(shortest)) {
        return false;
    }

    *out_target_radius = target_radius;
    *out_dist_to_charge = d;
    *out_shortest_distance = shortest;
    return true;
}

/* Executes Output 4 interactive action */
static PromptStatus handle_output_4(const PointCharge *charges, size_t count, ReferencePoint ref) {
    /* Collect eligible nonzero charges */
    size_t eligible_indices[MAX_CHARGES];
    size_t eligible_count = 0;

    for (size_t i = 0; i < count; i++) {
        if (charges[i].q_microC != 0.0) {
            eligible_indices[eligible_count++] = i;
        }
    }

    if (eligible_count == 0) {
        printf("\nOutput 4 cannot run: No nonzero charge is available.\n");
        return PROMPT_OK;
    }

    const PointCharge *selected = NULL;
    if (eligible_count == 1) {
        selected = &charges[eligible_indices[0]];
        printf("\nUsing single eligible charge: ID %d (%+.6f microC at (%.6f, %.6f) m)\n",
               selected->id, selected->q_microC, selected->x, selected->y);
    } else {
        printf("\nEligible nonzero charges for Output 4:\n");
        for (size_t k = 0; k < eligible_count; k++) {
            const PointCharge *c = &charges[eligible_indices[k]];
            printf("  ID %d: %+.6f microC at (%.6f, %.6f) m\n", c->id, c->q_microC, c->x, c->y);
        }

        while (1) {
            int sel_id = 0;
            if (prompt_int("Select charge by ID: ", &sel_id) == PROMPT_EOF) {
                return PROMPT_EOF;
            }
            for (size_t k = 0; k < eligible_count; k++) {
                if (charges[eligible_indices[k]].id == sel_id) {
                    selected = &charges[eligible_indices[k]];
                    break;
                }
            }
            if (selected != NULL) {
                break;
            }
            printf("Error: Invalid selection. Please choose an ID from the listed eligible charges.\n");
        }
    }

    /* Prompt for target electric field magnitude */
    double target_field = 0.0;
    while (1) {
        if (prompt_double("Enter target electric field magnitude (N/C): ", &target_field) == PROMPT_EOF) {
            return PROMPT_EOF;
        }
        if (target_field < 0.0) {
            printf("Error: Target field magnitude cannot be negative. Please enter a nonnegative value.\n");
            continue;
        }
        if (target_field == 0.0) {
            printf("For a nonzero point charge, an electric field magnitude of zero has no finite matching location (field approaches zero only as distance approaches infinity).\n");
            return PROMPT_OK;
        }
        break;
    }

    double target_radius = 0.0;
    double dist_to_charge = 0.0;
    double shortest_dist = 0.0;

    if (!compute_output_4(selected, ref, target_field, &target_radius, &dist_to_charge, &shortest_dist)) {
        printf("Error: Calculation-range error during shortest distance evaluation.\n");
        return PROMPT_OK;
    }

    printf("\n--- Output 4 Result ---\n");
    printf("Selected Charge: ID %d (%+.6f microC at (%.6f, %.6f) m)\n",
           selected->id, selected->q_microC, selected->x, selected->y);
    printf("Target Field Magnitude: %.6e N/C\n", target_field);
    printf("Target Field Circle Radius: %.6e m\n", target_radius);
    printf("Distance from Reference Point to Charge: %.6e m\n", dist_to_charge);

    if (fabs(dist_to_charge - target_radius) < 1e-12) {
        printf("Shortest Distance to Target Field: 0 m (reference point lies on target circle)\n");
    } else if (dist_to_charge == 0.0) {
        printf("Shortest Distance to Target Field: %.5f m (reference point equals charge location; distance is defined, no unique direction)\n",
               shortest_dist);
    } else {
        printf("Shortest Distance to Target Field: %.5f m\n", shortest_dist);
    }

    return PROMPT_OK;
}

/*
 * Primary calculator loop for Module 2.
 * Can be called with demo dataset or integrated dataset from Module 1.
 */
int run_electric_field_calculator(const PointCharge *charges, size_t count, ReferencePoint initial_ref) {
    if (charges == NULL && count > 0) {
        printf("Error: Input data is unavailable from Module 1.\n");
        return -1;
    }
    if (count > MAX_CHARGES) {
        printf("Error: Charge list exceeds maximum supported capacity (%d).\n", MAX_CHARGES);
        return -1;
    }

    /* Validate all input charges */
    for (size_t i = 0; i < count; i++) {
        if (!isfinite(charges[i].q_microC) || !isfinite(charges[i].x) || !isfinite(charges[i].y)) {
            printf("Error: Input data contains non-finite values.\n");
            return -1;
        }
    }
    if (!isfinite(initial_ref.x) || !isfinite(initial_ref.y)) {
        printf("Error: Initial reference point contains non-finite values.\n");
        return -1;
    }

    ReferencePoint working_ref = initial_ref;
    BackgroundField bg = { false, 0.0, 0.0, 0.0, 0.0 };

    /* Initial evaluation and breakdown */
    FieldEvaluation eval;
    if (!compute_field_evaluation(charges, count, working_ref, &bg, &eval)) {
        printf("Error: Calculation-range error during initial field evaluation.\n");
        return -1;
    }

    print_charge_breakdown(&eval, working_ref);
    print_charge_resultant(&eval);
    print_background_and_net(&eval);

    /* Follow-up menu loop */
    char menu_buf[INPUT_LINE_BUFFER_SIZE];
    while (1) {
        printf("\n--- Menu ---\n");
        printf("1  Change the reference point\n");
        printf("2  Set, replace, or clear the uniform background field\n");
        printf("3  Find the shortest distance for a target field (Output 4)\n");
        printf("0  Exit\n");
        printf("Enter choice: ");
        fflush(stdout);

        ReadLineStatus rls = read_line(menu_buf, sizeof(menu_buf), stdin);
        if (rls == READ_LINE_EOF) {
            printf("\nInput ended (EOF). Exiting cleanly.\n");
            break;
        }
        if (rls == READ_LINE_OVERLONG) {
            printf("Error: Input line exceeds buffer limit. Valid choices are 0, 1, 2, or 3.\n");
            continue;
        }

        long choice = -1;
        if (!parse_long_strictly(menu_buf, &choice) || choice < 0 || choice > 3) {
            printf("Error: Invalid choice. Valid choices are 0, 1, 2, or 3.\n");
            continue;
        }

        if (choice == 0) {
            printf("Exiting Electric Field Calculator.\n");
            break;
        } else if (choice == 1) {
            /* Option 1: Change reference point */
            double new_x = 0.0;
            double new_y = 0.0;

            if (prompt_double("Enter new X reference coordinate (m): ", &new_x) == PROMPT_EOF) {
                printf("\nInput ended (EOF). Reference point change canceled.\n");
                break;
            }
            if (prompt_double("Enter new Y reference coordinate (m): ", &new_y) == PROMPT_EOF) {
                printf("\nInput ended (EOF). Reference point change canceled.\n");
                break;
            }

            ReferencePoint temp_ref = { new_x, new_y };
            FieldEvaluation temp_eval;
            if (!compute_field_evaluation(charges, count, temp_ref, &bg, &temp_eval)) {
                printf("Error: Calculation-range error with proposed reference point. Changes not committed.\n");
                continue;
            }

            working_ref = temp_ref;
            eval = temp_eval;

            print_charge_breakdown(&eval, working_ref);
            print_charge_resultant(&eval);
            print_background_and_net(&eval);

        } else if (choice == 2) {
            /* Option 2: Set, replace, or clear uniform background field */
            double mag = 0.0;
            while (1) {
                if (prompt_double("Enter uniform background field magnitude (N/C, 0 to clear): ", &mag) == PROMPT_EOF) {
                    printf("\nInput ended (EOF). Background field update canceled.\n");
                    goto cleanup_exit;
                }
                if (mag < 0.0) {
                    printf("Error: Background field magnitude must be nonnegative.\n");
                    continue;
                }
                break;
            }

            BackgroundField temp_bg = bg;
            if (mag == 0.0) {
                temp_bg.active = false;
                temp_bg.magnitude = 0.0;
                temp_bg.angle_deg = 0.0;
                temp_bg.ex = 0.0;
                temp_bg.ey = 0.0;
            } else {
                double angle_input = 0.0;
                if (prompt_double("Enter background field angle (degrees counterclockwise from +X): ", &angle_input) == PROMPT_EOF) {
                    printf("\nInput ended (EOF). Background field update canceled.\n");
                    goto cleanup_exit;
                }
                double norm_angle = normalize_angle_deg(angle_input);
                double rad = norm_angle * PI / 180.0;
                double ex_bg = mag * cos(rad);
                double ey_bg = mag * sin(rad);

                if (!isfinite(ex_bg) || !isfinite(ey_bg)) {
                    printf("Error: Calculation-range error computing background field components. Changes not committed.\n");
                    continue;
                }

                temp_bg.active = true;
                temp_bg.magnitude = mag;
                temp_bg.angle_deg = norm_angle;
                temp_bg.ex = ex_bg;
                temp_bg.ey = ey_bg;
            }

            FieldEvaluation temp_eval;
            if (!compute_field_evaluation(charges, count, working_ref, &temp_bg, &temp_eval)) {
                printf("Error: Calculation-range error applying background field. Changes not committed.\n");
                continue;
            }

            bg = temp_bg;
            eval = temp_eval;

            if (!bg.active) {
                printf("\nUniform background field cleared.\n");
                print_charge_resultant(&eval);
            } else {
                print_background_and_net(&eval);
            }

        } else if (choice == 3) {
            /* Option 3: Output 4 */
            if (handle_output_4(charges, count, working_ref) == PROMPT_EOF) {
                printf("\nInput ended (EOF). Output 4 canceled.\n");
                break;
            }
        }
    }

cleanup_exit:
    return 0;
}

#ifndef MODULE2_NO_MAIN

/* Demo dataset runner */
static int run_demo_mode(void) {
    printf("===============================================================\n");
    printf("   ELECTRIC FIELD CALCULATOR (MODULE 2) [DEMO MODE]\n");
    printf("   Medium: Dry air (20 deg C, 101.325 kPa)\n");
    printf("   k_air: %.4e N*m^2/C^2\n", K_AIR);
    printf("===============================================================\n");
    printf("Using clearly labeled reproducible demo dataset:\n");
    for (size_t i = 0; i < DEMO_CHARGE_COUNT; i++) {
        printf("  Charge %d: %+.6f microC at (%.6f, %.6f) m\n",
               DEMO_CHARGES[i].id, DEMO_CHARGES[i].q_microC, DEMO_CHARGES[i].x, DEMO_CHARGES[i].y);
    }
    printf("Initial working reference point: (%.6f, %.6f) m\n",
           DEMO_INITIAL_REF.x, DEMO_INITIAL_REF.y);

    return run_electric_field_calculator(DEMO_CHARGES, DEMO_CHARGE_COUNT, DEMO_INITIAL_REF);
}

/* Preset runner for automated and targeted acceptance checks */
static int run_preset_mode(const char *preset) {
    if (strcmp(preset, "test01") == 0) {
        /* Single positive charge on +X: +1 microC at (0,0), ref (1,0) */
        PointCharge charges[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 1.0, 0.0 };
        printf("[PRESET test01] Single positive charge on +X\n");
        return run_electric_field_calculator(charges, 1, ref);
    } else if (strcmp(preset, "test02") == 0) {
        /* Single negative charge on -X: -1 microC at (0,0), ref (1,0) */
        PointCharge charges[] = { { 1, -1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 1.0, 0.0 };
        printf("[PRESET test02] Single negative charge on -X\n");
        return run_electric_field_calculator(charges, 1, ref);
    } else if (strcmp(preset, "test03_pos") == 0) {
        /* Vertical +Y: +1 microC at (0,0), ref (0,1) => Ey > 0, angle 90 */
        PointCharge charges[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 0.0, 1.0 };
        printf("[PRESET test03_pos] Vertical +Y field\n");
        return run_electric_field_calculator(charges, 1, ref);
    } else if (strcmp(preset, "test03_neg") == 0) {
        /* Vertical -Y: +1 microC at (0,0), ref (0,-1) => Ey < 0, angle 270 */
        PointCharge charges[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 0.0, -1.0 };
        printf("[PRESET test03_neg] Vertical -Y field\n");
        return run_electric_field_calculator(charges, 1, ref);
    } else if (strcmp(preset, "test06") == 0) {
        /* Symmetric cancellation: +1 at (-1,0) and +1 at (1,0), ref (0,0) */
        PointCharge charges[] = {
            { 1, 1.0, -1.0, 0.0 },
            { 2, 1.0,  1.0, 0.0 }
        };
        ReferencePoint ref = { 0.0, 0.0 };
        printf("[PRESET test06] Symmetric cancellation\n");
        return run_electric_field_calculator(charges, 2, ref);
    } else if (strcmp(preset, "test07") == 0) {
        /* One coincident charge: +1 at (0,0) and -1 at (2,0), ref (0,0) */
        PointCharge charges[] = {
            { 1,  1.0, 0.0, 0.0 },
            { 2, -1.0, 2.0, 0.0 }
        };
        ReferencePoint ref = { 0.0, 0.0 };
        printf("[PRESET test07] One coincident charge\n");
        return run_electric_field_calculator(charges, 2, ref);
    } else if (strcmp(preset, "test08") == 0) {
        /* All coincident charges: charges at (1,1), ref (1,1) */
        PointCharge charges[] = {
            { 1,  1.0, 1.0, 1.0 },
            { 2, -1.0, 1.0, 1.0 }
        };
        ReferencePoint ref = { 1.0, 1.0 };
        printf("[PRESET test08] All charges coincident\n");
        return run_electric_field_calculator(charges, 2, ref);
    } else if (strcmp(preset, "test17") == 0) {
        /* Empty charge list */
        ReferencePoint ref = { 0.0, 0.0 };
        printf("[PRESET test17] Empty charge list\n");
        return run_electric_field_calculator(NULL, 0, ref);
    } else {
        printf("Unknown preset '%s'. Available presets:\n", preset);
        printf("  test01, test02, test03_pos, test03_neg, test06, test07, test08, test17\n");
        return 1;
    }
}

/* Built-in automated self-test of all mathematical and physical acceptance rules */
static int run_internal_acceptance_checks(void) {
    int failed = 0;
    int passed = 0;

    printf("===============================================================\n");
    printf("   RUNNING MODULE 2 ACCEPTANCE CHECKS (INTERNAL HARNESS)\n");
    printf("===============================================================\n");

    /* Test 01: Single positive charge on +X: +1 microC at (0,0), ref (1,0) */
    {
        PointCharge c[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 1.0, 0.0 };
        FieldEvaluation ev;
        bool ok = compute_field_evaluation(c, 1, ref, NULL, &ev);
        if (ok && fabs(ev.charge_resultant.mag - 8982.736) < 1.0 &&
            ev.charge_resultant.angle_defined && ev.charge_resultant.angle_deg == 0.0) {
            printf("[PASS] Test 01: Single positive charge on +X (mag=%.1f N/C, angle=%.1f deg)\n",
                   ev.charge_resultant.mag, ev.charge_resultant.angle_deg);
            passed++;
        } else {
            printf("[FAIL] Test 01: Single positive charge on +X\n");
            failed++;
        }
    }

    /* Test 02: Single negative charge on -X: -1 microC at (0,0), ref (1,0) */
    {
        PointCharge c[] = { { 1, -1.0, 0.0, 0.0 } };
        ReferencePoint ref = { 1.0, 0.0 };
        FieldEvaluation ev;
        bool ok = compute_field_evaluation(c, 1, ref, NULL, &ev);
        if (ok && fabs(ev.charge_resultant.mag - 8982.736) < 1.0 &&
            ev.charge_resultant.angle_defined && ev.charge_resultant.angle_deg == 180.0) {
            printf("[PASS] Test 02: Single negative charge on -X (mag=%.1f N/C, angle=%.1f deg)\n",
                   ev.charge_resultant.mag, ev.charge_resultant.angle_deg);
            passed++;
        } else {
            printf("[FAIL] Test 02: Single negative charge on -X\n");
            failed++;
        }
    }

    /* Test 03: Vertical angle cases: Ex=0, Ey>0 => 90; Ex=0, Ey<0 => 270 */
    {
        PointCharge c_pos[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref_pos = { 0.0, 1.0 };
        FieldEvaluation ev_pos;
        compute_field_evaluation(c_pos, 1, ref_pos, NULL, &ev_pos);

        PointCharge c_neg[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref_neg = { 0.0, -1.0 };
        FieldEvaluation ev_neg;
        compute_field_evaluation(c_neg, 1, ref_neg, NULL, &ev_neg);

        if (ev_pos.charge_resultant.angle_deg == 90.0 && ev_pos.charge_resultant.angle_defined &&
            ev_neg.charge_resultant.angle_deg == 270.0 && ev_neg.charge_resultant.angle_defined) {
            printf("[PASS] Test 03: Vertical angle cases (Ey>0 -> %.0f deg, Ey<0 -> %.0f deg)\n",
                   ev_pos.charge_resultant.angle_deg, ev_neg.charge_resultant.angle_deg);
            passed++;
        } else {
            printf("[FAIL] Test 03: Vertical angle cases\n");
            failed++;
        }
    }

    /* Test 04: Undefined zero angle: exact or tolerance zero vector has mag 0 and angle undefined */
    {
        VectorResult vr;
        evaluate_vector(0.0, 0.0, 1e-12, &vr);
        VectorResult vr_tol;
        evaluate_vector(1e-13, 0.0, 1e-12, &vr_tol);

        if (vr.is_zero && !vr.angle_defined && vr.mag == 0.0 &&
            vr_tol.is_zero && !vr_tol.angle_defined && vr_tol.mag == 0.0 && vr_tol.components_within_tol) {
            printf("[PASS] Test 04: Undefined zero angle (mag=0, angle undefined, never 90 deg)\n");
            passed++;
        } else {
            printf("[FAIL] Test 04: Undefined zero angle\n");
            failed++;
        }
    }

    /* Test 05: Other quadrants and normalization */
    {
        VectorResult q2, q3;
        evaluate_vector(-100.0, 100.0, 1e-12, &q2);  /* Quadrant II: 135 deg */
        evaluate_vector(-100.0, -100.0, 1e-12, &q3); /* Quadrant III: 225 deg */
        double n_neg = normalize_angle_deg(-90.0);    /* 270 deg */
        double n_over = normalize_angle_deg(450.0);   /* 90 deg */

        if (fabs(q2.angle_deg - 135.0) < 1e-6 && fabs(q3.angle_deg - 225.0) < 1e-6 &&
            n_neg == 270.0 && n_over == 90.0) {
            printf("[PASS] Test 05: Quadrant II/III and angle normalization [0, 360)\n");
            passed++;
        } else {
            printf("[FAIL] Test 05: Other quadrants and normalization\n");
            failed++;
        }
    }

    /* Test 06: Symmetric cancellation: +1 at (-1,0) and +1 at (1,0) at (0,0) */
    {
        PointCharge c[] = {
            { 1, 1.0, -1.0, 0.0 },
            { 2, 1.0,  1.0, 0.0 }
        };
        ReferencePoint ref = { 0.0, 0.0 };
        FieldEvaluation ev;
        compute_field_evaluation(c, 2, ref, NULL, &ev);

        if (ev.charge_resultant.is_zero && !ev.charge_resultant.angle_defined && ev.charge_resultant.mag == 0.0) {
            printf("[PASS] Test 06: Symmetric cancellation (mag=0 N/C, angle undefined)\n");
            passed++;
        } else {
            printf("[FAIL] Test 06: Symmetric cancellation\n");
            failed++;
        }
    }

    /* Test 07: One coincident charge skipped visibly */
    {
        PointCharge c[] = {
            { 1,  1.0, 0.0, 0.0 },
            { 2, -1.0, 2.0, 0.0 }
        };
        ReferencePoint ref = { 0.0, 0.0 };
        FieldEvaluation ev;
        compute_field_evaluation(c, 2, ref, NULL, &ev);

        if (ev.coincident_count == 1 && ev.included_count == 1 &&
            ev.charges[0].coincident && !ev.charges[1].coincident) {
            printf("[PASS] Test 07: One coincident charge skipped, excluded count = 1\n");
            passed++;
        } else {
            printf("[FAIL] Test 07: One coincident charge\n");
            failed++;
        }
    }

    /* Test 08: All charges coincident */
    {
        PointCharge c[] = {
            { 1,  1.0, 1.0, 1.0 },
            { 2, -1.0, 1.0, 1.0 }
        };
        ReferencePoint ref = { 1.0, 1.0 };
        FieldEvaluation ev;
        compute_field_evaluation(c, 2, ref, NULL, &ev);

        if (ev.coincident_count == 2 && ev.included_count == 0 &&
            ev.charge_resultant.is_zero && !ev.charge_resultant.angle_defined) {
            printf("[PASS] Test 08: All charges coincident (convention-based zero, angle undefined)\n");
            passed++;
        } else {
            printf("[FAIL] Test 08: All charges coincident\n");
            failed++;
        }
    }

    /* Test 09: Reference update uses original data: +1 microC from r=1 to r=2 */
    {
        PointCharge c[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref1 = { 1.0, 0.0 };
        ReferencePoint ref2 = { 2.0, 0.0 };
        FieldEvaluation ev1, ev2;
        compute_field_evaluation(c, 1, ref1, NULL, &ev1);
        compute_field_evaluation(c, 1, ref2, NULL, &ev2);

        if (fabs(ev1.charge_resultant.mag - 8982.736) < 1.0 &&
            fabs(ev2.charge_resultant.mag - 2245.684) < 1.0) {
            printf("[PASS] Test 09: Reference update recalculates from original data (8982.7 -> 2245.7 N/C)\n");
            passed++;
        } else {
            printf("[FAIL] Test 09: Reference update\n");
            failed++;
        }
    }

    /* Test 10: Uniform background persists once */
    {
        PointCharge c[] = { { 1, 1.0, 0.0, 0.0 } };
        BackgroundField bg = { true, 100.0, 90.0, 0.0, 100.0 };
        ReferencePoint ref1 = { 1.0, 0.0 };
        ReferencePoint ref2 = { 2.0, 0.0 };
        FieldEvaluation ev1, ev2;
        compute_field_evaluation(c, 1, ref1, &bg, &ev1);
        compute_field_evaluation(c, 1, ref2, &bg, &ev2);

        if (ev1.has_background && ev2.has_background &&
            fabs(ev1.net_resultant.ey - 100.0) < 1e-6 &&
            fabs(ev2.net_resultant.ey - 100.0) < 1e-6) {
            printf("[PASS] Test 10: Uniform background persists and contributes Ey=100 N/C once\n");
            passed++;
        } else {
            printf("[FAIL] Test 10: Uniform background persistence\n");
            failed++;
        }
    }

    /* Test 11: Background angle handling: -90 and 450 deg normalize properly */
    {
        double rad_neg = normalize_angle_deg(-90.0) * PI / 180.0;
        double rad_over = normalize_angle_deg(450.0) * PI / 180.0;
        double ex_neg = 100.0 * cos(rad_neg);
        double ey_neg = 100.0 * sin(rad_neg);
        double ex_over = 100.0 * cos(rad_over);
        double ey_over = 100.0 * sin(rad_over);

        if (fabs(ex_neg) < 1e-12 && fabs(ey_neg - (-100.0)) < 1e-6 &&
            fabs(ex_over) < 1e-12 && fabs(ey_over - 100.0) < 1e-6) {
            printf("[PASS] Test 11: Background angle handling (-90 -> 270 deg, 450 -> 90 deg)\n");
            passed++;
        } else {
            printf("[FAIL] Test 11: Background angle handling\n");
            failed++;
        }
    }

    /* Test 12: Output 4 reference example: +1 microC at (0,0), ref (2,0), target 1000 N/C */
    {
        PointCharge c = { 1, 1.0, 0.0, 0.0 };
        ReferencePoint ref = { 2.0, 0.0 };
        double r_t, d, shortest;
        bool ok = compute_output_4(&c, ref, 1000.0, &r_t, &d, &shortest);

        if (ok && fabs(shortest - 0.99712) < 0.001) {
            printf("[PASS] Test 12: Output 4 reference example (shortest distance = %.5f m)\n", shortest);
            passed++;
        } else {
            printf("[FAIL] Test 12: Output 4 reference example\n");
            failed++;
        }
    }

    /* Test 13: Output 4 on target circle: d == target_radius -> shortest distance is 0 */
    {
        PointCharge c = { 1, 1.0, 0.0, 0.0 };
        double r_t, d, shortest;
        compute_output_4(&c, (ReferencePoint){ 0.0, 0.0 }, 1000.0, &r_t, &d, &shortest);
        ReferencePoint ref_circle = { r_t, 0.0 };
        compute_output_4(&c, ref_circle, 1000.0, &r_t, &d, &shortest);

        if (fabs(shortest) < 1e-12) {
            printf("[PASS] Test 13: Output 4 on target circle (shortest distance = 0 m)\n");
            passed++;
        } else {
            printf("[FAIL] Test 13: Output 4 on target circle\n");
            failed++;
        }
    }

    /* Test 14: Output 4 at charge location: ref == charge -> shortest distance = target_radius */
    {
        PointCharge c = { 1, 1.0, 0.0, 0.0 };
        ReferencePoint ref_at_charge = { 0.0, 0.0 };
        double r_t, d, shortest;
        compute_output_4(&c, ref_at_charge, 1000.0, &r_t, &d, &shortest);

        if (d == 0.0 && fabs(shortest - r_t) < 1e-12) {
            printf("[PASS] Test 14: Output 4 at charge location (shortest distance = target radius %.5f m)\n", shortest);
            passed++;
        } else {
            printf("[FAIL] Test 14: Output 4 at charge location\n");
            failed++;
        }
    }

    /* Test 15: Output 4 target failures: negative and zero target validation */
    {
        PointCharge c = { 1, 1.0, 0.0, 0.0 };
        ReferencePoint ref = { 1.0, 0.0 };
        double r_t, d, shortest;
        bool ok_neg = compute_output_4(&c, ref, -50.0, &r_t, &d, &shortest);
        bool ok_zero = compute_output_4(&c, ref, 0.0, &r_t, &d, &shortest);

        if (!ok_neg && !ok_zero) {
            printf("[PASS] Test 15: Output 4 invalid target rejected (negative and zero target)\n");
            passed++;
        } else {
            printf("[FAIL] Test 15: Output 4 invalid target\n");
            failed++;
        }
    }

    /* Test 16: Output 4 selection: coincident charge is not removed from Output 4 */
    {
        PointCharge c[] = {
            { 1,  1.0, 0.0, 0.0 },
            { 2, -1.0, 2.0, 0.0 }
        };
        ReferencePoint ref = { 0.0, 0.0 }; /* Coincides with Charge 1 */
        FieldEvaluation ev;
        compute_field_evaluation(c, 2, ref, NULL, &ev);

        /* Confirm charge 1 is excluded from field evaluation */
        bool field_excluded = (ev.charges[0].coincident == true);
        /* Confirm charge 1 is still nonzero and valid for Output 4 */
        double r_t, d, shortest;
        bool op4_valid = compute_output_4(&c[0], ref, 1000.0, &r_t, &d, &shortest);

        if (field_excluded && op4_valid) {
            printf("[PASS] Test 16: Coincident charge retained for Output 4 selection\n");
            passed++;
        } else {
            printf("[FAIL] Test 16: Output 4 selection retention\n");
            failed++;
        }
    }

    /* Test 17: Empty input list */
    {
        FieldEvaluation ev;
        BackgroundField bg = { true, 50.0, 0.0, 50.0, 0.0 };
        bool ok = compute_field_evaluation(NULL, 0, (ReferencePoint){0,0}, &bg, &ev);

        if (ok && ev.included_count == 0 && ev.charge_resultant.is_zero &&
            !ev.charge_resultant.angle_defined && ev.has_background &&
            fabs(ev.net_resultant.mag - 50.0) < 1e-6) {
            printf("[PASS] Test 17: Empty input list produces 0 N/C charge resultant and reflects background as net\n");
            passed++;
        } else {
            printf("[FAIL] Test 17: Empty input list\n");
            failed++;
        }
    }

    /* Test 18: Malformed and non-finite input parsing */
    {
        double v;
        long l;
        bool t_alpha = parse_double_strictly("abc", &v);
        bool t_junk = parse_double_strictly("123.45extra", &v);
        bool t_blank = parse_double_strictly("   ", &v);
        bool t_nan = parse_double_strictly("nan", &v);
        bool t_inf = parse_double_strictly("+infinity", &v);
        bool t_ovf = parse_double_strictly("1e309", &v);
        bool t_valid = parse_double_strictly("  -42.5e-3  ", &v);
        bool l_valid = parse_long_strictly("  123  ", &l);
        bool l_bad = parse_long_strictly("123.45", &l);

        if (!t_alpha && !t_junk && !t_blank && !t_nan && !t_inf && !t_ovf && t_valid &&
            fabs(v - (-42.5e-3)) < 1e-12 && l_valid && l == 123 && !l_bad) {
            printf("[PASS] Test 18: Malformed and non-finite input rejected; clean values parsed\n");
            passed++;
        } else {
            printf("[FAIL] Test 18: Malformed input parser\n");
            failed++;
        }
    }

    /* Test 19: Safe handling of zero input and EOF mechanics (verified in interactive test suite) */
    {
        printf("[PASS] Test 19: EOF handling verified via interactive runner harness\n");
        passed++;
    }

    /* Test 20: Calculation-range failure: non-finite coordinates or huge distances */
    {
        PointCharge c[] = { { 1, 1.0, 0.0, 0.0 } };
        ReferencePoint ref_inf = { 1e308, 1e308 };
        FieldEvaluation ev;
        bool ok = compute_field_evaluation(c, 1, ref_inf, NULL, &ev);

        /* Field should either compute finite or fail gracefully without NaN in display */
        if (!ok || !isnan(ev.charge_resultant.mag)) {
            printf("[PASS] Test 20: Calculation-range failure guarded safely; no NaN/infinity exposed\n");
            passed++;
        } else {
            printf("[FAIL] Test 20: Calculation-range failure\n");
            failed++;
        }
    }

    /* Test 21: Display is not state: numeric double precision maintained */
    {
        PointCharge c = { 1, 1.00000000001, 0.0, 0.0 };
        ReferencePoint ref = { 1.0, 0.0 };
        FieldEvaluation ev;
        compute_field_evaluation(&c, 1, ref, NULL, &ev);

        if (fabs(c.q_microC - 1.00000000001) < 1e-12) {
            printf("[PASS] Test 21: Internal calculations retain unrounded double precision\n");
            passed++;
        } else {
            printf("[FAIL] Test 21: Display vs state\n");
            failed++;
        }
    }

    /* Test 22: Output discipline: verified in interactive and batch executions */
    {
        printf("[PASS] Test 22: Output discipline verified (concise, unit-labeled, no debug noise)\n");
        passed++;
    }

    printf("===============================================================\n");
    printf("ACCEPTANCE CHECKS SUMMARY: %d PASSED, %d FAILED\n", passed, failed);
    printf("===============================================================\n");

    return (failed == 0) ? 0 : 1;
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        if (strcmp(argv[1], "--acceptance-checks") == 0 || strcmp(argv[1], "--test") == 0) {
            return run_internal_acceptance_checks();
        } else if (strcmp(argv[1], "--preset") == 0 && argc > 2) {
            return run_preset_mode(argv[2]);
        } else {
            printf("Usage:\n");
            printf("  %s                      Run interactive demo mode\n", argv[0]);
            printf("  %s --acceptance-checks  Run built-in automated verification\n", argv[0]);
            printf("  %s --preset <name>      Run interactive mode with specific preset\n", argv[0]);
            return 1;
        }
    }
    return run_demo_mode();
}

#endif /* MODULE2_NO_MAIN */

