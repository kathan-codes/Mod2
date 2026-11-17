/*
 * Module 2: Electric Field Calculator (ISO C99, dry air at 20 °C, 101.325 kPa)
 * Zero tolerances: vector <= 1e-12 N/C; charge sum <= max(1e-12, 1e-12*sum_mags);
 * net sum <= max(1e-12, 1e-12*(sum_mags + bg_mag)).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>

#define EPS0  8.8541878188e-12
#define EPSR  1.0005364
#define PI    3.14159265358979323846
#define K_AIR (1.0 / (4.0 * PI * EPS0 * EPSR))
#define MAX_C 256
#define BUF_SZ 256

typedef struct { int id; double q_microC, x, y; } PointCharge;
typedef struct { double x, y; } ReferencePoint;
typedef struct { double ex, ey, mag, angle_deg; bool is_zero, angle_defined, within_tol; } VectorResult;
typedef struct { int id; double q_microC, x, y; bool coincident; VectorResult field; } EvaluatedCharge;
typedef struct { bool active; double magnitude, angle_deg, ex, ey; } BackgroundField;
typedef struct {
    EvaluatedCharge charges[MAX_C];
    size_t total_count, coincident_count, included_count;
    VectorResult charge_resultant, background, net_resultant;
    bool has_background;
} FieldEvaluation;

static const PointCharge DEMO_C[] = { { 1, 1.0, 0, 0 }, { 2, -1.0, 2, 0 } };
static const ReferencePoint DEMO_REF = { 1.0, 1.0 };

static const char *skip_sp(const char *s) { while (s && isspace((unsigned char)*s)) s++; return s; }
static bool ci_match(const char *s, const char *p) {
    while (*p) if (tolower((unsigned char)*s++) != tolower((unsigned char)*p++)) return false;
    return true;
}
static int read_line(char *b, size_t sz, FILE *fp) {
    if (!sz || !fgets(b, (int)sz, fp)) return -1;
    size_t n = strlen(b);
    if (n && b[n - 1] == '\n') { b[--n] = '\0'; if (n && b[n - 1] == '\r') b[--n] = '\0'; return 0; }
    if (!feof(fp)) { int ch; while ((ch = fgetc(fp)) != '\n' && ch != EOF); return 1; }
    return 0;
}
static bool parse_double_strictly(const char *b, double *out) {
    b = skip_sp(b); if (!b || !*b) return false;
    const char *chk = (*b == '+' || *b == '-') ? b + 1 : b;
    if (ci_match(chk, "nan") || ci_match(chk, "inf")) return false;
    errno = 0; char *e = NULL; double v = strtod(b, &e);
    if (e == b || errno == ERANGE || !isfinite(v) || *skip_sp(e)) return false;
    *out = v; return true;
}
static bool parse_long_strictly(const char *b, long *out) {
    b = skip_sp(b); if (!b || !*b) return false;
    errno = 0; char *e = NULL; long v = strtol(b, &e, 10);
    if (e == b || errno == ERANGE || *skip_sp(e)) return false;
    *out = v; return true;
}
static int prompt_double(const char *msg, double *out) {
    char b[BUF_SZ];
    for (;;) {
        printf("%s", msg); fflush(stdout);
        int st = read_line(b, sizeof(b), stdin);
        if (st < 0) return -1;
        if (st > 0 || !parse_double_strictly(b, out)) {
            printf(st > 0 ? "Error: Input line exceeds buffer limit. Please enter a valid number.\n" : "Error: Invalid numeric input. Please enter a finite real number.\n");
            continue;
        }
        return 0;
    }
}

static double norm_angle(double d) {
    d = fmod(d, 360.0); if (d < 0.0) d += 360.0;
    return (d >= 360.0 || fabs(d - 360.0) < 1e-12 || fabs(d) < 1e-12) ? 0.0 : d;
}
static bool eval_vec(double ex, double ey, double tol, VectorResult *out) {
    if (!isfinite(ex) || !isfinite(ey)) return false;
    double m = hypot(ex, ey);
    if (!isfinite(m)) return false;
    if (m <= tol) *out = (VectorResult){ ex, ey, 0, 0, true, false, (m > 0) };
    else *out = (VectorResult){ ex, ey, m, (ex == 0.0) ? ((ey > 0.0) ? 90.0 : 270.0) : norm_angle(atan2(ey, ex) * 180.0 / PI), false, true, false };
    return true;
}

static bool compute_field_evaluation(const PointCharge *c, size_t count, ReferencePoint ref,
                                     const BackgroundField *bg, FieldEvaluation *out) {
    if (count > MAX_C) return false;
    memset(out, 0, sizeof(*out));
    size_t n = 0;
    while (c && n < count && c[n].q_microC != 0.0) n++;
    out->total_count = n;
    double sx = 0, sy = 0, sm = 0;
    for (size_t i = 0; i < n; i++) {
        out->charges[i].id = (c[i].id > 0) ? c[i].id : (int)(i + 1);
        out->charges[i].q_microC = c[i].q_microC; out->charges[i].x = c[i].x; out->charges[i].y = c[i].y;
        if (c[i].x == ref.x && c[i].y == ref.y) {
            out->charges[i].coincident = true; out->coincident_count++;
            out->charges[i].field.is_zero = true; continue;
        }
        out->included_count++;
        double dx = ref.x - c[i].x, dy = ref.y - c[i].y, r = hypot(dx, dy);
        if (!isfinite(r) || r == 0.0) return false;
        double f = (K_AIR * (c[i].q_microC * 1e-6) / r) / (r * r);
        if (!isfinite(f) || !eval_vec(f * dx, f * dy, 1e-12, &out->charges[i].field)) return false;
        sx += f * dx; sy += f * dy; sm += hypot(f * dx, f * dy);
    }
    if (!out->included_count) out->charge_resultant = (VectorResult){ 0, 0, 0, 0, true, false, false };
    else if (!eval_vec(sx, sy, fmax(1e-12, 1e-12 * sm), &out->charge_resultant)) return false;
    if (bg && bg->active) {
        out->has_background = true;
        out->background = (VectorResult){ bg->ex, bg->ey, bg->magnitude, bg->angle_deg, bg->magnitude == 0, bg->magnitude > 0, false };
        if (!eval_vec(sx + bg->ex, sy + bg->ey, fmax(1e-12, 1e-12 * (sm + bg->magnitude)), &out->net_resultant)) return false;
    }
    return true;
}

static void print_vec(const char *ind, const VectorResult *v) {
    printf("%sEx: %+.6e N/C, Ey: %+.6e N/C%s\n", ind, v->ex, v->ey, v->within_tol ? " (components within zero tolerance)" : "");
    if (v->is_zero) printf("%sMagnitude: 0 N/C, Angle: undefined\n", ind);
    else printf("%sMagnitude: %.6e N/C, Angle: %.2f degrees\n", ind, v->mag, v->angle_deg);
}

static void print_all(const FieldEvaluation *ev, ReferencePoint ref, bool show_breakdown) {
    if (show_breakdown) {
        printf("\n--- Point Charge Field Contributions (Reference Point: (%.6f, %.6f) m) ---\n", ref.x, ref.y);
        if (!ev->total_count) printf("No point charges present in dataset.\n");
        for (size_t i = 0; i < ev->total_count; i++) {
            const EvaluatedCharge *ec = &ev->charges[i];
            if (ec->coincident) printf("Charge %d (%+.6f microC at (%.6f, %.6f) m): skipped because it coincides with the reference point; its field is undefined here.\n", ec->id, ec->q_microC, ec->x, ec->y);
            else { printf("Charge %d (%+.6f microC at (%.6f, %.6f) m):\n", ec->id, ec->q_microC, ec->x, ec->y); print_vec("  ", &ec->field); }
        }
        printf("Excluded charges: %zu\n", ev->coincident_count);
    }
    if (!ev->included_count) {
        printf("\n--- Charge-Only Resultant Field%s ---\nNote: %s\n",
               ev->coincident_count ? " (Convention-based sum from 0 included charges)" : "",
               ev->coincident_count ? "All point charges coincide with the reference point; the physical point-charge field at this reference point is undefined." : "No charges in dataset.");
    } else if (ev->coincident_count) {
        printf("\n--- Charge-Only Resultant Field (Convention-based sum over %zu included charges) ---\nNote: %zu charge(s) coincided with reference point; physical total point-charge field is undefined at this point.\n", ev->included_count, ev->coincident_count);
    } else printf("\n--- Charge-Only Resultant Field ---\n");
    print_vec("", &ev->charge_resultant);

    if (ev->has_background) {
        printf("\n--- Uniform Background Field ---\n"); print_vec("", &ev->background);
        printf("\n--- Net Resultant Field%s ---\n%s", ev->coincident_count ? " (Convention-based sum over included charges + background)" : "",
               ev->coincident_count ? "Note: Physical total field is undefined at this point due to coincident charge(s).\n" : "");
        print_vec("", &ev->net_resultant);
    }
}

static bool compute_output_4(const PointCharge *c, ReferencePoint ref, double e_tgt, double *r, double *d, double *s) {
    if (!c || e_tgt <= 0.0 || c->q_microC == 0.0) return false;
    double r_sq = (K_AIR * fabs(c->q_microC * 1e-6)) / e_tgt;
    if (!isfinite(r_sq) || r_sq < 0.0) return false;
    *r = sqrt(r_sq); *d = hypot(ref.x - c->x, ref.y - c->y); *s = fabs(*d - *r);
    return isfinite(*r) && isfinite(*d) && isfinite(*s);
}

static int handle_output_4(const PointCharge *c, size_t count, ReferencePoint ref) {
    size_t elig[MAX_C], n = 0;
    for (size_t i = 0; i < count; i++) if (c[i].q_microC != 0.0) elig[n++] = i;
    if (!n) { printf("\nOutput 4 cannot run: No nonzero charge is available.\n"); return 0; }
    const PointCharge *sel = &c[elig[0]];
    if (n == 1) printf("\nUsing single eligible charge: ID %d (%+.6f microC at (%.6f, %.6f) m)\n", sel->id, sel->q_microC, sel->x, sel->y);
    else {
        printf("\nEligible nonzero charges for Output 4:\n");
        for (size_t k = 0; k < n; k++) printf("  ID %d: %+.6f microC at (%.6f, %.6f) m\n", c[elig[k]].id, c[elig[k]].q_microC, c[elig[k]].x, c[elig[k]].y);
        char b[BUF_SZ];
        for (sel = NULL; !sel;) {
            printf("Select charge by ID: "); fflush(stdout);
            int st = read_line(b, sizeof(b), stdin);
            if (st < 0) return -1;
            if (st > 0) { printf("Error: Input line exceeds buffer limit. Please enter a valid integer.\n"); continue; }
            long id;
            if (parse_long_strictly(b, &id)) {
                for (size_t k = 0; k < n; k++) if (c[elig[k]].id == (int)id) { sel = &c[elig[k]]; break; }
            }
            if (!sel) printf("Error: Invalid selection. Please choose an ID from the listed eligible charges.\n");
        }
    }
    double e_tgt;
    for (;;) {
        if (prompt_double("Enter target electric field magnitude (N/C): ", &e_tgt) < 0) return -1;
        if (e_tgt < 0.0) { printf("Error: Target field magnitude cannot be negative. Please enter a nonnegative value.\n"); continue; }
        if (e_tgt == 0.0) { printf("For a nonzero point charge, an electric field magnitude of zero has no finite matching location (field approaches zero only as distance approaches infinity).\n"); return 0; }
        break;
    }
    double r, d, s;
    if (!compute_output_4(sel, ref, e_tgt, &r, &d, &s)) { printf("Error: Calculation-range error during shortest distance evaluation.\n"); return 0; }
    printf("\n--- Output 4 Result ---\nSelected Charge: ID %d (%+.6f microC at (%.6f, %.6f) m)\nTarget Field Magnitude: %.6e N/C\nTarget Field Circle Radius: %.6e m\nDistance from Reference Point to Charge: %.6e m\n", sel->id, sel->q_microC, sel->x, sel->y, e_tgt, r, d);
    if (fabs(d - r) < 1e-9) printf("Shortest Distance to Target Field: 0 m (reference point lies on target circle)\n");
    else if (d == 0.0) printf("Shortest Distance to Target Field: %.5f m (reference point equals charge location; distance is defined, no unique direction)\n", s);
    else printf("Shortest Distance to Target Field: %.5f m\n", s);
    return 0;
}

int run_electric_field_calculator(const PointCharge *charges, size_t count, ReferencePoint initial_ref) {
    if (!charges && count > 0) { printf("Error: Input data is unavailable from Module 1.\n"); return -1; }
    if (count > MAX_C) { printf("Error: Charge list exceeds maximum supported capacity (%d).\n", MAX_C); return -1; }
    size_t n = 0;
    while (charges && n < count && charges[n].q_microC != 0.0) n++;
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(charges[i].q_microC) || !isfinite(charges[i].x) || !isfinite(charges[i].y)) {
            printf("Error: Input data contains non-finite values.\n"); return -1;
        }
    }
    if (!isfinite(initial_ref.x) || !isfinite(initial_ref.y)) { printf("Error: Initial reference point contains non-finite values.\n"); return -1; }

    ReferencePoint working_ref = initial_ref;
    BackgroundField bg = { false, 0, 0, 0, 0 };
    FieldEvaluation eval;
    if (!compute_field_evaluation(charges, n, working_ref, &bg, &eval)) {
        printf("Error: Calculation-range error during initial field evaluation.\n"); return -1;
    }
    print_all(&eval, working_ref, true);
    char buf[BUF_SZ];
    for (;;) {
        printf("\n--- Menu ---\n1  Change the reference point\n2  Set, replace, or clear the uniform background field\n3  Find the shortest distance for a target field (Output 4)\n0  Exit\nEnter choice: ");
        fflush(stdout);
        int rls = read_line(buf, sizeof(buf), stdin);
        if (rls < 0) { printf("\nInput ended (EOF). Exiting cleanly.\n"); break; }
        if (rls > 0) { printf("Error: Input line exceeds buffer limit. Valid choices are 0, 1, 2, or 3.\n"); continue; }
        long choice;
        if (!parse_long_strictly(buf, &choice) || choice < 0 || choice > 3) {
            printf("Error: Invalid choice. Valid choices are 0, 1, 2, or 3.\n"); continue;
        }
        if (choice == 0) { printf("Exiting Electric Field Calculator.\n"); break; }
        if (choice == 1) {
            double nx, ny;
            if (prompt_double("Enter new X reference coordinate (m): ", &nx) < 0 || prompt_double("Enter new Y reference coordinate (m): ", &ny) < 0) {
                printf("\nInput ended (EOF). Reference point change canceled.\n"); break;
            }
            FieldEvaluation t_eval;
            if (!compute_field_evaluation(charges, n, (ReferencePoint){nx, ny}, &bg, &t_eval)) {
                printf("Error: Calculation-range error with proposed reference point. Changes not committed.\n"); continue;
            }
            working_ref = (ReferencePoint){nx, ny}; eval = t_eval;
            print_all(&eval, working_ref, true);
        } else if (choice == 2) {
            double mag, ang = 0;
            for (;;) {
                if (prompt_double("Enter uniform background field magnitude (N/C, 0 to clear): ", &mag) < 0) {
                    printf("\nInput ended (EOF). Background field update canceled.\n"); return 0;
                }
                if (mag >= 0.0) break;
                printf("Error: Background field magnitude must be nonnegative.\n");
            }
            if (mag > 0.0 && prompt_double("Enter background field angle (degrees counterclockwise from +X): ", &ang) < 0) {
                printf("\nInput ended (EOF). Background field update canceled.\n"); return 0;
            }
            double na = norm_angle(ang), rad = na * PI / 180.0;
            BackgroundField t_bg = { mag > 0.0, mag, na, mag * cos(rad), mag * sin(rad) };
            FieldEvaluation t_eval;
            if (!compute_field_evaluation(charges, n, working_ref, &t_bg, &t_eval)) {
                printf("Error: Calculation-range error applying background field. Changes not committed.\n"); continue;
            }
            bg = t_bg; eval = t_eval;
            if (!bg.active) printf("\nUniform background field cleared.\n");
            print_all(&eval, working_ref, false);
        } else if (choice == 3) {
            if (handle_output_4(charges, n, working_ref) < 0) {
                printf("\nInput ended (EOF). Output 4 canceled.\n"); break;
            }
        }
    }
    return 0;
}

#ifndef MODULE2_NO_MAIN

static int run_demo_mode(void) {
    printf("===============================================================\n   ELECTRIC FIELD CALCULATOR (MODULE 2) [DEMO MODE]\n   Medium: Dry air (20 deg C, 101.325 kPa)\n   k_air: %.4e N*m^2/C^2\n===============================================================\nUsing clearly labeled reproducible demo dataset:\n", K_AIR);
    for (size_t i = 0; i < 2; i++) printf("  Charge %d: %+.6f microC at (%.6f, %.6f) m\n", DEMO_C[i].id, DEMO_C[i].q_microC, DEMO_C[i].x, DEMO_C[i].y);
    printf("Initial working reference point: (%.6f, %.6f) m\n", DEMO_REF.x, DEMO_REF.y);
    return run_electric_field_calculator(DEMO_C, 2, DEMO_REF);
}

static int run_preset_mode(const char *p) {
    if (!strcmp(p, "test01")) return run_electric_field_calculator(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){1, 0});
    if (!strcmp(p, "test12")) return run_electric_field_calculator(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){2, 0});
    return 1;
}

static bool run_proc(const char *a, const char *in, char *out, size_t sz, int *c) {
    FILE *f = fopen("_i.tmp", "wb"); if (!f) return false;
    if (in) { fputs(in, f); }
    fclose(f);
    char cmd[256]; snprintf(cmd, sizeof(cmd), "electric_field_calculator.exe %s < _i.tmp > _o.tmp 2>&1", a ? a : "");
    if (c) *c = system(cmd);
    f = fopen("_o.tmp", "rb"); if (!f) { remove("_i.tmp"); return false; }
    size_t n = fread(out, 1, sz - 1, f); out[n] = '\0'; fclose(f);
    remove("_i.tmp"); remove("_o.tmp"); return true;
}

#define CHK(cond, n, desc) do { \
    if (cond) { printf("[PASS] Test %02d: %s\n", n, desc); passed++; } \
    else      { printf("[FAIL] Test %02d: %s\n", n, desc); failed++; } \
} while (0)

static int run_internal_acceptance_checks(void) {
    int failed = 0, passed = 0, code;
    FieldEvaluation ev, ev2, ev_cp, ev_cn; VectorResult vr, vr2;
    char out[4096], lin[512]; double rt, dt, st, v18; long l18;
    PointCharge c_canc[] = { {1, 1, -1, 0}, {2, 1, 1, 0} };

    printf("===============================================================\n   RUNNING MODULE 2 ACCEPTANCE CHECKS (23/23 C-ONLY HARNESS)\n===============================================================\n");
    compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){1, 0}, NULL, &ev); CHK(fabs(ev.charge_resultant.mag - 8982.736) < 1.0 && ev.charge_resultant.angle_deg == 0.0, 1, "Single positive charge on +X (mag=8982.7 N/C, angle=0.0 deg)");
    compute_field_evaluation(&(PointCharge){1, -1, 0, 0}, 1, (ReferencePoint){1, 0}, NULL, &ev); CHK(fabs(ev.charge_resultant.mag - 8982.736) < 1.0 && ev.charge_resultant.angle_deg == 180.0, 2, "Single negative charge on -X (mag=8982.7 N/C, angle=180.0 deg)");
    compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){0, 1}, NULL, &ev); compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){0, -1}, NULL, &ev2);
    compute_field_evaluation(c_canc, 2, (ReferencePoint){0, 1}, NULL, &ev_cp); compute_field_evaluation(c_canc, 2, (ReferencePoint){0, -1}, NULL, &ev_cn);
    CHK(ev.charges[0].field.angle_deg == 90.0 && ev2.charges[0].field.angle_deg == 270.0 && ev_cp.charge_resultant.angle_deg == 90.0 && ev_cn.charge_resultant.angle_deg == 270.0, 3, "Vertical angle cases (Ey>0 -> 90 deg, Ey<0 -> 270 deg)");
    eval_vec(0, 0, 1e-12, &vr); eval_vec(1e-13, 0, 1e-12, &vr2); CHK(vr.is_zero && !vr.angle_defined && vr2.is_zero && !vr2.angle_defined && vr2.within_tol, 4, "Undefined zero angle (exact & tolerance zero: mag 0, angle undefined)");
    eval_vec(-100, 100, 1e-12, &vr); eval_vec(-100, -100, 1e-12, &vr2); CHK(fabs(vr.angle_deg - 135.0) < 1e-6 && fabs(vr2.angle_deg - 225.0) < 1e-6 && norm_angle(-90) == 270.0 && norm_angle(450) == 90.0, 5, "Other quadrants and angle normalization [0, 360)");
    compute_field_evaluation(c_canc, 2, (ReferencePoint){0, 0}, NULL, &ev); CHK(ev.charge_resultant.is_zero && !ev.charge_resultant.angle_defined, 6, "Symmetric cancellation (mag=0 N/C, angle undefined)");
    compute_field_evaluation((PointCharge[]){{1, 1, 0, 0}, {2, -1, 2, 0}}, 2, (ReferencePoint){0, 0}, NULL, &ev); CHK(ev.coincident_count == 1 && ev.included_count == 1 && ev.charges[0].coincident && !ev.charges[1].coincident, 7, "One coincident charge skipped, excluded count = 1");
    compute_field_evaluation((PointCharge[]){{1, 1, 1, 1}, {2, -1, 1, 1}}, 2, (ReferencePoint){1, 1}, NULL, &ev); CHK(ev.coincident_count == 2 && ev.included_count == 0 && ev.charge_resultant.is_zero, 8, "All charges coincident (convention-based zero, angle undefined)");
    compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){1, 0}, NULL, &ev); compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){2, 0}, NULL, &ev2); CHK(fabs(ev.charge_resultant.mag - 8982.736) < 1.0 && fabs(ev2.charge_resultant.mag - 2245.684) < 1.0, 9, "Reference update recalculates from original data");
    compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){1, 0}, &(BackgroundField){true, 100, 90, 0, 100}, &ev); compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){2, 0}, &(BackgroundField){true, 100, 90, 0, 100}, &ev2); CHK(fabs(ev.net_resultant.ey - 100.0) < 1e-6 && fabs(ev2.net_resultant.ey - 100.0) < 1e-6, 10, "Uniform background field superposition");
    double ex11 = 100.0 * cos(norm_angle(-45.0) * PI / 180.0), ey11 = 100.0 * sin(norm_angle(-45.0) * PI / 180.0); compute_field_evaluation(&(PointCharge){1, 1, 0, 0}, 1, (ReferencePoint){1, 0}, &(BackgroundField){true, 100.0, 315.0, ex11, ey11}, &ev); CHK(fabs(ex11 - 70.710678) < 1e-4 && fabs(ey11 + 70.710678) < 1e-4 && fabs(ev.charges[0].field.mag - 8982.736) < 1.0, 11, "Background angle normalization (-45 -> 315 deg)");
    compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){2, 0}, 1000.0, &rt, &dt, &st); CHK(fabs(st - 0.99712) < 0.001, 12, "Output 4 reference example (shortest distance = 0.99712 m)");
    compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){rt, 0}, 1000.0, &rt, &dt, &st); CHK(fabs(st) < 1e-12, 13, "Output 4 on target circle (shortest distance = 0 m)");
    compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){0, 0}, 1000.0, &rt, &dt, &st); CHK(dt == 0.0 && fabs(st - rt) < 1e-12, 14, "Output 4 at charge location (distance = radius, no unique dir)");
    CHK(!compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){1, 0}, -50.0, &rt, &dt, &st) && !compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){1, 0}, 0.0, &rt, &dt, &st) && !compute_output_4(&(PointCharge){2, 0, 0, 0}, (ReferencePoint){1, 0}, 1000.0, &rt, &dt, &st), 15, "Output 4 invalid targets (negative, zero, zero charge)");
    compute_field_evaluation((PointCharge[]){{1, 1, 0, 0}, {2, -1, 2, 0}}, 2, (ReferencePoint){0, 0}, NULL, &ev); CHK(ev.charges[0].coincident && compute_output_4(&(PointCharge){1, 1, 0, 0}, (ReferencePoint){0, 0}, 1000.0, &rt, &dt, &st), 16, "Output 4 allows coincident charge");
    compute_field_evaluation(NULL, 0, (ReferencePoint){0, 0}, &(BackgroundField){true, 50, 0, 50, 0}, &ev); CHK(ev.included_count == 0 && ev.charge_resultant.is_zero && fabs(ev.net_resultant.mag - 50.0) < 1e-6, 17, "Empty charge list produces 0 N/C resultant");
    memset(lin, '9', 300); lin[300] = '\n'; lin[301] = '0'; lin[302] = '\n'; lin[303] = '\0'; run_proc(NULL, lin, out, sizeof(out), &code); CHK(!parse_double_strictly("abc", &v18) && !parse_double_strictly("nan", &v18) && !parse_double_strictly("+inf", &v18) && !parse_double_strictly("1e309", &v18) && parse_double_strictly("-42.5e-3", &v18) && parse_long_strictly("123", &l18) && strstr(out, "exceeds buffer limit"), 18, "Strict input parsing (NaN, inf, range, buffer limits)");
    CHK(run_proc(NULL, "", out, sizeof(out), &code) && !code && strstr(out, "Exiting cleanly") && run_proc(NULL, "1\n", out, sizeof(out), &code) && !code && strstr(out, "Reference point change canceled") && run_proc(NULL, "2\n", out, sizeof(out), &code) && !code && strstr(out, "Background field update canceled") && run_proc(NULL, "3\n", out, sizeof(out), &code) && !code && strstr(out, "Output 4 canceled"), 19, "EOF handling verified at all prompts without hanging");
    run_proc(NULL, "1\n1e308\n1e308\n0\n", out, sizeof(out), &code); CHK(!strstr(out, "nan") && !strstr(out, "NaN") && !strstr(out, "Infinity"), 20, "Calculation-range failure guarded safely; no NaN/infinity exposed");
    run_proc("--preset test01", "1\n2.0\n0.0\n1\n1.0\n0.0\n0\n", out, sizeof(out), &code); CHK(strstr(out, "8.982733") != NULL, 21, "Display is not state (internal calculations retain double precision)");
    run_proc(NULL, "0\n", out, sizeof(out), &code); CHK(!strstr(out, "DEBUG") && !strstr(out, "TODO") && strstr(out, "microC") && strstr(out, "N/C") && strstr(out, "degrees"), 22, "Output discipline verified (clean formatting, proper units)");
    FILE *f1 = fopen("run_all_acceptance_checks.py", "r"), *f2 = fopen("test_interactive.py", "r"); bool no_py = (!f1 && !f2); if (f1) fclose(f1); if (f2) fclose(f2); CHK(no_py, 23, "C-only deliverables verified (all test/utility code is C)");
    printf("===============================================================\nACCEPTANCE CHECKS SUMMARY: %d PASSED, %d FAILED\n===============================================================\n", passed, failed);
    return (failed == 0) ? 0 : 1;
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        if (!strcmp(argv[1], "--acceptance-checks") || !strcmp(argv[1], "--test")) return run_internal_acceptance_checks();
        if (!strcmp(argv[1], "--preset") && argc > 2) return run_preset_mode(argv[2]);
        printf("Usage: %s [--acceptance-checks | --preset <name>]\n", argv[0]); return 1;
    }
    return run_demo_mode();
}

#endif /* MODULE2_NO_MAIN */
