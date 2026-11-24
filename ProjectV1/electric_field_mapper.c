/* Module 4: 2D Electric Field Mapper
 * Electric Field & Potential Simulator
 *
 * Calculates the electric field of several point charges on a 2D grid
 * (Coulomb's law + superposition), prints an arrow map in the terminal
 * and exports an SVG picture.
 */
#include "electric_field_mapper.h"

/* ---------- Global variable DEFINITIONS (exactly once, no 'extern') ---------- */
int    W, H;
double Ex[MAXG][MAXG];
double Ey[MAXG][MAXG];
double E[MAXG][MAXG];
double ang[MAXG][MAXG];
int    bad[MAXG][MAXG];
double Emax = 0.0;

/* Points closer than this to a charge are skipped (field would blow up) */
#define MIN_DIST 0.5

/* ---------- Step 1: Coulomb's law + superposition ---------- */
void compute_field(Charge c[], int n)
{
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            double ex = 0.0, ey = 0.0;
            bad[j][i] = 0;

            for (int k = 0; k < n; k++) {
                double dx = i - c[k].x;
                double dy = j - c[k].y;
                double r2 = dx * dx + dy * dy;
                double r  = sqrt(r2);

                if (r < MIN_DIST) {      /* too close to a charge */
                    bad[j][i] = 1;
                    break;
                }
                double e = K * c[k].q / r2;   /* |E| from this charge     */
                ex += e * dx / r;             /* superposition: add parts */
                ey += e * dy / r;
            }
            Ex[j][i] = bad[j][i] ? 0.0 : ex;
            Ey[j][i] = bad[j][i] ? 0.0 : ey;
        }
    }
}

/* ---------- Step 2: magnitude, direction (degrees), maximum ---------- */
void compute_mag_angle(void)
{
    Emax = 0.0;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            if (bad[j][i]) {
                E[j][i] = 0.0;
                ang[j][i] = 0.0;
                continue;
            }
            E[j][i]   = sqrt(Ex[j][i] * Ex[j][i] + Ey[j][i] * Ey[j][i]);
            ang[j][i] = atan2(Ey[j][i], Ex[j][i]) * 180.0 / PI;
            if (ang[j][i] < 0) ang[j][i] += 360.0;
            if (E[j][i] > Emax) Emax = E[j][i];
        }
    }
}

/* Pick an arrow for an angle (8 directions) */
static const char *arrow_for(double a)
{
    static const char *arrows[8] = {
        "\xE2\x86\x92", /* -> */  "\xE2\x86\x97", /* NE */
        "\xE2\x86\x91", /* ^  */  "\xE2\x86\x96", /* NW */
        "\xE2\x86\x90", /* <- */  "\xE2\x86\x99", /* SW */
        "\xE2\x86\x93", /* v  */  "\xE2\x86\x98"  /* SE */
    };
    int idx = (int)((a + 22.5) / 45.0) % 8;
    return arrows[idx];
}

/* Pick an ANSI color by field strength (relative to Emax; sqrt scale) */
static const char *color_for(double e)
{
    double t = (Emax > 0) ? sqrt(e / Emax) : 0.0;
    if (t > 0.75) return "\x1b[31m";   /* red     = strong */
    if (t > 0.50) return "\x1b[33m";   /* yellow         */
    if (t > 0.25) return "\x1b[32m";   /* green          */
    return "\x1b[34m";                 /* blue    = weak   */
}

/* ---------- Step 3: terminal map ---------- */
void print_map(Charge c[], int n)
{
    printf("\nElectric field map (red = strong, blue = weak)\n");
    printf("+ positive charge, - negative charge\n\n");

    for (int j = H - 1; j >= 0; j--) {          /* y grows upward */
        for (int i = 0; i < W; i++) {
            int shown = 0;
            for (int k = 0; k < n; k++) {
                if ((int)(c[k].x + 0.5) == i && (int)(c[k].y + 0.5) == j) {
                    printf("%s%c\x1b[0m ", c[k].q > 0 ? "\x1b[1;31m" : "\x1b[1;34m",
                           c[k].q > 0 ? '+' : '-');
                    shown = 1;
                    break;
                }
            }
            if (shown) continue;
            if (bad[j][i]) { printf("  "); continue; }
            printf("%s%s\x1b[0m ", color_for(E[j][i]), arrow_for(ang[j][i]));
        }
        printf("\n");
    }
    printf("\nEmax = %.3e N/C\n", Emax);
}

/* ---------- Step 4: SVG export ---------- */
void save_svg(Charge c[], int n)
{
    const int CELL = 20;                       /* pixels per grid cell */
    FILE *f = fopen("field_map.svg", "w");
    if (!f) {
        printf("Could not create field_map.svg\n");
        return;
    }

    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n",
            W * CELL, H * CELL);
    fprintf(f, "<rect width=\"100%%\" height=\"100%%\" fill=\"white\"/>\n");

    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            if (bad[j][i] || E[j][i] <= 0) continue;

            double t   = sqrt(E[j][i] / Emax);          /* 0..1 strength */
            double len = 4.0 + 6.0 * t;                  /* half length   */
            double cx  = i * CELL + CELL / 2.0;
            double cy  = (H - 1 - j) * CELL + CELL / 2.0; /* flip y      */
            double a   = ang[j][i] * PI / 180.0;
            double ux  = cos(a), uy = -sin(a);           /* SVG y points down */

            double x1 = cx - ux * len, y1 = cy - uy * len;
            double x2 = cx + ux * len, y2 = cy + uy * len;
            int r = (int)(255 * t), b = (int)(255 * (1 - t));

            fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" "
                       "stroke=\"rgb(%d,0,%d)\" stroke-width=\"1.5\"/>\n",
                    x1, y1, x2, y2, r, b);
            /* arrow head */
            double hl = 4.0;
            fprintf(f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f\" "
                       "fill=\"rgb(%d,0,%d)\"/>\n",
                    x2, y2,
                    x2 - ux * hl - uy * hl * 0.5, y2 - uy * hl + ux * hl * 0.5,
                    x2 - ux * hl + uy * hl * 0.5, y2 - uy * hl - ux * hl * 0.5,
                    r, b);
        }
    }

    for (int k = 0; k < n; k++) {
        double cx = c[k].x * CELL + CELL / 2.0;
        double cy = (H - 1 - c[k].y) * CELL + CELL / 2.0;
        fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"7\" fill=\"%s\"/>\n",
                cx, cy, c[k].q > 0 ? "red" : "blue");
    }

    fprintf(f, "</svg>\n");
    fclose(f);
    printf("Saved field_map.svg\n");
}

#ifndef INTEGRATED_BUILD
/* ---------- main ---------- */
int main(void)
{
    Charge c[10];
    int n;

    W = 30;
    H = 15;

    printf("Number of charges (1-10): ");
    if (scanf("%d", &n) != 1 || n < 1 || n > 10) {
        printf("Invalid number of charges.\n");
        return 1;
    }
    for (int k = 0; k < n; k++) {
        double q_uC;
        printf("Charge %d - q (microcoulombs), x (0-%d), y (0-%d): ", k + 1, W - 1, H - 1);
        if (scanf("%lf %lf %lf", &q_uC, &c[k].x, &c[k].y) != 3) {
            printf("Invalid input.\n");
            return 1;
        }
        c[k].q = q_uC * 1e-6;       /* microcoulombs -> coulombs */
    }

    compute_field(c, n);
    compute_mag_angle();
    print_map(c, n);
    save_svg(c, n);
    return 0;
}
#endif
