#ifndef SIMPLE_FIELD_MAP_H
#define SIMPLE_FIELD_MAP_H

/* Module 4: 2D Electric Field Mapper
 * Header file: declarations only (no function bodies, no variable definitions).
 * NOTE: values/types marked GUESS must be checked against simple_field_map.c
 */

/* ---------- Standard headers ---------- */
#include <stdio.h>   /* printf, FILE, fopen, fprintf, fclose */
#include <math.h>    /* sqrt, atan2 (link with -lm)          */

/* ---------- Constants ---------- */
#define MAXG 100              /* maximum grid size                     */
#define K    8.98271086e9     /* Coulomb constant for dry air (N*m^2/C^2) */
#define PI   3.14159265358979323846

/* ---------- Data type ---------- */
typedef struct {
    double q;   /* charge value  (GUESS: double) */
    double x;   /* x-coordinate  (GUESS: double) */
    double y;   /* y-coordinate  (GUESS: double) */
} Charge;

/* ---------- Global variables (defined once in simple_field_map.c) ---------- */
extern int    W, H;               /* grid width and height                   */
extern double Ex[MAXG][MAXG];     /* GUESS: electric field x-component       */
extern double Ey[MAXG][MAXG];     /* GUESS: electric field y-component       */
extern double E[MAXG][MAXG];      /* GUESS: field magnitude                  */
extern double ang[MAXG][MAXG];    /* GUESS: field direction (degrees)        */
extern int    bad[MAXG][MAXG];    /* GUESS: invalid grid points              */
extern double Emax;               /* maximum field strength                  */

/* ---------- Function prototypes ---------- */
void compute_field(Charge c[], int n);   /* GUESS: void return type */
void compute_mag_angle(void);
void print_map(Charge c[], int n);
void save_svg(Charge c[], int n);

#endif /* SIMPLE_FIELD_MAP_H */
