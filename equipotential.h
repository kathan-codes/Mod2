#ifndef EQUIPOTENTIAL_H
#define EQUIPOTENTIAL_H

#include "raylib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Screen and simulation world constants */
#define SCREEN_W 1000
#define SCREEN_H 700
#ifndef MAX_CHARGES
#define MAX_CHARGES 100
#endif

#define WORLD_W 12.0f
#define WORLD_H 8.4f

/* Charge state variables */
extern float charges[MAX_CHARGES];
extern int chargeCount;

/* Function declarations */
Vector2 ChargePosition(int index, int columns, int rows);
Vector2 WorldToScreen(Vector2 p);
void RunEquipotentialSimulation(void);

#ifdef __cplusplus
}
#endif

#endif /* EQUIPOTENTIAL_H */
