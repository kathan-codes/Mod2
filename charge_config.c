#include "charge_config.h"
#include <math.h>
#include <stdio.h>

/*
 * MODULE 1: CHARGE CONFIGURATION MANAGER
 *
 * The input sequence and original variable names from propro.c are retained:
 * n, q, x, y, and gridSize.
 *
 * The main addition is that the input is placed in a function and copied
 * into ChargeConfig, allowing other modules to receive the same data.
 */

int readChargeConfiguration(ChargeConfig *config) {
  int n, gridSize;
  int i;
  double q[MAX_CHARGES];
  double x[MAX_CHARGES], y[MAX_CHARGES];

  if (config == NULL) {
    printf("Invalid configuration pointer.\n");
    return 0;
  }

  /* Number of charges */
  printf("Enter number of charges: ");

  if (scanf("%d", &n) != 1) {
    printf("Invalid input for number of charges.\n");
    return 0;
  }

  if (n <= 0 || n > MAX_CHARGES) {
    printf("Number of charges must be between 1 and %d.\n", MAX_CHARGES);
    return 0;
  }

  printf("\nNOTE: Enter numbers only; units are already specified.\n");
  printf("Scientific notation is accepted (for example, 1e-9 = 1 x 10^-9).\n");

  /* Charge values and coordinates */
  for (i = 0; i < n; i++) {
    printf("\nCharge %d:\n", i + 1);

    printf("Enter charge in Coulombs (C) [e.g. 1e-9, -2e-6]: ");
    if (scanf("%lf", &q[i]) != 1 || !isfinite(q[i])) {
      printf("Invalid charge value. Enter a finite number.\n");
      return 0;
    }

    printf("Enter x-coordinate in metres (m) [e.g. 0.5, -1.2]: ");
    if (scanf("%lf", &x[i]) != 1 || !isfinite(x[i])) {
      printf("Invalid x coordinate. Enter a finite number.\n");
      return 0;
    }

    printf("Enter y-coordinate in metres (m) [e.g. 0.5, -1.2]: ");
    if (scanf("%lf", &y[i]) != 1 || !isfinite(y[i])) {
      printf("Invalid y coordinate. Enter a finite number.\n");
      return 0;
    }
  }

  /* Grid size */
  printf("\nEnter grid size (1 to %d): ", MAX_GRID_SIZE);

  if (scanf("%d", &gridSize) != 1) {
    printf("Invalid grid size input.\n");
    return 0;
  }

  if (gridSize < 1 || gridSize > MAX_GRID_SIZE) {
    printf("Grid size must be between 1 and %d.\n", MAX_GRID_SIZE);
    return 0;
  }

  /*
   * Copy the values into the shared structure.
   * The names n, q, x, y, and gridSize above are unchanged.
   */
  config->n = n;
  config->gridSize = gridSize;

  for (i = 0; i < n; i++) {
    config->q[i] = q[i];
    config->x[i] = x[i];
    config->y[i] = y[i];
  }

  return 1;
}

void displayChargeConfiguration(const ChargeConfig *config) {
  int i;

  if (config == NULL)
    return;

  printf("\n========== CHARGE CONFIGURATION ==========\n");
  printf("Number of charges: %d\n", config->n);

  for (i = 0; i < config->n; i++) {
    printf("Charge %d: q = %.6g C, x = %.6g m, y = %.6g m\n", i + 1,
           config->q[i], config->x[i], config->y[i]);
  }

  printf("Grid size: %d x %d\n", config->gridSize, config->gridSize);
  printf("==========================================\n");
}

#include "electric_field.h"

int main(void) {
  ChargeConfig config;
  if (readChargeConfiguration(&config)) {
    displayChargeConfiguration(&config);
    run_electric_field_calculator(&config);
  }
  return 0;
}
