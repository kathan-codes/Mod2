#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "charge_config.h"

#define COULOMB_CONSTANT 8.9875517923e9

/*
 * Return values used by calculatePotentialAtPoint():
 *  1 = success
 *  0 = point coincides with a point charge (potential is undefined)
 * -1 = invalid input or numerical overflow/error
 */
#define POTENTIAL_ERROR (-1)
#define POTENTIAL_UNDEFINED 0
#define POTENTIAL_SUCCESS 1

int calculatePotentialAtPoint(
    const ChargeConfig *config,
    double px,
    double py,
    double *potential
);

/*
 * Calculate the potential energy of an external test charge q0 placed
 * at (px, py): U = q0 * V(px, py).
 * Uses the same return values as calculatePotentialAtPoint().
 */
int calculateTestChargePotentialEnergy(
    const ChargeConfig *config,
    double px,
    double py,
    double testCharge,
    double *potentialEnergy
);

/*
 * Alias for calculateTestChargePotentialEnergy.
 * Calculates potential energy at point (px, py) for a given test charge.
 */
int calculatePotentialEnergyAtPoint(
    const ChargeConfig *config,
    double px,
    double py,
    double testCharge,
    double *potentialEnergy
);

/*
 * Calculate the total electrostatic interaction energy of the configured
 * point charges, counting each pair once:
 * U = sum(i < j) [k * qi * qj / rij].
 * Uses the same return values as calculatePotentialAtPoint().
 */
int calculateTotalPotentialEnergy(
    const ChargeConfig *config,
    double *totalPotentialEnergy
);

/*
 * Alias for calculateTotalPotentialEnergy.
 * Calculates the total electrostatic potential energy of the system of charges.
 */
int calculateSystemPotentialEnergy(
    const ChargeConfig *config,
    double *systemPotentialEnergy
);

#endif
