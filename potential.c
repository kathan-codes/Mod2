#include <stdio.h>
#include <math.h>
#include "potential.h"

/* MODULE 3: ELECTRIC POTENTIAL CALCULATOR */

int calculatePotentialAtPoint(
    const ChargeConfig *config,
    double px,
    double py,
    double *potential)
{
    int i;
    double total = 0.0;

    if (config == NULL || potential == NULL ||
        config->n < 1 || config->n > MAX_CHARGES ||
        !isfinite(px) || !isfinite(py))
        return POTENTIAL_ERROR;

    for (i = 0; i < config->n; i++)
    {
        double dx, dy, distance, contribution;

        /* Validate before doing arithmetic with the charge data. */
        if (!isfinite(config->q[i]) ||
            !isfinite(config->x[i]) ||
            !isfinite(config->y[i]))
            return POTENTIAL_ERROR;

        dx = px - config->x[i];
        dy = py - config->y[i];
        distance = hypot(dx, dy);

        /* A point charge has undefined potential at its exact location. */
        if (distance == 0.0)
            return POTENTIAL_UNDEFINED;

        /* Overflow or invalid arithmetic is a numerical error, not undefined physics. */
        if (!isfinite(distance))
            return POTENTIAL_ERROR;

        contribution = COULOMB_CONSTANT * config->q[i] / distance;

        if (!isfinite(contribution) || !isfinite(total + contribution))
            return POTENTIAL_ERROR;

        total += contribution;
    }

    *potential = total;
    return POTENTIAL_SUCCESS;
}


/* MODULE 3: POTENTIAL ENERGY OF AN EXTERNAL TEST CHARGE */

int calculateTestChargePotentialEnergy(
    const ChargeConfig *config,
    double px,
    double py,
    double testCharge,
    double *potentialEnergy)
{
    double potential;
    double energy;
    int status;

    if (potentialEnergy == NULL || !isfinite(testCharge))
        return POTENTIAL_ERROR;

    /* Reuse the existing validated potential calculation. */
    status = calculatePotentialAtPoint(config, px, py, &potential);
    if (status != POTENTIAL_SUCCESS)
        return status;

    energy = testCharge * potential;
    if (!isfinite(energy))
        return POTENTIAL_ERROR;

    *potentialEnergy = energy;
    return POTENTIAL_SUCCESS;
}


/* MODULE 3: TOTAL INTERACTION ENERGY OF THE CHARGE CONFIGURATION */

int calculateTotalPotentialEnergy(
    const ChargeConfig *config,
    double *totalPotentialEnergy)
{
    int i, j;
    double total = 0.0;

    if (config == NULL || totalPotentialEnergy == NULL ||
        config->n < 1 || config->n > MAX_CHARGES)
        return POTENTIAL_ERROR;

    /* Validate all charge data before calculating pair interactions. */
    for (i = 0; i < config->n; i++)
    {
        if (!isfinite(config->q[i]) ||
            !isfinite(config->x[i]) ||
            !isfinite(config->y[i]))
            return POTENTIAL_ERROR;
    }

    /* Each distinct pair is included exactly once (j starts at i + 1). */
    for (i = 0; i < config->n - 1; i++)
    {
        for (j = i + 1; j < config->n; j++)
        {
            double dx = config->x[i] - config->x[j];
            double dy = config->y[i] - config->y[j];
            double distance = hypot(dx, dy);
            double pairEnergy;

            /* Coincident ideal point charges have undefined interaction energy. */
            if (distance == 0.0)
                return POTENTIAL_UNDEFINED;

            if (!isfinite(distance))
                return POTENTIAL_ERROR;

            pairEnergy = COULOMB_CONSTANT * config->q[i] * config->q[j]
                         / distance;

            if (!isfinite(pairEnergy) || !isfinite(total + pairEnergy))
                return POTENTIAL_ERROR;

            total += pairEnergy;
        }
    }

    *totalPotentialEnergy = total;
    return POTENTIAL_SUCCESS;
}
