/*
 * ============================================================================
 * Central Integration Controller: main.c
 * Electrostatic Simulation Suite (IC103 Programming Project)
 *
 * Coordinates All Five Modules:
 *   Module 1: Charge Configuration Manager (charge_config.h/.c)
 *   Module 2: Electric Field Calculator     (electric_field.h/.c)
 *   Module 3: Equipotential Visualizer      (equipotential.h/.c)
 *   Module 4: 2D Electric Field Mapper      (electric_field_mapper.h/.c)
 *   Module 5: Results & Report Manager      (results_manager.h/.c)
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Module Header Integrations */
#include "charge_config.h"
#include "electric_field.h"
#include "equipotential.h"
#include "electric_field_mapper.h"
#include "results_manager.h"

/* Helper to flush leftover characters from stdin */
static void flush_stdin(void)
{
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {
        /* discard */
    }
}

/*
 * Synchronize charge configuration from Module 1 into Module 3 (Equipotential).
 */
static void sync_config_to_equipotential(const ChargeConfig *config)
{
    if (config == NULL || config->n <= 0) {
        chargeCount = 0;
        return;
    }

    chargeCount = config->n;
    for (int i = 0; i < config->n; i++) {
        charges[i] = (float)config->q[i];
    }
}

/*
 * Synchronize charge configuration from Module 1 into Module 4 (Field Mapper).
 */
static void sync_config_to_mapper(const ChargeConfig *config, Charge *mapper_charges, int *charge_count)
{
    if (config == NULL || config->n <= 0) {
        *charge_count = 0;
        return;
    }

    *charge_count = config->n;
    W = (config->gridSize > 0 && config->gridSize <= MAXG) ? config->gridSize : 25;
    H = (config->gridSize > 0 && config->gridSize <= MAXG) ? config->gridSize : 25;

    for (int i = 0; i < config->n; i++) {
        mapper_charges[i].q = config->q[i];
        mapper_charges[i].x = config->x[i];
        mapper_charges[i].y = config->y[i];
    }
}

/*
 * Synchronize simulation state across Modules 1, 2, and 4 into Module 5 (ReportData).
 */
static void sync_state_to_results(const ChargeConfig *config,
                                  int single_point_computed,
                                  int grid_computed,
                                  ReportData *report_data)
{
    if (config == NULL || report_data == NULL) {
        return;
    }

    results_init(report_data);

    /* Generate simulation identifier and title */
    time_t now = time(NULL);
    snprintf(report_data->simulation_id, sizeof(report_data->simulation_id),
             "SIM_%08lX", (unsigned long)now);
    snprintf(report_data->configuration_name, sizeof(report_data->configuration_name),
             "Electrostatic System (%d Charges, %dx%d Grid)",
             config->n, config->gridSize, config->gridSize);

    /* Ingest charge configuration from Module 1 */
    report_data->number_of_charges = config->n;
    for (int i = 0; i < config->n; i++) {
        report_data->charges[i].q = config->q[i];
        report_data->charges[i].x = config->x[i];
        report_data->charges[i].y = config->y[i];
    }

    /* Ingest single-point field results from Module 2 */
    if (single_point_computed) {
        report_data->has_single_point = 1;
        report_data->reference_x = referenceX;
        report_data->reference_y = referenceY;
        report_data->charge_Ex = resultantField.Ex;
        report_data->charge_Ey = resultantField.Ey;
        report_data->charge_field_magnitude = resultantField.magnitude;
        report_data->charge_field_angle = resultantField.angle;
        report_data->charge_field_defined = resultantField.defined;

        report_data->external_field_active = externalField.defined;
        report_data->external_Ex = externalField.Ex;
        report_data->external_Ey = externalField.Ey;

        report_data->net_Ex = netField.Ex;
        report_data->net_Ey = netField.Ey;
        report_data->net_field_magnitude = netField.magnitude;
        report_data->net_field_angle = netField.angle;
        report_data->net_field_defined = netField.defined;
    } else {
        report_data->has_single_point = 0;
    }

    /* Ingest grid calculation results from Module 4 and calculate potentials */
    if (grid_computed) {
        report_data->has_grid = 1;
        report_data->grid_width = (W <= REPORT_MAX_GRID) ? W : REPORT_MAX_GRID;
        report_data->grid_height = (H <= REPORT_MAX_GRID) ? H : REPORT_MAX_GRID;

        double k_coulomb = calculate_coulomb_constant();

        for (int j = 0; j < report_data->grid_height; j++) {
            for (int i = 0; i < report_data->grid_width; i++) {
                report_data->grid[j][i].x = (double)i;
                report_data->grid[j][i].y = (double)j;
                report_data->grid[j][i].Ex = Ex[j][i];
                report_data->grid[j][i].Ey = Ey[j][i];
                report_data->grid[j][i].magnitude = E[j][i];
                report_data->grid[j][i].angle = ang[j][i];
                report_data->grid[j][i].field_valid = !bad[j][i];

                if (bad[j][i]) {
                    report_data->grid[j][i].potential = 0.0;
                    report_data->grid[j][i].potential_valid = 0;
                } else {
                    double v_total = 0.0;
                    for (int c_idx = 0; c_idx < config->n; c_idx++) {
                        double dx = (double)i - config->x[c_idx];
                        double dy = (double)j - config->y[c_idx];
                        double r = hypot(dx, dy);
                        if (r > 1e-9) {
                            v_total += (k_coulomb * config->q[c_idx]) / r;
                        }
                    }
                    report_data->grid[j][i].potential = v_total;
                    report_data->grid[j][i].potential_valid = 1;
                }
            }
        }
    } else {
        report_data->has_grid = 0;
    }
}

/*
 * Automated End-to-End Pipeline:
 * Coordinates all five modules in logical execution sequence.
 */
static void run_complete_pipeline(ChargeConfig *config,
                                  int *is_configured,
                                  int *single_point_computed,
                                  int *grid_computed,
                                  ReportData *report_data)
{
    printf("\n=============================================================\n");
    printf("   STARTING END-TO-END ELECTROSTATIC SIMULATION PIPELINE     \n");
    printf("=============================================================\n");

    /* Step 1: Charge Configuration */
    if (!(*is_configured)) {
        printf("\n[Step 1/5] Ingesting Charge Configuration (Module 1)...\n");
        if (!readChargeConfiguration(config)) {
            printf("Pipeline aborted: Failed to configure charges.\n");
            return;
        }
        *is_configured = 1;
    } else {
        printf("\n[Step 1/5] Using existing Charge Configuration (Module 1):\n");
        displayChargeConfiguration(config);
    }

    /* Step 2: Electric Field Calculation */
    printf("\n[Step 2/5] Electric Field Calculation (Module 2)...\n");
    double rx = 0.0, ry = 0.0;
    printf("Enter observation / reference point X coordinate (m): ");
    if (scanf("%lf", &rx) != 1) rx = 0.0;
    printf("Enter observation / reference point Y coordinate (m): ");
    if (scanf("%lf", &ry) != 1) ry = 0.0;
    flush_stdin();

    clear_external_field();
    calculate_electric_field(config, rx, ry);
    display_results();
    *single_point_computed = 1;

    /* Step 3: Equipotential Simulation */
    printf("\n[Step 3/5] Equipotential Simulation (Module 3)...\n");
    sync_config_to_equipotential(config);

    char launch_eq = 'y';
    printf("Launch graphical Equipotential window (Module 3)? (y/n) [default: y]: ");
    int ch = getchar();
    if (ch != '\n' && ch != EOF) {
        launch_eq = (char)ch;
        flush_stdin();
    }
    if (launch_eq == 'y' || launch_eq == 'Y') {
        printf("Opening Equipotential window (press ESC or close window to return)...\n");
        RunEquipotentialSimulation();
        printf("Equipotential window closed. Control returned to pipeline.\n");
    } else {
        printf("Equipotential graphical window skipped by user.\n");
    }

    /* Step 4: 2D Electric Field Mapping */
    printf("\n[Step 4/5] 2D Electric Field Mapping (Module 4)...\n");
    Charge mapper_charges[MAX_CHARGES];
    int mapper_count = 0;
    sync_config_to_mapper(config, mapper_charges, &mapper_count);

    compute_field(mapper_charges, mapper_count);
    compute_mag_angle();
    print_map(mapper_charges, mapper_count);
    save_svg(mapper_charges, mapper_count);
    *grid_computed = 1;

    /* Step 5: Results Consolidation and Report Management */
    printf("\n[Step 5/5] Consolidating Simulation Results (Module 5)...\n");
    sync_state_to_results(config, *single_point_computed, *grid_computed, report_data);
    displaySimulationResults(report_data);

    /* Open Results Manager Interactive Menu */
    printf("\nEntering Results and Report Manager (Module 5)...\n");
    runResultsManager(report_data);

    printf("\n=============================================================\n");
    printf("     END-TO-END SIMULATION PIPELINE COMPLETED SUCCESSFULLY   \n");
    printf("=============================================================\n");
}

/*
 * Main Application Entry Point
 */
int main(void)
{
    ChargeConfig config;
    memset(&config, 0, sizeof(config));

    ReportData *report_data = (ReportData *)calloc(1, sizeof(ReportData));
    if (report_data == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for ReportData.\n");
        return 1;
    }

    int is_configured = 0;
    int single_point_computed = 0;
    int grid_computed = 0;
    int choice = -1;

    printf("=============================================================\n");
    printf("    WELCOME TO THE ELECTROSTATIC SIMULATION SUITE (IC103)   \n");
    printf("=============================================================\n");

    while (1) {
        printf("\n======================= MAIN MENU ===========================\n");
        printf("Status: %s", is_configured ? "CHARGES CONFIGURED" : "NO CHARGES CONFIGURED");
        if (is_configured) {
            printf(" (%d charges, %dx%d grid)\n", config.n, config.gridSize, config.gridSize);
        } else {
            printf("\n");
        }
        printf("-------------------------------------------------------------\n");
        printf("1. Run Complete End-to-End Simulation Pipeline (Modules 1-5)\n");
        printf("2. Module 1: Configure Charges (Charge Configuration Manager)\n");
        printf("3. Module 2: Electric Field Calculator (Point & Superposition)\n");
        printf("4. Module 3: Equipotential Simulation (Raylib Visualizer)\n");
        printf("5. Module 4: 2D Electric Field Mapper (Terminal Map & SVG)\n");
        printf("6. Module 5: Results & Report Manager (View, TXT/CSV Export)\n");
        printf("7. View Current Configuration & System Status\n");
        printf("0. Exit Application\n");
        printf("=============================================================\n");
        printf("Enter your choice (0-7): ");

        if (scanf("%d", &choice) != 1) {
            flush_stdin();
            printf("Invalid input. Please enter a number between 0 and 7.\n");
            continue;
        }
        flush_stdin();

        switch (choice) {
            case 1:
                run_complete_pipeline(&config,
                                      &is_configured,
                                      &single_point_computed,
                                      &grid_computed,
                                      report_data);
                break;

            case 2:
                printf("\n--- Module 1: Charge Configuration Manager ---\n");
                if (readChargeConfiguration(&config)) {
                    is_configured = 1;
                    displayChargeConfiguration(&config);
                    /* Reset dependent calculation flags */
                    single_point_computed = 0;
                    grid_computed = 0;
                } else {
                    printf("Configuration failed or was cancelled.\n");
                }
                break;

            case 3:
                printf("\n--- Module 2: Electric Field Calculator ---\n");
                if (!is_configured) {
                    printf("No charges configured yet. Please run Module 1 first.\n");
                    break;
                }
                run_electric_field_calculator(&config);
                single_point_computed = 1;
                break;

            case 4:
                printf("\n--- Module 3: Equipotential Simulation ---\n");
                if (!is_configured) {
                    printf("Warning: No charges configured in Module 1.\n");
                    printf("Module 3 will allow entering charges directly for standalone visualization.\n");
                    chargeCount = 0;
                } else {
                    sync_config_to_equipotential(&config);
                }
                printf("Launching Equipotential window (press ESC or close window to return)...\n");
                RunEquipotentialSimulation();
                printf("Equipotential simulation finished.\n");
                break;

            case 5:
                printf("\n--- Module 4: 2D Electric Field Mapper ---\n");
                if (!is_configured) {
                    printf("No charges configured yet. Please run Module 1 first.\n");
                    break;
                }
                {
                    Charge mapper_charges[MAX_CHARGES];
                    int mapper_count = 0;
                    sync_config_to_mapper(&config, mapper_charges, &mapper_count);

                    printf("Calculating 2D vector field on %dx%d grid...\n", W, H);
                    compute_field(mapper_charges, mapper_count);
                    compute_mag_angle();
                    print_map(mapper_charges, mapper_count);
                    save_svg(mapper_charges, mapper_count);
                    grid_computed = 1;
                }
                break;

            case 6:
                printf("\n--- Module 5: Results & Report Manager ---\n");
                if (!is_configured) {
                    printf("No simulation data available. Please configure charges first.\n");
                    break;
                }
                sync_state_to_results(&config, single_point_computed, grid_computed, report_data);
                runResultsManager(report_data);
                break;

            case 7:
                printf("\n--- System Status & Active Configuration ---\n");
                if (!is_configured) {
                    printf("Charges: Not configured.\n");
                } else {
                    displayChargeConfiguration(&config);
                    printf("\nCalculation Status:\n");
                    printf("  - Module 2 Single-Point Field: %s\n",
                           single_point_computed ? "Calculated" : "Not yet run");
                    if (single_point_computed) {
                        printf("    Reference Point: (%.4f, %.4f) m | Net |E|: %.4e N/C\n",
                               referenceX, referenceY, netField.magnitude);
                    }
                    printf("  - Module 4 2D Field Grid:      %s\n",
                           grid_computed ? "Calculated" : "Not yet run");
                    if (grid_computed) {
                        printf("    Grid Dimensions: %d x %d | Maximum |E|: %.4e N/C\n",
                               W, H, Emax);
                    }
                }
                break;

            case 0:
                printf("\nCleaning up resources and exiting...\n");
                if (report_data != NULL) {
                    free(report_data);
                    report_data = NULL;
                }
                printf("Thank you for using the Electrostatic Simulation Suite. Goodbye!\n");
                return 0;

            default:
                printf("Invalid choice. Please select an option from 0 to 7.\n");
                break;
        }
    }

    /* Fallback cleanup */
    if (report_data != NULL) {
        free(report_data);
        report_data = NULL;
    }
    return 0;
}
