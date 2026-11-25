/*
 * ============================================================================
 * Central Integration Entry Point: main.c
 * Electrostatic Simulation Suite (Integrated Project)
 *
 * Coordinates All Six Modules:
 *   Module 1: Charge Configuration Manager (charge_config.h / charge_config.c)
 *   Module 2: Electric Field Calculator     (electric_field.h / electric_field.c)
 *   Module 3: Electric Potential Calculator (potential.h / potential.c)
 *   Module 4: 2D Electric Field Mapper      (electric_field_mapper.h / electric_field_mapper.c)
 *   Module 5: Equipotential Visualizer      (equipotential.h / equipotential.c)
 *   Module 6: Results & Report Manager      (results_manager.h / results_manager.c)
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#ifdef _WIN32
#ifndef CP_UTF8
#define CP_UTF8 65001
#endif
unsigned int __stdcall SetConsoleOutputCP(unsigned int wCodePageID);
unsigned int __stdcall SetConsoleCP(unsigned int wCodePageID);
#endif

/* Include All Six Module Header Files */
#include "charge_config.h"
#include "electric_field.h"
#include "potential.h"
#include "electric_field_mapper.h"
#include "equipotential.h"
#include "results_manager.h"

/* Helper function: flush leftover characters from stdin */
static void flush_stdin(void)
{
    int c;
    while ((c = getchar()) != '\n' && c != EOF)
    {
        /* discard */
    }
}

/*
 * Synchronize charge configuration from Module 1 into Module 4 (Field Mapper).
 */
static void sync_config_to_mapper(const ChargeConfig *config, Charge *mapper_charges, int *charge_count)
{
    if (config == NULL || config->n <= 0)
    {
        *charge_count = 0;
        return;
    }

    *charge_count = config->n;
    W = (config->gridSize > 0 && config->gridSize <= MAXG) ? config->gridSize : 25;
    H = (config->gridSize > 0 && config->gridSize <= MAXG) ? config->gridSize : 25;

    for (int i = 0; i < config->n; i++)
    {
        mapper_charges[i].q = config->q[i];
        mapper_charges[i].x = config->x[i];
        mapper_charges[i].y = config->y[i];
    }
}

/*
 * Synchronize charge configuration from Module 1 into Module 5 (Equipotential).
 */
static void sync_config_to_equipotential(const ChargeConfig *config)
{
    if (config == NULL || config->n <= 0)
    {
        chargeCount = 0;
        return;
    }

    chargeCount = config->n;
    for (int i = 0; i < config->n; i++)
    {
        charges[i] = (float)config->q[i];
    }
}

/*
 * Synchronize simulation state across Modules 1, 2, 3, and 4 into Module 6 (ReportData).
 */
static void sync_state_to_results(const ChargeConfig *config,
                                  int single_point_computed,
                                  int grid_computed,
                                  ReportData *report_data)
{
    if (config == NULL || report_data == NULL)
    {
        return;
    }

    results_init(report_data);

    /* Generate simulation identifier and timestamped title */
    time_t now = time(NULL);
    snprintf(report_data->simulation_id, sizeof(report_data->simulation_id),
             "SIM_%08lX", (unsigned long)now);
    snprintf(report_data->configuration_name, sizeof(report_data->configuration_name),
             "Electrostatic System (%d Charges, %dx%d Grid)",
             config->n, config->gridSize, config->gridSize);

    /* Ingest charge configuration from Module 1 */
    report_data->number_of_charges = config->n;
    for (int i = 0; i < config->n; i++)
    {
        report_data->charges[i].q = config->q[i];
        report_data->charges[i].x = config->x[i];
        report_data->charges[i].y = config->y[i];
    }

    /* Ingest single-point field results from Module 2 */
    if (single_point_computed)
    {
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
    }
    else
    {
        report_data->has_single_point = 0;
    }

    /* Ingest grid calculation results from Module 4 and calculate potentials via Module 3 */
    if (grid_computed)
    {
        report_data->has_grid = 1;
        report_data->grid_width = (W <= REPORT_MAX_GRID) ? W : REPORT_MAX_GRID;
        report_data->grid_height = (H <= REPORT_MAX_GRID) ? H : REPORT_MAX_GRID;

        for (int j = 0; j < report_data->grid_height; j++)
        {
            for (int i = 0; i < report_data->grid_width; i++)
            {
                report_data->grid[j][i].x = (double)i;
                report_data->grid[j][i].y = (double)j;
                report_data->grid[j][i].Ex = Ex[j][i];
                report_data->grid[j][i].Ey = Ey[j][i];
                report_data->grid[j][i].magnitude = E[j][i];
                report_data->grid[j][i].angle = ang[j][i];
                report_data->grid[j][i].field_valid = !bad[j][i];

                /* Calculate potential using Module 3 */
                double v_val = 0.0;
                int pot_status = calculatePotentialAtPoint(config, (double)i, (double)j, &v_val);
                if (pot_status == POTENTIAL_SUCCESS)
                {
                    report_data->grid[j][i].potential = v_val;
                    report_data->grid[j][i].potential_valid = 1;
                }
                else
                {
                    report_data->grid[j][i].potential = 0.0;
                    report_data->grid[j][i].potential_valid = 0;
                }
            }
        }
    }
    else
    {
        report_data->has_grid = 0;
    }
}

/*
 * Automated End-to-End Pipeline:
 * Coordinates all six modules in logical execution sequence.
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

    /* Step 1: Charge Configuration (Module 1) */
    if (!(*is_configured))
    {
        printf("\n[Step 1/6] Ingesting Charge Configuration (Module 1)...\n");
        if (!readChargeConfiguration(config))
        {
            printf("Pipeline aborted: Failed to configure charges.\n");
            return;
        }
        *is_configured = 1;
    }
    else
    {
        printf("\n[Step 1/6] Using existing Charge Configuration (Module 1):\n");
        displayChargeConfiguration(config);
    }

    /* Step 2: Electric Field Calculation (Module 2) */
    printf("\n[Step 2/6] Electric Field Calculation (Module 2)...\n");
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

    /* Step 3: Electric Potential Calculation (Module 3) */
    printf("\n[Step 3/6] Electric Potential Calculation (Module 3)...\n");
    double point_pot = 0.0;
    int pot_status = calculatePotentialAtPoint(config, rx, ry, &point_pot);
    if (pot_status == POTENTIAL_SUCCESS)
    {
        printf("Electric Potential at (%.4f, %.4f) m: %.6e V\n", rx, ry, point_pot);
    }
    else if (pot_status == POTENTIAL_UNDEFINED)
    {
        printf("Electric Potential at (%.4f, %.4f) m: Undefined (coincides with charge)\n", rx, ry);
    }
    else
    {
        printf("Electric Potential calculation error at (%.4f, %.4f) m.\n", rx, ry);
    }

    /* Step 4: 2D Electric Field Mapping (Module 4) */
    printf("\n[Step 4/6] 2D Electric Field Mapping (Module 4)...\n");
    Charge mapper_charges[MAX_CHARGES];
    int mapper_count = 0;
    sync_config_to_mapper(config, mapper_charges, &mapper_count);

    compute_field(mapper_charges, mapper_count);
    compute_mag_angle();
    print_map(mapper_charges, mapper_count);
    save_svg(mapper_charges, mapper_count);
    *grid_computed = 1;

    /* Step 5: Equipotential Simulation (Module 5) */
    printf("\n[Step 5/6] Equipotential Simulation (Module 5)...\n");
    sync_config_to_equipotential(config);

    char launch_eq = 'y';
    printf("Launch graphical Equipotential window (Module 5)? (y/n) [default: y]: ");
    int ch = getchar();
    if (ch != '\n' && ch != EOF)
    {
        launch_eq = (char)ch;
        flush_stdin();
    }
    if (launch_eq == 'y' || launch_eq == 'Y')
    {
        printf("Opening Equipotential window (press ESC or close window to return)...\n");
        RunEquipotentialSimulation();
        printf("Equipotential window closed. Control returned to pipeline.\n");
    }
    else
    {
        printf("Equipotential graphical window skipped by user.\n");
    }

    /* Step 6: Results Consolidation and Report Management (Module 6) */
    printf("\n[Step 6/6] Consolidating Simulation Results (Module 6)...\n");
    sync_state_to_results(config, *single_point_computed, *grid_computed, report_data);
    displaySimulationResults(report_data);
    displayStatistics(report_data);

    /* Automatically export demonstration reports */
    exportTXT(report_data, "complete_demo_report.txt");
    exportCSV(report_data, "complete_demo_grid.csv");

    printf("\nEntering Results and Report Manager (Module 6)...\n");
    runResultsManager(report_data);

    printf("\n=============================================================\n");
    printf("     END-TO-END SIMULATION PIPELINE COMPLETED SUCCESSFULLY   \n");
    printf("=============================================================\n");
}

/*
 * Central Main Entry Point
 */
int main(void)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    ChargeConfig config;
    memset(&config, 0, sizeof(config));

    ReportData *report_data = (ReportData *)calloc(1, sizeof(ReportData));
    if (report_data == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed for ReportData.\n");
        return 1;
    }

    int is_configured = 0;
    int single_point_computed = 0;
    int grid_computed = 0;
    int choice = -1;

    printf("=============================================================\n");
    printf("    WELCOME TO THE ELECTROSTATIC SIMULATION SUITE (IC103)    \n");
    printf("    Integrating All Six Project Modules into One Platform    \n");
    printf("=============================================================\n");

    while (1)
    {
        printf("\n======================= MAIN MENU ===========================\n");
        printf("Status: %s", is_configured ? "CHARGES CONFIGURED" : "NO CHARGES CONFIGURED");
        if (is_configured)
        {
            printf(" (%d charges, %dx%d grid)\n", config.n, config.gridSize, config.gridSize);
        }
        else
        {
            printf("\n");
        }
        printf("-------------------------------------------------------------\n");
        printf("1. Charge Configuration          (Module 1)\n");
        printf("2. Electric Field Calculation    (Module 2)\n");
        printf("3. Electric Potential Calculation (Module 3)\n");
        printf("4. 2D Electric Field Mapping     (Module 4)\n");
        printf("5. Equipotential Mapping         (Module 5)\n");
        printf("6. Results Management            (Module 6)\n");
        printf("7. Run Complete Demonstration    (All 6 Modules in Sequence)\n");
        printf("8. Exit Application\n");
        printf("=============================================================\n");
        printf("Enter your choice (1-8): ");

        if (scanf("%d", &choice) != 1)
        {
            if (feof(stdin))
            {
                break;
            }
            flush_stdin();
            printf("Invalid input. Please enter a number between 1 and 8.\n");
            continue;
        }
        flush_stdin();

        switch (choice)
        {
            case 1:
                printf("\n--- Option 1: Charge Configuration (Module 1) ---\n");
                if (readChargeConfiguration(&config))
                {
                    is_configured = 1;
                    displayChargeConfiguration(&config);
                    single_point_computed = 0;
                    grid_computed = 0;
                }
                else
                {
                    printf("Configuration failed or was cancelled.\n");
                }
                break;

            case 2:
                printf("\n--- Option 2: Electric Field Calculation (Module 2) ---\n");
                if (!is_configured)
                {
                    printf("No charges configured yet. Please configure charges first (Option 1).\n");
                    break;
                }
                run_electric_field_calculator(&config);
                single_point_computed = 1;
                break;

            case 3:
                printf("\n--- Option 3: Electric Potential Calculation (Module 3) ---\n");
                if (!is_configured)
                {
                    printf("No charges configured yet. Please configure charges first (Option 1).\n");
                    break;
                }
                {
                    double px = 0.0, py = 0.0, val = 0.0;
                    printf("Enter point X coordinate (m): ");
                    if (scanf("%lf", &px) != 1) px = 0.0;
                    printf("Enter point Y coordinate (m): ");
                    if (scanf("%lf", &py) != 1) py = 0.0;
                    flush_stdin();

                    int st = calculatePotentialAtPoint(&config, px, py, &val);
                    if (st == POTENTIAL_SUCCESS)
                    {
                        printf("\n>>> Potential at (%.4f, %.4f) m = %.6e V\n", px, py, val);
                    }
                    else if (st == POTENTIAL_UNDEFINED)
                    {
                        printf("\n>>> Point coincides with a charge. Potential is mathematically undefined.\n");
                    }
                    else
                    {
                        printf("\n>>> Calculation error or invalid point.\n");
                    }
                }
                break;

            case 4:
                printf("\n--- Option 4: 2D Electric Field Mapping (Module 4) ---\n");
                if (!is_configured)
                {
                    printf("No charges configured yet. Please configure charges first (Option 1).\n");
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

            case 5:
                printf("\n--- Option 5: Equipotential Mapping (Module 5) ---\n");
                if (!is_configured)
                {
                    printf("Warning: No charges configured in Module 1.\n");
                    printf("Equipotential module will prompt for standalone charge values.\n");
                    chargeCount = 0;
                }
                else
                {
                    sync_config_to_equipotential(&config);
                }
                printf("Launching Equipotential window (press ESC or close window to return)...\n");
                RunEquipotentialSimulation();
                printf("Equipotential simulation finished.\n");
                break;

            case 6:
                printf("\n--- Option 6: Results Management (Module 6) ---\n");
                if (!is_configured)
                {
                    printf("No simulation data available. Please configure charges first (Option 1).\n");
                    break;
                }
                sync_state_to_results(&config, single_point_computed, grid_computed, report_data);
                runResultsManager(report_data);
                break;

            case 7:
                run_complete_pipeline(&config,
                                      &is_configured,
                                      &single_point_computed,
                                      &grid_computed,
                                      report_data);
                break;

            case 8:
                printf("\nCleaning up resources and exiting application...\n");
                if (report_data != NULL)
                {
                    free(report_data);
                    report_data = NULL;
                }
                printf("Electrostatic Simulation Suite terminated successfully. Goodbye!\n");
                return 0;

            default:
                printf("Invalid choice (%d). Please select an option between 1 and 8.\n", choice);
                break;
        }
    }

    /* Fallback cleanup on EOF */
    if (report_data != NULL)
    {
        free(report_data);
        report_data = NULL;
    }
    return 0;
}
