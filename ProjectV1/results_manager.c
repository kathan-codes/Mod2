#include "results_manager.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define REPORT_INDEX_FILE "saved_reports.txt"

static int valid_report_data(const ReportData *data)
{
    if (data == NULL)
    {
        printf("No simulation data was supplied.\n");
        return 0;
    }

    if (data->number_of_charges < 1 ||
        data->number_of_charges > REPORT_MAX_CHARGES)
    {
        printf("Invalid number of charges in report data.\n");
        return 0;
    }

    if (data->has_grid &&
        (data->grid_width < 1 || data->grid_width > REPORT_MAX_GRID ||
         data->grid_height < 1 || data->grid_height > REPORT_MAX_GRID))
    {
        printf("Invalid grid dimensions in report data.\n");
        return 0;
    }

    return 1;
}

static void print_timestamp(FILE *out)
{
    time_t now = time(NULL);
    struct tm *local = localtime(&now);
    char buffer[64];

    if (local != NULL &&
        strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", local) > 0)
    {
        fprintf(out, "%s", buffer);
    }
    else
    {
        fprintf(out, "Unavailable");
    }
}

void results_init(ReportData *data)
{
    if (data == NULL)
    {
        return;
    }

    memset(data, 0, sizeof(*data));
    strcpy(data->simulation_id, "SIM_UNNAMED");
    strcpy(data->configuration_name, "Unspecified configuration");
}

void displaySimulationResults(const ReportData *data)
{
    int i;

    if (!valid_report_data(data))
    {
        return;
    }

    printf("\n====================================================\n");
    printf("              CONSOLIDATED SIMULATION RESULTS\n");
    printf("====================================================\n");
    printf("Simulation ID      : %s\n", data->simulation_id);
    printf("Configuration      : %s\n", data->configuration_name);
    printf("Number of charges  : %d\n", data->number_of_charges);

    if (data->number_of_charges > 0)
    {
        printf("\nCHARGE CONFIGURATION\n");
        printf("----------------------------------------------------\n");
        for (i = 0; i < data->number_of_charges; ++i)
        {
            printf("Charge %d: q = %+.6e C, position = (%+.4f, %+.4f) m\n",
                   i + 1, data->charges[i].q,
                   data->charges[i].x, data->charges[i].y);
        }
    }

    if (data->has_single_point)
    {
        printf("\nSINGLE-POINT ELECTRIC FIELD RESULTS\n");
        printf("----------------------------------------------------\n");
        printf("Reference point : (%+.4f, %+.4f) m\n",
               data->reference_x, data->reference_y);
        printf("Charge-only Ex  : %+.6e N/C\n", data->charge_Ex);
        printf("Charge-only Ey  : %+.6e N/C\n", data->charge_Ey);

        if (data->charge_field_defined)
        {
            printf("Charge-only |E| : %.6e N/C\n",
                   data->charge_field_magnitude);
            printf("Charge-only angle: %.3f degrees\n",
                   data->charge_field_angle);
        }
        else
        {
            printf("Charge-only direction: undefined/zero field\n");
        }

        if (data->external_field_active)
        {
            printf("External Ex     : %+.6e N/C\n", data->external_Ex);
            printf("External Ey     : %+.6e N/C\n", data->external_Ey);
        }

        printf("Net Ex          : %+.6e N/C\n", data->net_Ex);
        printf("Net Ey          : %+.6e N/C\n", data->net_Ey);

        if (data->net_field_defined)
        {
            printf("Net |E|         : %.6e N/C\n",
                   data->net_field_magnitude);
            printf("Net angle       : %.3f degrees\n",
                   data->net_field_angle);
        }
        else
        {
            printf("Net field direction: undefined/zero field\n");
        }
    }
    else
    {
        printf("\nSingle-point electric field results: not supplied.\n");
    }

    if (data->has_grid)
    {
        printf("\nGRID DATA\n");
        printf("----------------------------------------------------\n");
        printf("Grid dimensions  : %d x %d\n",
               data->grid_width, data->grid_height);
        printf("Grid locations   : %d\n",
               data->grid_width * data->grid_height);
    }
    else
    {
        printf("\nGrid data: not supplied.\n");
    }

    printf("====================================================\n");
}

static void write_report_body(FILE *fp, const ReportData *data)
{
    int i;

    fprintf(fp, "====================================================\n");
    fprintf(fp, "             ELECTROSTATIC SIMULATION REPORT\n");
    fprintf(fp, "====================================================\n\n");

    fprintf(fp, "SIMULATION INFORMATION\n");
    fprintf(fp, "----------------------------------------------------\n");
    fprintf(fp, "Simulation ID : %s\n", data->simulation_id);
    fprintf(fp, "Configuration : %s\n", data->configuration_name);
    fprintf(fp, "Generated     : ");
    print_timestamp(fp);
    fprintf(fp, "\n\n");

    fprintf(fp, "CHARGE CONFIGURATION\n");
    fprintf(fp, "----------------------------------------------------\n");
    fprintf(fp, "Number of charges: %d\n", data->number_of_charges);
    for (i = 0; i < data->number_of_charges; ++i)
    {
        fprintf(fp, "Charge %d: q = %+.9e C, x = %+.6f m, y = %+.6f m\n",
                i + 1, data->charges[i].q,
                data->charges[i].x, data->charges[i].y);
    }

    fprintf(fp, "\nSINGLE-POINT ELECTRIC FIELD\n");
    fprintf(fp, "----------------------------------------------------\n");
    if (data->has_single_point)
    {
        fprintf(fp, "Reference point: (%+.6f, %+.6f) m\n",
                data->reference_x, data->reference_y);
        fprintf(fp, "Charge-only Ex: %+.9e N/C\n", data->charge_Ex);
        fprintf(fp, "Charge-only Ey: %+.9e N/C\n", data->charge_Ey);
        if (data->charge_field_defined)
        {
            fprintf(fp, "Charge-only magnitude: %.9e N/C\n",
                    data->charge_field_magnitude);
            fprintf(fp, "Charge-only angle: %.6f degrees\n",
                    data->charge_field_angle);
        }
        else
        {
            fprintf(fp, "Charge-only field: Undefined\n");
        }

        if (data->external_field_active)
        {
            fprintf(fp, "External Ex: %+.9e N/C\n", data->external_Ex);
            fprintf(fp, "External Ey: %+.9e N/C\n", data->external_Ey);
        }

        fprintf(fp, "Net Ex: %+.9e N/C\n", data->net_Ex);
        fprintf(fp, "Net Ey: %+.9e N/C\n", data->net_Ey);
        fprintf(fp, "Net magnitude: %.9e N/C\n",
                data->net_field_magnitude);

        if (data->net_field_defined)
        {
            fprintf(fp, "Net angle: %.6f degrees\n",
                    data->net_field_angle);
        }
        else
        {
            fprintf(fp, "Net angle: undefined\n");
        }
    }
    else
    {
        fprintf(fp, "Not supplied by the electric-field module.\n");
    }

    fprintf(fp, "\nGRID SUMMARY\n");
    fprintf(fp, "----------------------------------------------------\n");
    if (data->has_grid)
    {
        fprintf(fp, "Grid dimensions: %d x %d\n",
                data->grid_width, data->grid_height);
        fprintf(fp, "Total grid locations: %d\n",
                data->grid_width * data->grid_height);
        fprintf(fp, "Invalid field points are marked in the CSV as field_valid=0.\n");
    }
    else
    {
        fprintf(fp, "Grid data not supplied.\n");
    }

    fprintf(fp, "\n====================================================\n");
    fprintf(fp, "                    END OF REPORT\n");
    fprintf(fp, "====================================================\n");
}

int exportTXT(const ReportData *data, const char *filename)
{
    FILE *fp;

    if (!valid_report_data(data) || filename == NULL || filename[0] == '\0')
    {
        printf("Cannot export TXT: invalid data or filename.\n");
        return 0;
    }

    fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Could not create %s.\n", filename);
        return 0;
    }

    write_report_body(fp, data);

    if (fclose(fp) != 0)
    {
        printf("An error occurred while closing %s.\n", filename);
        return 0;
    }

    if (!registerSavedReport(filename))
    {
        printf("Warning: Report saved, but could not add it to the saved reports list.\n");
    }

    printf("Text report saved to %s\n", filename);
    return 1;
    printf("Text report saved to %s\n", filename);
    return 1;
}

int exportCSV(const ReportData *data, const char *filename)
{
    FILE *fp;
    int i, j;

    if (!valid_report_data(data) || filename == NULL || filename[0] == '\0')
    {
        printf("Cannot export CSV: invalid data or filename.\n");
        return 0;
    }

    if (!data->has_grid)
    {
        printf("No grid data is available to export.\n");
        return 0;
    }

    fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Could not create %s.\n", filename);
        return 0;
    }

    fprintf(fp, "x_m,y_m,Ex_N_per_C,Ey_N_per_C,E_magnitude_N_per_C,"
                "angle_degrees,potential_V,field_valid,potential_valid\n");

    for (j = 0; j < data->grid_height; ++j)
    {
        for (i = 0; i < data->grid_width; ++i)
        {
            const ReportGridPoint *p = &data->grid[j][i];

            fprintf(fp, "%.9g,%.9g", p->x, p->y);

            if (p->field_valid)
            {
                fprintf(fp, ",%.9g,%.9g,%.9g,%.9g",
                        p->Ex, p->Ey, p->magnitude, p->angle);
            }
            else
            {
                fprintf(fp, ",,,,");
            }

            if (p->potential_valid)
            {
                fprintf(fp, ",%.9g", p->potential);
            }
            else
            {
                fprintf(fp, ",");
            }

            fprintf(fp, ",%d,%d\n",
                    p->field_valid ? 1 : 0,
                    p->potential_valid ? 1 : 0);
        }
    }

    if (fclose(fp) != 0)
    {
        printf("An error occurred while closing %s.\n", filename);
        return 0;
    }

    if (!registerSavedReport(filename))
    {
        printf("Warning: CSV saved, but could not add it to the saved reports list.\n");
    }

    printf("Grid data exported to %s\n", filename);
    return 1;
    printf("Grid data exported to %s\n", filename);
    return 1;
}

void displayStatistics(const ReportData *data)
{
    double minimum = 0.0, maximum = 0.0, sum = 0.0;
    size_t valid_count = 0;
    int i, j;

    if (!valid_report_data(data))
    {
        return;
    }

    if (!data->has_grid)
    {
        printf("No grid data is available for statistical analysis.\n");
        return;
    }

    for (j = 0; j < data->grid_height; ++j)
    {
        for (i = 0; i < data->grid_width; ++i)
        {
            const ReportGridPoint *p = &data->grid[j][i];

            if (!p->field_valid || !isfinite(p->magnitude))
            {
                continue;
            }

            if (valid_count == 0)
            {
                minimum = maximum = p->magnitude;
            }
            else
            {
                if (p->magnitude < minimum)
                    minimum = p->magnitude;
                if (p->magnitude > maximum)
                    maximum = p->magnitude;
            }

            sum += p->magnitude;
            ++valid_count;
        }
    }

    printf("\n========== GRID FIELD STATISTICS ==========\n");
    printf("Valid grid points : %zu\n", valid_count);
    printf("Invalid grid points: %d\n",
           data->grid_width * data->grid_height - (int)valid_count);

    if (valid_count == 0)
    {
        printf("No valid field values are available.\n");
    }
    else
    {
        printf("Minimum |E|      : %.6e N/C\n", minimum);
        printf("Maximum |E|      : %.6e N/C\n", maximum);
        printf("Average |E|      : %.6e N/C\n",
               sum / (double)valid_count);
    }
    printf("===========================================\n");
}

int registerSavedReport(const char *filename)
{
    FILE *fp;

    if (filename == NULL || filename[0] == '\0')
    {
        return 0;
    }

    fp = fopen(REPORT_INDEX_FILE, "a");
    if (fp == NULL)
    {
        return 0;
    }
    if (strchr(filename, '\n') != NULL || strchr(filename, '\r') != NULL)
    {
        fclose(fp);
        return 0;
    }

    fprintf(fp, "%s\n", filename);
    return fclose(fp) == 0;
}

void viewSavedReports(void)
{
    FILE *fp = fopen(REPORT_INDEX_FILE, "r");
    char line[256];
    int found = 0;

    printf("\n========== SAVED REPORTS ==========\n");

    if (fp == NULL)
    {
        printf("No report index exists yet. Export a report first.\n");
        printf("===================================\n");
        return;
    }

    while (fgets(line, sizeof(line), fp) != NULL)
    {
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n')
        {
            line[len - 1] = '\0';
        }
        if (line[0] != '\0')
        {
            printf("%d. %s\n", ++found, line);
        }
    }

    fclose(fp);

    if (!found)
    {
        printf("No saved reports are listed.\n");
    }
    printf("===================================\n");
}

void runResultsManager(const ReportData *data)
{
    int choice;
    char filename[160];

    if (!valid_report_data(data))
    {
        return;
    }

    for (;;)
    {
        printf("\n========================================\n");
        printf("       RESULTS AND REPORT MANAGER\n");
        printf("========================================\n");
        printf("1. View consolidated simulation results\n");
        printf("2. View grid statistics\n");
        printf("3. Export TXT report\n");
        printf("4. Export grid data as CSV\n");
        printf("5. View saved report filenames\n");
        printf("0. Return to main menu\n");
        printf("----------------------------------------\n");
        printf("Enter choice: ");

        if (scanf("%d", &choice) != 1)
        {
            int ch;
            while ((ch = getchar()) != '\n' && ch != EOF)
            {
                /* discard invalid input */
            }
            printf("Invalid input. Enter a menu number.\n");
            continue;
        }

        switch (choice)
        {
        case 1:
            displaySimulationResults(data);
            break;

        case 2:
            displayStatistics(data);
            break;

        case 3:
            printf("Enter TXT filename (e.g. simulation_001.txt): ");
            if (scanf("%159s", filename) == 1)
            {
                exportTXT(data, filename);
            }
            else
            {
                printf("Invalid filename.\n");
            }
            break;

        case 4:
            printf("Enter CSV filename (e.g. simulation_001.csv): ");
            if (scanf("%159s", filename) == 1)
            {
                exportCSV(data, filename);
            }
            else
            {
                printf("Invalid filename.\n");
            }
            break;

        case 5:
            viewSavedReports();
            break;

        case 0:
            return;

        default:
            printf("Invalid choice. Select 0 to 5.\n");
        }
    }
}
