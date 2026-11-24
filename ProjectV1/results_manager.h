#ifndef RESULTS_MANAGER_H
#define RESULTS_MANAGER_H

#include <stddef.h>

#define REPORT_MAX_CHARGES 100
#define REPORT_MAX_GRID 100

typedef struct
{
    double q;
    double x;
    double y;
} ReportCharge;

typedef struct
{
    double x;
    double y;
    double Ex;
    double Ey;
    double magnitude;
    double angle;
    double potential;
    int field_valid;
    int potential_valid;
} ReportGridPoint;

typedef struct
{
    char simulation_id[64];
    char configuration_name[100];

    int number_of_charges;
    ReportCharge charges[REPORT_MAX_CHARGES];

    int has_single_point;
    double reference_x;
    double reference_y;
    double charge_Ex;
    double charge_Ey;
    double charge_field_magnitude;
    double charge_field_angle;
    int charge_field_defined;

    int external_field_active;
    double external_Ex;
    double external_Ey;

    double net_Ex;
    double net_Ey;
    double net_field_magnitude;
    double net_field_angle;
    int net_field_defined;

    int grid_width;
    int grid_height;
    ReportGridPoint grid[REPORT_MAX_GRID][REPORT_MAX_GRID];
    int has_grid;
} ReportData;

void results_init(ReportData *data);

void displaySimulationResults(const ReportData *data);

int exportTXT(const ReportData *data, const char *filename);

int exportCSV(const ReportData *data, const char *filename);

void displayStatistics(const ReportData *data);

void viewSavedReports(void);

int registerSavedReport(const char *filename);

void runResultsManager(const ReportData *data);

#endif
