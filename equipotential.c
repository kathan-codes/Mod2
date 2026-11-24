#include "equipotential.h"
#include <math.h>
#include <stdio.h>

#ifndef MAX_CHARGES
#define MAX_CHARGES 100
#endif

#define WORLD_W 12.0f
#define WORLD_H 8.4f

float charges[MAX_CHARGES];
int chargeCount = 0;

Vector2 ChargePosition(int index, int columns, int rows)
{
    int row = index / columns;
    int column = index % columns;

    int itemsInRow = chargeCount - row * columns;
    if (itemsInRow > columns) itemsInRow = columns;

    float x = -WORLD_W/2.0f
            + WORLD_W * (column + 1) / (itemsInRow + 1);

    float y = WORLD_H/2.0f
            - WORLD_H * (row + 1) / (rows + 1);

    return (Vector2){x, y};
}

Vector2 WorldToScreen(Vector2 p)
{
    return (Vector2){
        (p.x + WORLD_W/2.0f) * SCREEN_W / WORLD_W,
        (WORLD_H/2.0f - p.y) * SCREEN_H / WORLD_H
    };
}

void RunEquipotentialSimulation(void)
{
    if (chargeCount <= 0) {
        printf("Number of charges (1-%d): ", MAX_CHARGES);

        if (scanf("%d", &chargeCount) != 1 ||
            chargeCount < 1 || chargeCount > MAX_CHARGES) {
            printf("Enter a whole number from 1 to %d.\n", MAX_CHARGES);
            chargeCount = 0;
            return;
        }

        for (int i = 0; i < chargeCount; i++) {
            printf("Enter value for charge %d (positive or negative): ", i + 1);

            if (scanf("%f", &charges[i]) != 1) {
                printf("Invalid charge value.\n");
                chargeCount = 0;
                return;
            }
        }
    }

    float maxCharge = 0.0f;
    for (int i = 0; i < chargeCount; i++) {
        if (fabsf(charges[i]) > maxCharge) {
            maxCharge = fabsf(charges[i]);
        }
    }

    int columns = (int)ceilf(sqrtf((float)chargeCount));
    int rows = (chargeCount + columns - 1) / columns;

    float cellW = WORLD_W / columns;
    float cellH = WORLD_H / rows;
    float largestRadius = 0.38f * fminf(cellW, cellH);

    const float ringFractions[] = {
        0.25f, 0.40f, 0.55f, 0.70f, 0.85f
    };
    const int ringCount =
        sizeof(ringFractions) / sizeof(ringFractions[0]);

#if !defined(NO_RAYLIB)
    InitWindow(SCREEN_W, SCREEN_H, "Equipotential Circles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground((Color){18, 20, 28, 255});

        DrawText("Equipotential circles for each charge individually",
                 18, 18, 20, RAYWHITE);
        DrawText("Charge centers are arranged automatically",
                 18, 44, 16, LIGHTGRAY);

        for (int i = 0; i < chargeCount; i++) {
            Vector2 worldPos = ChargePosition(i, columns, rows);
            Vector2 screenPos = WorldToScreen(worldPos);

            Color color;
            if (charges[i] > 0.0f) color = ORANGE;
            else if (charges[i] < 0.0f) color = SKYBLUE;
            else color = GRAY;

            if (maxCharge > 0.0f && charges[i] != 0.0f) {
                float chargeRatio = fabsf(charges[i]) / maxCharge;
                float pixelsPerUnit = SCREEN_W / WORLD_W;

                for (int ring = 0; ring < ringCount; ring++) {
                    float radiusWorld =
                        largestRadius * ringFractions[ring] * chargeRatio;
                    float radiusPixels = radiusWorld * pixelsPerUnit;

                    DrawCircleLines((int)screenPos.x, (int)screenPos.y,
                                    radiusPixels, color);
                }
            }

            DrawCircleV(screenPos, 7.0f, color);

            DrawText(TextFormat("q%d = %+.3g", i + 1, charges[i]),
                     (int)screenPos.x - 38,
                     (int)(screenPos.y - largestRadius *
                           (SCREEN_W / WORLD_W) - 24),
                     16, RAYWHITE);
        }

        DrawText("Orange: positive    Blue: negative    Gray: zero",
                 18, SCREEN_H - 30, 16, RAYWHITE);

        EndDrawing();
    }

    CloseWindow();
#else
    printf("\n--- Equipotential Circles Summary (Module 5) ---\n");
    printf("Simulated grid: %d columns x %d rows\n", columns, rows);
    for (int i = 0; i < chargeCount; i++) {
        Vector2 worldPos = ChargePosition(i, columns, rows);
        printf("Charge %d: q = %+.3e C at (%.2f, %.2f) m\n", i + 1, charges[i], worldPos.x, worldPos.y);
        if (maxCharge > 0.0f && charges[i] != 0.0f) {
            float chargeRatio = fabsf(charges[i]) / maxCharge;
            printf("  Equipotential ring radii:\n");
            for (int ring = 0; ring < ringCount; ring++) {
                float radiusWorld = largestRadius * ringFractions[ring] * chargeRatio;
                printf("    Ring %d (%.0f%%): %.4f m\n", ring + 1, ringFractions[ring] * 100.0f, radiusWorld);
            }
        }
    }
    printf("------------------------------------------------\n");
#endif
}

#if defined(STANDALONE_BUILD) && !defined(INTEGRATED_BUILD)
int main(void)
{
    RunEquipotentialSimulation();
    return 0;
}
#endif
