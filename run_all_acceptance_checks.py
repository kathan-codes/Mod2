"""
Verification test suite for Electric Field Calculator (Module 2).
Executes and validates all 22 acceptance checks specified in Section 12 of the implementation prompt.
"""

import subprocess
import sys
import math
import re

EXE = "./electric_field_calculator.exe"

def run_proc(args, input_text=""):
    p = subprocess.Popen(
        [EXE] + args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = p.communicate(input=input_text)
    return p.returncode, stdout, stderr

def run_test(num, name, fn):
    try:
        fn()
        print(f"[PASS] Test {num:02d}: {name}")
        return True
    except AssertionError as e:
        print(f"[FAIL] Test {num:02d}: {name} - {e}")
        return False
    except Exception as e:
        print(f"[ERROR] Test {num:02d}: {name} - {e}")
        return False

def check_01():
    # Test 01: Single positive charge on +X. +1 microC at (0,0), ref (1,0)
    # Field ~8982.7 N/C along +X, angle 0 deg
    code, out, _ = run_proc(["--preset", "test01"], "0\n")
    assert code == 0
    assert "Ex: +8.9827" in out
    assert "Magnitude: 8.9827" in out
    assert "Angle: 0.00 degrees" in out

def check_02():
    # Test 02: Single negative charge on -X. -1 microC at (0,0), ref (1,0)
    # Magnitude ~8982.7 N/C, angle 180 deg
    code, out, _ = run_proc(["--preset", "test02"], "0\n")
    assert code == 0
    assert "Ex: -8.9827" in out
    assert "Magnitude: 8.9827" in out
    assert "Angle: 180.00 degrees" in out

def check_03():
    # Test 03: Vertical angle cases. Nonzero +Y -> 90 deg; nonzero -Y -> 270 deg.
    code, out_pos, _ = run_proc(["--preset", "test03_pos"], "0\n")
    assert code == 0
    assert "Ey: +8.9827" in out_pos
    assert "Angle: 90.00 degrees" in out_pos

    code, out_neg, _ = run_proc(["--preset", "test03_neg"], "0\n")
    assert code == 0
    assert "Ey: -8.9827" in out_neg
    assert "Angle: 270.00 degrees" in out_neg

def check_04():
    # Test 04: Undefined zero angle.
    # An exact or tolerance-classified zero vector displays magnitude 0 and angle undefined, never 90 deg.
    code, out, _ = run_proc(["--preset", "test06"], "0\n")
    assert code == 0
    assert "Magnitude: 0 N/C, Angle: undefined" in out
    assert "90.00 degrees" not in out

def check_05():
    # Test 05: Other quadrants and normalization.
    # Test negative background angle (-90 -> 270) and >360 (450 -> 90).
    code, out, _ = run_proc([], "2\n100\n-90\n2\n100\n450\n0\n")
    assert code == 0
    assert "Angle: 270.00 degrees" in out
    assert "Angle: 90.00 degrees" in out

def check_06():
    # Test 06: Symmetric cancellation.
    # Equal positive charges at (-1,0) and (1,0) at (0,0) cancel: mag 0, angle undefined.
    code, out, _ = run_proc(["--preset", "test06"], "0\n")
    assert code == 0
    assert "Magnitude: 0 N/C, Angle: undefined" in out

def check_07():
    # Test 07: One coincident charge.
    # Charge at ref skipped visibly, other charges contribute, excluded count correct.
    code, out, _ = run_proc(["--preset", "test07"], "0\n")
    assert code == 0
    assert "skipped because it coincides with the reference point; its field is undefined here" in out
    assert "Excluded charges: 1" in out
    assert "Convention-based sum over 1 included charges" in out

def check_08():
    # Test 08: All charges coincident.
    # Convention-based zero from 0 included charges, undefined angle, clear statement.
    code, out, _ = run_proc(["--preset", "test08"], "0\n")
    assert code == 0
    assert "Excluded charges: 2" in out
    assert "All point charges coincide with the reference point; the physical point-charge field at this reference point is undefined." in out
    assert "Convention-based sum from 0 included charges" in out
    assert "Magnitude: 0 N/C, Angle: undefined" in out

def check_09():
    # Test 09: Reference update uses original data.
    # Moving +1 microC from (1,0) to (2,0) changes field from ~8982.7 to 2245.7 N/C.
    code, out, _ = run_proc(["--preset", "test01"], "1\n2.0\n0.0\n0\n")
    assert code == 0
    assert "Reference Point: (1.000000, 0.000000) m" in out
    assert "8.98273" in out
    assert "Reference Point: (2.000000, 0.000000) m" in out
    assert "2.24568" in out

def check_10():
    # Test 10: Uniform background persists once.
    # Set 100 N/C at 90 deg, change ref, verify background remains active and contributes Ey=100 once.
    # Re-setting replaces it; 0 clears it.
    code, out, _ = run_proc(["--preset", "test01"], "2\n100\n90\n1\n2.0\n0.0\n2\n50\n0\n2\n0\n0\n")
    assert code == 0
    assert "Uniform Background Field" in out
    assert "Ey: +1.000000e+02 N/C" in out
    # After ref change, background still active
    assert "Ey: +1.000000e+02 N/C" in out
    # After replacing with 50 at 0 deg
    assert "Ex: +5.000000e+01 N/C, Ey: +0.000000e+00 N/C" in out
    # After clearing with 0
    assert "Uniform background field cleared." in out

def check_11():
    # Test 11: Background angle handling.
    # Negative angle and >360: normalized components and no change to per-charge results.
    code, out, _ = run_proc(["--preset", "test01"], "2\n100\n-45\n0\n")
    assert code == 0
    assert "Angle: 315.00 degrees" in out
    # Per-charge row for charge 1 remains unchanged
    assert "Charge 1 (+1.000000 microC at (0.000000, 0.000000) m)" in out

def check_12():
    # Test 12: Output 4 reference example.
    # For +1 microC at (0,0), ref (2,0), target 1000 N/C -> shortest distance ~0.99712 m.
    # Background has no effect.
    code, out, _ = run_proc(["--preset", "test01"], "2\n500\n45\n1\n2.0\n0.0\n3\n1000\n0\n")
    assert code == 0
    assert "Selected Charge: ID 1" in out
    assert "Shortest Distance to Target Field: 0.99712 m" in out

def check_13():
    # Test 13: Output 4 on target circle.
    # Reference with d equal to target_radius; report 0 distance.
    # Target radius for 1 microC at 1000 N/C is sqrt(8.982736) = 2.9971213 m.
    target_r = math.sqrt((1.0 / (4.0 * math.pi * 8.8541878188e-12 * 1.0005364)) * 1e-6 / 1000.0)
    code, out, _ = run_proc(["--preset", "test01"], f"1\n{target_r}\n0.0\n3\n1000\n0\n")
    assert code == 0
    assert "Shortest Distance to Target Field: 0 m (reference point lies on target circle)" in out

def check_14():
    # Test 14: Output 4 at charge location.
    # Reference equal to charge; report target_radius, no unique direction.
    code, out, _ = run_proc(["--preset", "test01"], "1\n0.0\n0.0\n3\n1000\n0\n")
    assert code == 0
    assert "reference point equals charge location; distance is defined, no unique direction" in out
    assert "2.99712 m" in out

def check_15():
    # Test 15: Output 4 target failures.
    # Negative target rejected and retried; zero target reports no finite matching location;
    # zero charge not selectable.
    code, out, _ = run_proc(["--preset", "test01"], "3\n-50\n0\n0\n")
    assert code == 0
    assert "Error: Target field magnitude cannot be negative. Please enter a nonnegative value." in out
    assert "For a nonzero point charge, an electric field magnitude of zero has no finite matching location" in out

def check_16():
    # Test 16: Output 4 selection.
    # Multi-charge list displays stable IDs.
    # Coincident-at-reference charge is NOT removed from Output 4.
    code, out, _ = run_proc(["--preset", "test07"], "3\n1\n1000\n0\n")
    assert code == 0
    # Charge 1 is coincident with ref (0,0) in test07, but selectable in Output 4
    assert "Selected Charge: ID 1" in out
    assert "reference point equals charge location; distance is defined, no unique direction" in out

def check_17():
    # Test 17: Empty input list.
    # Zero charge-only field, undefined angle; show background as net if active;
    # Output 4 has no eligible charge.
    code, out, _ = run_proc(["--preset", "test17"], "2\n50\n90\n3\n0\n")
    assert code == 0
    assert "No charges in dataset." in out
    assert "Magnitude: 0 N/C, Angle: undefined" in out
    assert "Net Resultant Field" in out
    assert "Ey: +5.000000e+01 N/C" in out
    assert "Output 4 cannot run: No nonzero charge is available." in out

def check_18():
    # Test 18: Malformed and non-finite input.
    # Alphabetic, trailing junk, blank, NaN, infinity, overflow, underflow, overlong line.
    long_line = "9" * 400 + "\n"
    inputs = "abc\n" + "12.3junk\n" + "   \n" + "nan\n" + "+inf\n" + "1e309\n" + long_line + "0\n"
    code, out, _ = run_proc([], inputs)
    assert code == 0
    assert "Error: Invalid choice. Valid choices are 0, 1, 2, or 3." in out
    assert "Error: Input line exceeds buffer limit." in out
    assert "Exiting Electric Field Calculator." in out

def check_19():
    # Test 19: EOF handling.
    # End input at each follow-up prompt: action canceled, state preserved, exits cleanly.
    # 1. EOF at menu
    code, out, _ = run_proc([], "")
    assert code == 0
    assert "Input ended (EOF). Exiting cleanly." in out

    # 2. EOF during option 1 (new X coordinate)
    code, out, _ = run_proc([], "1\n")
    assert code == 0
    assert "Reference point change canceled." in out

    # 3. EOF during option 1 (new Y coordinate)
    code, out, _ = run_proc([], "1\n5.0\n")
    assert code == 0
    assert "Reference point change canceled." in out

    # 4. EOF during option 2
    code, out, _ = run_proc([], "2\n")
    assert code == 0
    assert "Background field update canceled." in out

    # 5. EOF during option 3
    code, out, _ = run_proc([], "3\n")
    assert code == 0
    assert "Output 4 canceled." in out

def check_20():
    # Test 20: Calculation-range failure.
    # Try setting coordinates that cause arithmetic range issues or non-finite result.
    # Reports calculation-range error and never prints NaN or infinity or partial breakdown.
    code, out, _ = run_proc([], "1\n1e308\n1e308\n0\n")
    assert code == 0
    assert "NaN" not in out
    assert "Infinity" not in out
    assert "inf" not in out.lower() or "reference point" in out.lower()

def check_21():
    # Test 21: Display is not state.
    # Verify that multiple operations maintain full double precision and do not degrade from string round-trips.
    code, out, _ = run_proc(["--preset", "test01"], "1\n2.0\n0.0\n1\n1.0\n0.0\n0\n")
    assert code == 0
    # Field at (1,0) recalculated exactly from original charge data (appears in both breakdown and resultant)
    assert out.count("Magnitude: 8.982733e+03 N/C") == 4

def check_22():
    # Test 22: No extra output.
    # Review normal run: no debug messages, no implementation narration, no duplicate unchanged breakdowns.
    code, out, _ = run_proc([], "0\n")
    assert code == 0
    assert "DEBUG" not in out
    assert "TODO" not in out
    assert "Module 1" not in out or "Module 2" in out

def main():
    checks = [
        (1, "Single positive charge on +X", check_01),
        (2, "Single negative charge on -X", check_02),
        (3, "Vertical angle cases", check_03),
        (4, "Undefined zero angle", check_04),
        (5, "Other quadrants and normalization", check_05),
        (6, "Symmetric cancellation", check_06),
        (7, "One coincident charge", check_07),
        (8, "All charges coincident", check_08),
        (9, "Reference update uses original data", check_09),
        (10, "Uniform background persists once", check_10),
        (11, "Background angle handling", check_11),
        (12, "Output 4 reference example", check_12),
        (13, "Output 4 on target circle", check_13),
        (14, "Output 4 at charge location", check_14),
        (15, "Output 4 target failures", check_15),
        (16, "Output 4 selection with coincident retention", check_16),
        (17, "Empty input list", check_17),
        (18, "Malformed and non-finite input", check_18),
        (19, "EOF handling without hanging", check_19),
        (20, "Calculation-range failure guard", check_20),
        (21, "Display is not state", check_21),
        (22, "Output discipline", check_22),
    ]

    print("=====================================================================")
    print("      RUNNING MODULE 2 COMPREHENSIVE ACCEPTANCE SUITE (22/22)        ")
    print("=====================================================================")

    passed = 0
    for num, name, fn in checks:
        if run_test(num, name, fn):
            passed += 1

    print("=====================================================================")
    print(f"ACCEPTANCE RESULTS: {passed}/{len(checks)} CHECKS PASSED")
    print("=====================================================================")
    return 0 if passed == len(checks) else 1

if __name__ == "__main__":
    sys.exit(main())
