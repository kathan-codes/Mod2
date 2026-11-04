import subprocess
import sys
import math

EXE = "./electric_field_calculator.exe"

def run_proc(args, input_text):
    p = subprocess.Popen(
        [EXE] + args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = p.communicate(input=input_text)
    return p.returncode, stdout, stderr

def test_startup_and_exit():
    code, out, _ = run_proc([], "0\n")
    assert code == 0
    assert "[DEMO MODE]" in out
    assert "Charge 1: +1.000000 microC at (0.000000, 0.000000) m" in out
    assert "Charge 2: -1.000000 microC at (2.000000, 0.000000) m" in out
    assert "Exiting Electric Field Calculator." in out
    print("PASS: test_startup_and_exit")

def test_change_ref_point():
    # Change ref to (1.0, 0.0)
    input_seq = "1\n1.0\n0.0\n0\n"
    code, out, _ = run_proc([], input_seq)
    assert code == 0
    assert "Reference Point: (1.000000, 0.000000) m" in out
    assert "Ex: +1.796547e+04 N/C" in out or "Ex: +1.796547e+004 N/C" in out
    assert "Angle: 0.00 degrees" in out
    print("PASS: test_change_ref_point")

def test_background_field_lifecycle():
    # Set background 100 N/C at 90 deg, then change ref, then clear background (0 mag)
    input_seq = "2\n100\n90\n1\n1.0\n0.0\n2\n0\n0\n"
    code, out, _ = run_proc([], input_seq)
    assert code == 0
    assert "Uniform Background Field" in out
    assert "Ey: +1.000000e+02 N/C" in out
    assert "Net Resultant Field" in out
    assert "Uniform background field cleared." in out
    print("PASS: test_background_field_lifecycle")

def test_output_4_multi_charge():
    # Option 3, select charge 1, target 1000 N/C
    input_seq = "3\n1\n1000\n0\n"
    code, out, _ = run_proc([], input_seq)
    assert code == 0
    assert "Output 4 Result" in out
    assert "Selected Charge: ID 1" in out
    assert "Target Field Magnitude: 1.000000e+03 N/C" in out
    assert "Shortest Distance to Target Field: 1.58291 m" in out
    print("PASS: test_output_4_multi_charge")

def test_output_4_failures():
    # Negative target re-prompted, then 0 returns to menu
    input_seq = "3\n1\n-50\n0\n0\n"
    code, out, _ = run_proc([], input_seq)
    assert code == 0
    assert "Error: Target field magnitude cannot be negative." in out
    assert "For a nonzero point charge, an electric field magnitude of zero has no finite matching location" in out
    print("PASS: test_output_4_failures")

def test_invalid_input_and_overlong():
    long_line = "x" * 300 + "\n"
    input_seq = "abc\n" + "99\n" + "nan\n" + "\n" + long_line + "0\n"
    code, out, _ = run_proc([], input_seq)
    assert code == 0
    assert "Error: Invalid choice. Valid choices are 0, 1, 2, or 3." in out
    assert "Error: Input line exceeds buffer limit." in out
    print("PASS: test_invalid_input_and_overlong")

def test_eof_handling():
    # EOF at menu
    code, out, _ = run_proc([], "")
    assert code == 0
    assert "Input ended (EOF). Exiting cleanly." in out

    # EOF during option 1
    code, out, _ = run_proc([], "1\n1.0\n")
    assert code == 0
    assert "Reference point change canceled." in out

    # EOF during option 2
    code, out, _ = run_proc([], "2\n")
    assert code == 0
    assert "Background field update canceled." in out

    # EOF during option 3
    code, out, _ = run_proc([], "3\n")
    assert code == 0
    assert "Output 4 canceled." in out
    print("PASS: test_eof_handling")

def test_preset_vertical():
    code, out, _ = run_proc(["--preset", "test03_pos"], "0\n")
    assert code == 0
    assert "Angle: 90.00 degrees" in out

    code, out, _ = run_proc(["--preset", "test03_neg"], "0\n")
    assert code == 0
    assert "Angle: 270.00 degrees" in out
    print("PASS: test_preset_vertical")

def test_preset_symmetric_cancellation():
    code, out, _ = run_proc(["--preset", "test06"], "0\n")
    assert code == 0
    assert "Magnitude: 0 N/C, Angle: undefined" in out
    print("PASS: test_preset_symmetric_cancellation")

def test_preset_coincident_and_output4():
    # Presets test07: charge 1 is coincident at ref (0,0)
    # Output 4 should STILL allow selecting charge 1!
    # Target 1000 N/C -> d=0 -> shortest distance = target radius (2.99712 m)
    input_seq = "3\n1\n1000\n0\n"
    code, out, _ = run_proc(["--preset", "test07"], input_seq)
    assert code == 0
    assert "Charge 1 (+1.000000 microC at (0.000000, 0.000000) m): skipped because it coincides with the reference point; its field is undefined here." in out
    assert "Excluded charges: 1" in out
    assert "Shortest Distance to Target Field: 2.99712 m (reference point equals charge location; distance is defined, no unique direction)" in out
    print("PASS: test_preset_coincident_and_output4")

def test_preset_all_coincident():
    code, out, _ = run_proc(["--preset", "test08"], "0\n")
    assert code == 0
    assert "All point charges coincide with the reference point; the physical point-charge field at this reference point is undefined." in out
    assert "Magnitude: 0 N/C, Angle: undefined" in out
    print("PASS: test_preset_all_coincident")

def test_preset_empty():
    code, out, _ = run_proc(["--preset", "test17"], "3\n0\n")
    assert code == 0
    assert "No charges in dataset." in out
    assert "Magnitude: 0 N/C, Angle: undefined" in out
    assert "Output 4 cannot run: No nonzero charge is available." in out
    print("PASS: test_preset_empty")

if __name__ == "__main__":
    test_startup_and_exit()
    test_change_ref_point()
    test_background_field_lifecycle()
    test_output_4_multi_charge()
    test_output_4_failures()
    test_invalid_input_and_overlong()
    test_eof_handling()
    test_preset_vertical()
    test_preset_symmetric_cancellation()
    test_preset_coincident_and_output4()
    test_preset_all_coincident()
    test_preset_empty()
    print("ALL INTERACTIVE END-TO-END TESTS PASSED!")

