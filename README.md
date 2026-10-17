# Gram-Schmidt, constellation diagrams and BER analysis (MATLAB)

MATLAB course project for ECE352s Analog and Digital Communication Systems, Ain Shams University, Spring 2026 (Team 11).

- **Part 1: Gram-Schmidt and constellations.** A small GUI takes M arbitrary signals, typed as MATLAB expressions in *t*, and finds the minimum set of N orthonormal basis functions with the Gram-Schmidt procedure. It plots the signals and basis functions, draws the 2-D or 3-D constellation diagram, and checks each signal's energy in the time domain against its energy from the constellation.
- **Part 2: BER over AWGN.** It simulates Unipolar NRZ and Manchester line coding over an AWGN channel with a correlation receiver, sweeps Eb/N0, and compares the simulated BER with theory: Q(√(Eb/N0)) for Unipolar NRZ and Q(√(2Eb/N0)) for Manchester, which is antipodal.

**Result:** the simulated BER points follow the theoretical curves. Manchester beats Unipolar NRZ by 3 dB because antipodal symbols are twice as far apart for the same bit energy. The plots are in the [project report](docs/project-report.pdf), figures 1–11.

## My role (Abdelrahman Sherif)

The whole project except the GUI: the Gram-Schmidt procedure, the constellation diagrams and energy checks (Part 1), and the line coding, AWGN channel, detection and theoretical/simulated BER analysis (Part 2). Teammates built and integrated the GUI and dialogs.

## Files

| File | Contents |
|---|---|
| [`part1_gram_schmidt_constellation.m`](part1_gram_schmidt_constellation.m) | Part 1: GUI, Gram-Schmidt, constellation plot, energy check. Line 2 (`GramSchmidt_GUI();`) was added after submission; without it the script never opened the GUI |
| [`part2_ber_line_coding.m`](part2_ber_line_coding.m) | Part 2: parameter dialogs, line coding, AWGN channel, detection, BER curves |
| [`docs/project-report.pdf`](docs/project-report.pdf) | Team report: method, test signal sets, theoretical derivations, results |
| [`development/`](development) | Earlier working versions of the scripts (`part1_optimized*`, `part2_optimized*`, `Part_2_Final`, `GramSchmidt.m`) |

## How to run

MATLAB R2020b or newer (uses `uifigure`) with the Communications Toolbox (for `qfunc`).

- **Part 1:** run `part1_gram_schmidt_constellation.m`, enter the time span, M and one signal expression per line (the defaults are a 4-signal example), then press the button.
- **Part 2:** run `part2_ber_line_coding.m`, then enter the number of bits and the Eb/N0 range and pick Unipolar NRZ, Manchester or both.

## Team 11

Yossef Medhat Ali, Khaled Abdelrahman Hussien, Abdelrahman Mohamed Aly, Mostafa Mohamed Abdullah, Ganna Abdelhamid Mohamed, Amr Hatem Sabry, Abdelrahman Sherif Medhat, Ali Amir Abouzaid, Mohamed Fouad Mahmoud.

Supervised by Prof. Bassant Abdelhamid, with TA Amr Hammouda.

## License

No license. All rights reserved by Team 11. The code and report are shared for viewing only.
