# Tunneling and zero-point energy in enzymatic C-H activation

Code for the paper "Hydrogen Tunneling and Zero-Point Energy in Enzymatic C-H Activation: A WKB Study of the Kinetic Isotope Effect" by Enes Yigit Cebeci.

The code models proton and deuteron tunneling through an Eckart barrier with the WKB approximation. It averages the transmission over a Boltzmann distribution to get rate constants, scans barrier height and width, and adds the zero-point energy difference between C-H and C-D to the kinetic isotope effect (KIE).

## Run it

    pip install -r requirements.txt
    ./run_all.sh

This runs the validation and all figure scripts and saves the figures in `figures/`. It takes about 25 minutes, most of it for Figure 5.

## Files

- `wkb_core.py`: Eckart potential, turning points, WKB transmission, thermal rate constant.
- `validate_wkb.py`: compares the WKB transmission with the exact solution for the symmetric Eckart barrier (32 cases, WKB is 94% to 99% of the exact value).
- `rate_sweep.py`: rate constants with and without ZPE in the barrier, Arrhenius fits, barrier width (FWHM).
- `figure1_electron_proton.py`, `figure2_enzyme.py`, `figure3_KIE.py`: transmission probabilities for different masses and barriers.
- `figure4_arrhenius_final.py`: tunneling-corrected Arrhenius plot for a 1.5 eV barrier with L = 0.25 A.
- `figure5.py`: KIE at 300 K for 117 barriers, tunneling alone and with ZPE in the barrier (takes about 20 minutes).
- `figure6_zpe.py`: Arrhenius plot, KIE against temperature and Table 1 (A_H/A_D, activation energies) with ZPE in the barrier.

## Main results

- The WKB code matches the exact transmission of the symmetric Eckart barrier to within 6% (32 cases).
- Tunneling alone gives a KIE above 7 at 300 K in only 12 of 117 barriers.
- With the zero-point energy (ZPE) of C-H and C-D included in the barrier, 80 of 117 barriers give a KIE between 7 and 100 at 300 K.
- A 0.9 eV barrier with L = 0.40 A (about 1.4 A wide) gives a KIE of 8.2 at 300 K and 21.6 at 220 K, with A_H/A_D = 0.75.
- A very narrow barrier (1.5 eV, L = 0.25 A) gives a KIE of 181 at 300 K, but an activation energy difference of 9.7 kcal/mol and A_H/A_D of about 1e-5.

## Limits

The model is one-dimensional with a symmetric, static barrier and the same prefactor for both isotopes. ZPE is included by lowering each isotope's barrier by the full zero-point energy of its stretch, which assumes the stretch is lost at the transition state. The C-H and C-D frequencies (3000 and 2250 cm^-1) are typical textbook values. The barrier parameters were scanned and not fitted to a specific enzyme.

## License

MIT. See `LICENSE`.
