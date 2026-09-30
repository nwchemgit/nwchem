# DFTD4 interface for NWChem

## Scope

This interface connects the popular DFTD4 (Version 4.2.0) dispersion corrections to molecular NWChem DFT energies, analytic gradients, and numerical Hessian corrections (and, concurrently, enables DFTD4 corrected geometry optimizations and harmonic vibrational frequency and normal-mode calculations). It calls the DFTD4 C API directly and does not start the `dftd4` command-line program.

The interface recognizes supported named NWChem functionals and native or LibXC exchange/correlation pairs. Native pairs with omitted coefficients use NWChem's default coefficient of 1.0; explicit pair coefficients must also both be 1.0. DFTD4 performs the final parameter lookup, so the interface follows the parameter set in its pinned release instead of maintaining a separate short allowlist. Bq centers are excluded from the DFTD4 structure and their derivative rows and columns remain zero. Existing `DISP` dispersion corrections cannot be combined with DFTD4, and this will trigger an error.

The interface cannot make an electronic functional available when NWChem or LibXC cannot evaluate it, and will not work if the requested functional is not found within the DFTD4 library.

## Build

Use C and Fortran compilers supported by NWChem and ensure that CMake, `curl`, `gzip`, and `tar` are available. Enable the interface in the same environment used for the NWChem build:

```bash
export USE_DFTD4=1
```

Like `USE_LIBXC=1`, `USE_DFTD4=1` downloads a pinned source release and builds a private copy within `src/libext`. The default version is 4.2.0 and can be changed at build time with `DFTD4_VERSION`.

Keep the existing MPI, Simint, LibXC, BLAS, and ScaLAPACK variables. `USE_DFTD4` and `USE_TBLITE` are mutually exclusive because both interfaces provide incompatible DFTD4 libraries. A build uses:

```bash
cd "$NWCHEM_TOP/src"
make -C libext/dftd4
make
```

No external DFTD4 installation path is required. The build downloads the tagged archive, lets DFTD4 fetch its pinned dependencies, and installs static libraries and headers under `src/libext/dftd4/install`.

### Library linkage

The build creates the static bridge archive `lib/LINUX64/libnwc_dftd4.a` and links it with the DFTD4 4.2.0 libraries installed under `src/libext/dftd4/install`. NWChem's tblite/XTB interface links a bundled DFTD4 3.3.0 library with conflicting symbols and a potentially incompatible ABI. The build therefore rejects configurations that enable both `USE_DFTD4` and `USE_TBLITE`.

DFTD4 uses NWChem's selected BLAS/LAPACK libraries and follows `BLAS_SIZE`: `BLAS_SIZE=8` enables DFTD4's ILP64 interface, while `BLAS_SIZE=4` uses LP64. The static bridge requires no DFTD4-specific shared library or runtime RPATH.

## Input syntax

Enable or disable the correction inside the DFT block:

```text
dft
  xc b3lyp
  dftd4 on
end
```

```text
dft
  dftd4 off
end
```

The default is off. Both forms require a build made with `USE_DFTD4=1`.

The correction uses rational damping with the Axilrod–Teller–Muto three-body
term enabled (`s9 = 1`) by default. There is currently no input keyword for
changing the ATM setting or its damping parameters.

### Functional selection

Named NWChem combinations are translated to DFTD4's canonical identifiers. The current mappings include B3LYP, B97, PBE0, TPSSh, MPW1B95, MPWB1K, B1B95, PW6B95, M06-L, M06, B97-D, BHLYP, B3P86, B3PW91, PBE, HSE03, SCAN, rSCAN, r2SCAN, omegaB97, omegaB97X-2008, MN12-SX, and r2SCAN0. A single native XC keyword is otherwise passed to DFTD4 for its own parameter lookup.

LibXC selections are obtained from NWChem's LibXC adapter rather than from input-text guessing. Both composite functionals and supported exchange/correlation pairs are accepted. Pair order is normalized for DFTD4-supported combinations such as BLYP, BP86, PBE, PBEsol, TPSS, revTPSS, SCAN, rSCAN, r2SCAN, M06-L, and related GGAs. Each LibXC term must have unit weight; scaled or arbitrary mixtures are rejected because they do not identify a published DFTD4 parameterization.

Supported native NWChem exchange/correlation pairs are normalized in the same way. For example, `xc xpbe96 cpbe96` and `xc xpbe96 1.0 cpbe96 1.0` select the same DFTD4 parameters as `xc pbe96`. Scaled native pairs are rejected.

Examples are:

```text
dft
  xc hyb_gga_xc_cam_b3lyp
  dftd4 on
end

dft
  xc gga_c_pbe gga_x_pbe
  dftd4 on
end
```

```text
start water_d4
geometry units angstrom noautoz
 O  0.000000  0.000000  0.000000
 H  0.000000  0.757160  0.586260
 H  0.000000 -0.757160  0.586260
end
basis spherical
 * library def2-svp
end
dft
 xc b3lyp
 dftd4 on
end
task dft energy
```

### Energy

The output uses a `DFTD4 Dispersion Correction` block. It identifies the mapped functional, reports the loaded rational-damping parameters (`s6`, `s8`, `s9`, `a1`, `a2`, and `alp`), and prints the energy to 12 decimal places.

Coordinates are converted to the internal NWChem Bohr representation before
the DFTD4 call. Thus the interface gives the same physical result whether the
input geometry is declared in angstroms or Bohr.

### Optimization

```text
dft
 xc pbe0
 dftd4 on
 convergence energy 1d-8 density 1d-7 gradient 1d-6
end
driver
 tight
end
task dft optimize
```

The joint nuclear gradient is the direct sum
`g_total = g_KS-DFT + g_DFTD4`. The DFTD4 analytic gradient is evaluated once
on the master Global Arrays rank, reduced through NWChem's existing GA path,
and added once to the ordinary Kohn–Sham DFT gradient. It is a post-SCF energy
correction and does not modify the Kohn–Sham Fock matrix or orbitals.
The `DFTD4 DISPERSION GRADIENT` block prints fixed-point components to 12 decimal places, matching the precision of the DFTD4 energy block. Bq centers have zero DFTD4 gradient rows.

### Hessian and frequencies

```text
dft
 xc r2scan
 dftd4 on
end
task dft frequencies
```

For an analytic molecular DFT Hessian, rank zero obtains the numerical DFTD4
Hessian from the 4.2.0 API and adds it to NWChem's Hessian file before
vibrational analysis. DFTD4 4.2.0 exposes no native analytic Hessian API; its
public Hessian routine is finite difference. All ranks participate in NWChem
geometry lifecycle and synchronization calls. If NWChem selects a numerical
Hessian, the analytic DFTD4 gradient is naturally included in each displaced
gradient.

## Parallel behavior

DFTD4 evaluation occurs only on Global Arrays rank zero; energy and gradient data are distributed through NWChem's existing GA conventions. Hessian file modification occurs only on rank zero, while all ranks participate in collective geometry setup, cleanup, and synchronization. Run NWChem with the process layout appropriate for the selected MPI and Global Arrays configuration.

## Validation

- S66 energies: cases 01, 23, 47, and 66 for r2SCAN, B3LYP, PBE0, and PW6B95 matched independent standalone DFTD4 4.2.0 references to at most `1.1e-16` Hartree when identical Bohr coordinates were supplied.
- Gradient: analytic DFTD4 gradient agreed with central differences of DFTD4 energy within `2.16e-13` Hartree/Bohr.
- Hessian: the propagated correction agreed with central differences of analytic DFTD4 gradients within `4.48e-11` Hartree/Bohr².
- Parallel consistency: gradients were identical across the tested process counts; Hessian corrections differed by at most `1.0e-11` Hartree/Bohr², the precision of NWChem's formatted Hessian file.
- Frequencies: complete B3LYP frequency calculations at multiple process counts included the DFTD4 Hessian contribution.
- Bq: adding a ghost center left the DFTD4 energy bit-for-bit unchanged.
- Conflict handling: `dftd4 on` with `disp vdw 4` failed with the intended double-dispersion diagnostic.

## Source layout

- `src/libext/dftd4/nw_dftd4.c`: narrow C API bridge and integer-width conversion.
- `src/libext/dftd4/nw_dftd4_params.F90`: read-only access to the parameters loaded from the pinned DFTD4 module.
- `src/libext/dftd4/build_dftd4.sh`: pinned source download and private CMake build.
- `src/libext/dftd4/GNUmakefile`: dependency orchestration and static bridge-archive build.
- `src/nwdft/xc/xc_dftd4.F`: functional mapping and energy, gradient, Bq, and Hessian handling.
- `src/nwdft/libxc/nwchem_libxc_util.F`: LibXC functional-name handoff.
- `src/nwdft/input_dft/dft_input.F`: `dftd4 on|off` parsing.
- Molecular DFT SCF, gradient, and task Hessian call sites contain the optional integration hooks.

## Development attribution

This interface was developed with assistance from OpenAI Codex. All AI generated code was reviewed, validated, and tested by MPH.
