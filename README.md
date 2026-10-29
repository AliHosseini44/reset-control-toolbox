# Reset Control Toolbox

MATLAB implementations of analysis and design tools for reset control systems, accompanying the PhD thesis.

This repository provides a compact toolbox and a worked example showing how the methods developed throughout the thesis can be applied to a common closed-loop reset-control problem. The example considers a mass plant with time delay and demonstrates controller design, frequency-domain analysis, stability assessment, nonlinearity quantification, and shaping-filter design. The `ResetScope.m` function brings the main closed-loop analysis tools together in a single interface.

## Main features

The repository includes tools for:

- CgLp controller design
- Higher-order sinusoidal-input describing functions (HOSIDFs)
- Closed-loop pseudo-sensitivity and higher-order sinusoidal-input sensitivity functions (HOSISFs)
- Frequency-domain stability assessment
- Nonlinearity quantification using $\sigma_2$
- Reset-event reliability assessment using $\sigma_t$ and $\sigma_d$
- Shaping-filter design for reset control systems
- Unified closed-loop analysis through `ResetScope.m`

## Repository structure

```text
reset-control-toolbox/
│
├── README.md
│
├── examples/
│   └── Example_script.m
│
└── functions/
    ├── calc_omega_f.m
    ├── computeResetHOSIDF_SF.m
    ├── computePseudoSens_SF.m
    ├── ResetScope.m
    ├── convertToLure_SF.m
    ├── theta_deg.m
    └── design_Cs_maxPhase.m
```

## MATLAB functions

### `calc_omega_f.m`

Calculates the CgLp filter frequency required to obtain a specified phase contribution at the desired bandwidth.

### `computeResetHOSIDF_SF.m`

Computes the higher-order sinusoidal-input describing functions of a reset element, including the effect of a shaping filter in the reset-triggering path.

### `computePseudoSens_SF.m`

Computes the closed-loop HOSISFs, pseudo-sensitivity magnitude, and number of reset events for a zero-crossing reset control system. This is the **focused option** when sensitivity analysis is needed without the additional stability and performance/reliability metrics calculated by `ResetScope.m`.

### `ResetScope.m`

Provides a **unified interface** for closed-loop frequency-domain analysis of reset control systems. It computes:

- The pseudo-sensitivity magnitude and number of reset events per input period
- The closed-loop HOSISFs (`Swz`, `Swy`, and `Ln`)
- The nonlinearity metric $\sigma_2$
- The reset-trigger reliability metrics $\sigma_t$ and $\sigma_d$
- The stability angle $\theta_{\mathcal{N}}$ and a logical flag, `stable_reset`, indicating whether the implemented sufficient frequency-domain stability condition is satisfied

The integrated stability test is intended for the first-order reset-element setting considered in the thesis. A failed sufficient-condition test does **not** by itself establish instability.

**Feedthrough convention:** If the reset element's physical direct-feedthrough term is already represented by the parallel filter `frfCpar`, pass `Dr = 0` to `ResetScope` to avoid double-counting it, as done in the worked example.

### `convertToLure_SF.m`

Converts the considered closed-loop reset-control configuration into the Lur'e representation used by the frequency-domain analysis functions.

### `theta_deg.m`

Calculates the frequency-dependent angle $\theta_{\mathcal{N}}$ used in the frequency-domain stability condition for reset control systems. It remains available as a standalone utility, although the worked example obtains this angle directly from `ResetScope.m`.

### `design_Cs_maxPhase.m`

Designs a shaping filter $C_s$ subject to the shaping-filter constraints developed in the thesis while maximizing the achievable phase contribution around a selected frequency.

## Choosing the analysis function

Use `computePseudoSens_SF.m` when you only need pseudo-sensitivity, HOSISFs, or the number of reset events. Use `ResetScope.m` when you also need the nonlinearity/reliability metrics and the stability assessment. The worked example uses `ResetScope.m` to avoid repeating the same closed-loop analysis with several separate function calls.

For example, after defining the reset element, frequency grid, and component frequency responses:

```matlab
[pseudoMag,nrResets,Swz,Swy,Ln,sigma_2,sigma_d,sigma_t,theta_N,stable_reset] = ...
    ResetScope(f,sysR,0,A_rho,frfCs,frfCpre,frfCpar,frfCpos, ...
               frfPlant,Fs,nrHOSIDFsMax,overSampleRatio);
```

Here, `0` is used for `Dr` because the reset direct-feedthrough term is already included in `frfCpar`.

For sensitivity-only analysis, the alternative is:

```matlab
[pseudoMag,nrResets,Swz,Swy,Ln] = computePseudoSens_SF( ...
    f,sysR,A_rho,frfCs,frfCpre,frfCpar,frfCpos, ...
    frfPlant,Fs,nrHOSIDFsMax,overSampleRatio);
```

## Worked example

The file `examples/Example_script.m` demonstrates the complete workflow.

The example includes:

1. Definition of a mass plant with time delay.
2. Design of a conventional PID controller as a linear baseline.
3. Design of an additional PI+CgLp controller.
4. Calculation of first- and higher-order reset HOSIDFs.
5. Calculation of closed-loop pseudo-sensitivity and HOSISFs using `ResetScope.m`.
6. Frequency-domain stability assessment using $\theta_{\mathcal{N}}$ returned by `ResetScope.m`.
7. Evaluation of the nonlinearity metric $\sigma_2$ using `ResetScope.m`.
8. Nonlinearity shaping using complementary pre- and post-filters.
9. Evaluation of the reliability metrics $\sigma_t$ and $\sigma_d$ using `ResetScope.m`.
10. Design of a shaping filter $C_s$ and assessment of its effect on reset behavior.

The example also includes a commented call to `computePseudoSens_SF.m` for users interested only in sensitivity analysis.

## Requirements

The code was developed in MATLAB and uses:

- Control System Toolbox
- Optimization Toolbox

The Optimization Toolbox is required for `design_Cs_maxPhase.m`.

## How to run

Clone or download this repository and add the `functions` folder to the MATLAB path.

For example:

```matlab
addpath('functions')
```

Then run:

```matlab
run('examples/Example_script.m')
```

The example script uses the frequency grid

```matlab
f = 0.1:0.1:2000;
```

in Hz for the frequency-domain calculations.

## Relation to the thesis

This repository accompanies Chapter 9 (the MATLAB toolbox chapter) of the PhD thesis and brings together methods developed in different chapters into one worked example. In particular, `ResetScope.m` provides a common interface for closed-loop sensitivities, stability assessment, nonlinearity quantification, and reset-trigger reliability analysis.

The toolbox is intended primarily to illustrate how the proposed analysis and design tools can be used in practice. For the theoretical derivations, assumptions, and detailed interpretation of the methods, please refer to the corresponding chapters and publications.

## Attribution

The function `computePseudoSens_SF.m` is adapted from the work of L. F. van Eijk, D. Kostić, and S. H. HosseinNia on frequency-response analysis of general zero-crossing reset control systems. `ResetScope.m` builds on the same frequency-domain sensitivity-analysis foundation and integrates the additional metrics and stability assessment presented in the thesis.

The original work also provides an implementation of the pseudo-sensitivity calculation for reset-control configurations without the shaping-filter extension used in this repository.

Users of this toolbox are encouraged to cite both the corresponding original work and the relevant thesis/publications when using these functions.

## References

1. A. Hosseini, PhD thesis,  
   Delft University of Technology, 2026.

2. L. F. van Eijk, D. Kostić, and S. H. HosseinNia,  
   *Frequency Response Analysis of General Zero-Crossing Reset Control Systems*,  
   IEEE Control Systems Letters, 2025.

Additional publications associated with the individual methods are listed in the thesis.

## Citation

If you use this toolbox in academic work, please cite the corresponding thesis and the relevant publications associated with the method being used.

A `CITATION.cff` file can also be included in this repository to provide citation information directly through GitHub.

## License

This repository is distributed under the license provided in the `LICENSE` file.

## Contact

For questions regarding the toolbox or the methods implemented here, please contact:

**Ali Hosseini**  
Delft University of Technology  
s.a.hosseini@tudelft.nl
