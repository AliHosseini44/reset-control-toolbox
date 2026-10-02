# Reset Control Toolbox

MATLAB implementations of analysis and design tools for reset control systems, accompanying the PhD thesis.

This repository provides a compact toolbox and a worked example showing how the methods developed throughout the thesis can be applied to a common closed-loop reset-control problem. The example considers a mass plant with time delay and demonstrates controller design, frequency-domain analysis, stability assessment, nonlinearity quantification, and shaping-filter design.

## Main features

The repository includes tools for:

- CgLp controller design
- Higher-order sinusoidal-input describing functions (HOSIDFs)
- Closed-loop pseudo-sensitivity analysis
- Frequency-domain stability assessment
- Nonlinearity quantification using $\sigma_2$
- Reset-event reliability assessment using $\sigma_t$ and $\sigma_d$
- Shaping-filter design for reset control systems

## Repository structure


reset-control-toolbox/
│
├── README.md
|
│
├── examples/
│   └── Example_script.m
│
└── functions/
    ├── calc_omega_f.m
    ├── computeResetHOSIDF_SF.m
    ├── computePseudoSens_SF.m
    ├── computePseudoSens_SF_withMetrics.m
    ├── convertToLure_SF.m
    ├── theta_deg.m
    └── design_Cs_maxPhase.m


## MATLAB functions

### `calc_omega_f.m`

Calculates the CgLp filter frequency required to obtain a specified phase contribution at the desired bandwidth.

### `computeResetHOSIDF_SF.m`

Computes the higher-order sinusoidal-input describing functions of a reset element, including the effect of a shaping filter in the reset-triggering path.

### `computePseudoSens_SF.m`

Computes the closed-loop higher-order sinusoidal-input sensitivity functions and the pseudo-sensitivity magnitude of a zero-crossing reset control system.

### `computePseudoSens_SF_withMetrics.m`

Extends the pseudo-sensitivity calculation with the nonlinearity and reset-reliability metrics $\sigma_2$, $\sigma_t$, and $\sigma_d$.

### `convertToLure_SF.m`

Converts the considered closed-loop reset-control configuration into the Lur'e representation used by the frequency-domain analysis functions.

### `theta_deg.m`

Calculates the frequency-dependent angle $\theta_{\mathcal{N}}$ used in the frequency-domain stability condition for reset control systems.

### `design_Cs_maxPhase.m

Designs a shaping filter $C_s$ subject to the shaping-filter constraints developed in the thesis while maximizing the achievable phase contribution around a selected frequency.

## Worked example

The file


examples/Example_script.m


demonstrates the complete workflow.

The example includes:

1. Definition of a mass plant with time delay.
2. Design of a conventional PID controller as a linear baseline.
3. Design of an additional PI+CgLp controller.
4. Calculation of first- and higher-order reset HOSIDFs.
5. Calculation of the closed-loop pseudo-sensitivity.
6. Stability assessment using $\theta_{\mathcal{N}}$.
7. Evaluation of the nonlinearity metric $\sigma_2$.
8. Nonlinearity shaping using complementary pre- and post-filters.
9. Evaluation of the reliability metrics $\sigma_t$ and $\sigma_d$.
10. Design of a shaping filter $C_s$ and assessment of its effect on the reset behavior.

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

This repository accompanies the toolbox chapter of the PhD thesis and brings together methods developed in different chapters into one worked example.

The toolbox is intended primarily to illustrate how the proposed analysis and design tools can be used in practice. For the theoretical derivations, assumptions, and detailed interpretation of the methods, please refer to the corresponding chapters and publications.

## Attribution

The function `computePseudoSens_SF.m` is adapted from the work of L. F. van Eijk, D. Kostić, and S. H. HosseinNia on frequency-response analysis of general zero-crossing reset control systems.

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
