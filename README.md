# Differential-Evolution-based Spectral Correction (DESC) Method for Earthquake Stress Drop Estimation

Author: Jiewen Zhang (jayxps.work@gmail.com)

This document describes the complete workflow for estimating earthquake stress drops using the DESC package.

Note that all the output files are saved as MATLAB .mat binary format, which can be easily accessed in both Python and MATLAB.

During the program execution, you may notice some dependencies/packages/modules are missing. Please follow some common pipelines on the Internet to install them.

---

## Overview

The DESC method estimates earthquake stress drops by combining spectral decomposition with differential evolution optimization. The workflow consists of three main stages:

1. **Data Preparation** - Convert SAC files to MATLAB format and compute spectra
2. **Spectral Decomposition** - Separate observed spectra into source, site, and path effects
3. **Stress Drop Estimation** - Estimate stress drops using the DESC method


Required Python packages:
- NumPy
- SciPy
- ObsPy
- Matplotlib
- JobLib (For parallel processing)
- Multitaper (Credit: Germán A. Prieto, https://github.com/gaprieto/multitaper)


References:

Zhang J, Yang H, Zi J, et al. An improved estimation of stress drop and its application on induced earthquakes in the Weiyuan Shale Gas Field in China. Geophysical Journal International, 2024, 236(3): 1785-1803.

Zhang J, Yang H. Improved source parameter estimation of earthquakes in the 2019 Ridgecrest sequence based on a global‐optimization algorithm and their implications on fault behaviors. Bulletin of the Seismological Society of America, 2025, 115(3): 805-824.

Jiang Y, Zhang J, Zi J, et al. Source mechanism and rupture directivity of small earthquakes in the Changning region, China, using a dense array data. arXiv preprint arXiv:2512.11326, 2025.

Prieto G A. The multitaper spectrum analysis package in Python. Seismological Society of America, 2022, 93(3): 1922-1929.

Storn R, Price K. Differential evolution–a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 1997, 11(4): 341-359.


Contributors:
Jiewen Zhang, Xiaowei Chen, Rachel E. Abercrombie, Hongfeng Yang, Youjie Jiang

Special thanks to Prof. Zhigang Peng who provided the module reading SAC files.

---

## Directory Structure

**Note 1**: the Python and MATLAB packages follow the same naming patterns other than the suffixes
- Python: .ipynb (main programs) and .py (subroutines) for Python (Tested for Python 3.8 and higher)
- MATLAB: .m for MATLAB (Tested for MATLAB 2020 and higher)

**Note 2**: You are recommended to use the parallel version to accelerate preprocessing in case of a large dataset.

**Note 3**: If you compare the absolute stress drops calculated using the Python and MATLAB versions, you may notice very minor difference; this is because when calculating the spectra, the Python version is using a different package than the MATLAB one, which introduce slight differences in parameter selections. Such difference is too minor compared to the uncertainty and variability levels.

```
DESC_Python/
├── data/
│   ├── sac/                          # Raw seismic waveforms (SAC format) (Make sure all events have 3C data)
│   ├── matfile/                      # Event data in MATLAB format
│   └── matfile_spec/                 # Event data with computed spectra
├── preprocessing/
│   ├── catalog                       # Event catalog
│   ├── station.txt                   # Station information
│   ├── specinp                       # Processing parameters
│   ├── sac_mat.ipynb                 # SAC → MAT conversion
│   ├── sac_mat_par.ipynb             # SAC → MAT conversion with parallel processing
│   ├── compute_spec.ipynb            # Spectrum computation
│   └── compute_spec_par.ipynb        # Spectrum computation with parallel processing
├── analysis/
│   ├── spectral_decomposition.ipynb  # Spectral decomposition
│   ├── DESC_inversion.ipynb          # Stress drop estimation
│   ├── evspec_*.mat                  # Event spectra from decomposition
│   ├── P.mat                         # Processing parameters
│   └── maginfo.mat                   # Magnitude calibration info
└── src/
    ├── update_eqinfo.py               # Update event catalog
    ├── update_data_spec.py            # Compute spectra
    ├── specprocess_*.py               # Processing modules
    ├── sourcepara.py                  # Source parameter estimation
    └── DE_Module/                     # Differential Evolution optimizer
```

A quick run using samples provided:
- **Step 1**: preprocessing/sac_mat(_par)
- **Step 2**: preprocessing/compute_spec(_par)
- **Step 3**: analysis/spectral_decomposition
- **Step 4**: analysis/DESC_inversion

---

## Stage 1: Data Preparation

### Step 1.1: Prepare Event Catalog

**Input:**
- Event origin times, locations, magnitudes

**Process:**
Build yourself a catalog file (`preprocessing/catalog`) following the format below:
```
YEAR  MONTH  DAY  HOUR  MINUTE  SECOND  MAGNITUDE  LATITUDE  LONGTIUDE  DEPTH  ID
2019     11    1     1      58   22.99       0.30   29.6120   -95.2213   2.91   1
2019     11    1    23      23   57.89       0.40   29.6150   -95.2252   3.18   2
...
```

### Step 1.2: Prepare Station Information

**Input:**
- Station locations

**Process:**
Build yourself a station list `preprocessing/station.txt` following the format below:
```
Net_Sta_Name Latitude Longitude
     YXYX001  29.5000  -95.2000
     YXYX002  29.6000  -95.3000
 
    NET: YX
    STATION: YX001

...
```

### Step 1.3: Organize SAC Files

**Input:**
- Seismic waveforms in SAC format

**Structure:**
```
data/sac/
├── 1/          # Event ID 1
│   ├── {SAC file for YX001}
│   ├── {SAC file for YX002}
│   └── ...
├── 2/          # Event ID 2
│   ├── {SAC file for YX001}
│   ├── {SAC file for YX002}
│   └── ...
```

### Step 1.4: Convert SAC to MATLAB Format

#### Note: Inspect the SAC headers in the sample dataset first, and match your own with them

**Notebook:** `preprocessing/sac_mat.ipynb`
              `preprocessing/sac_mat_par.ipynb`
**MATLAB**    `preprocessing/sac_mat.m`
              `preprocessing/sac_mat_par.m`

**Process:**
1. Read catalog and station information
2. For each event:
   - Load SAC files for all stations
   - Extract waveform data and headers
   - Calculate event-station distances
   - Save to `data/matfile/YYYY/MM/evXXXXX.mat`

**Output:**
- `data/matfile/eqinfo.mat` - Event catalog in MATLAB format
- `data/matfile/YYYY/MM/evXXXXX.mat` - Waveform data and metadata for each event

### Step 1.5: Compute Spectra

**Notebook:** `preprocessing/compute_spec.ipynb`
              `preprocessing/compute_spec_par.ipynb`
**MATLAB:**    `preprocessing/compute_spec.m`
              `preprocessing/compute_spec_par.m`

**Process:**
1. Load `eqinfo.mat`
2. For each event:
   - Read waveform data
   - Compute displacement spectra using multitaper method
   - Save the displacement spectra along with the original data file

**Input data:**
- `data/matfile/YYYY/MM/evXXXXX.mat` - Event data

**Output data:**
- `data/matfile_spec/YYYY/MM/evXXXXX.mat` - Event data with spectra

---

## Stage 2: Spectral Decomposition

**Notebook:** `analysis/spectral_decomposition.ipynb`
**MATLAB:**    `analysis/spectral_decomposition.m`

**Goal:** Separate observed spectra into source, site, and path components.

### Step 2.1: Configure Processing Parameters

**File:** `preprocessing/specinp`

### Step 2.2: Data Organization and Filtering

**Process:**
1. Read parameters from `specinp`
2. Load all event spectra from `data/matfile_spec/`
3. Apply criteria in `specinp`

**Key Functions:**
- `specprocess_readpara()` - Read configuration
- `specprocess_organize()` - Organize and filter data

**Output:**
- `spec.mat`: log10 Spectra for spectral decomposition (all the spectra below are log10 converted)

### Step 2.3: Iterative Spectral Inversion

**Mathematical Model:**

For each spectrum observation:
```
log10(A_ijk) = log10(E_i) + log10(S_j) + log10(D_k) + ε
```

Where:
- `A_ijk` = Observed spectral amplitude
- `E_i` = Event term (source spectrum)
- `S_j` = Station term (site amplification)
- `D_k` = Distance term (geometric spreading)
- `ε` = Residual error

**Process:**
1. Initialize all terms to zero
2. Iteratively solve for event terms, station terms and path (distance) terms

**Key Functions:**
- `specprocess_inversion()` - Main iterative solver

**Output:**
- `evspec` - Event spectra (source term)
- `stspec` - Station spectra (site term)
- `distspec` - Distance spectra (path term)
- `freq` - Frequency array

### Step 2.4: Quality Control

**Filters Applied:**
1. Minimum station coverage: Events with ≥ `nspecmin` stations (refer to Shearer et al., 2006)
2. Spectral quality: Remove events with flat/unrealistic spectra
3. Magnitude consistency: Events with reasonable magnitude estimates

### Step 2.5: Magnitude Calibration

**Key Functions:**
- `specprocess_specmag()` - Magnitude calculation
- `specprocess_plotmag()` - Diagnostic plots

**Output:**
- `maginfo` - Magnitude calibration information

### Step 2.6: Save Decomposition Results

**Files Saved:**
- `evspec_XXXXX.mat` - Event spectra and metadata (case suffix: XXXXX)
- `spec_XXXXX.mat` - All decomposed spectra
- `P.mat` - Processing parameters
- `maginfo.mat` - Magnitude calibration

**Event Spectrum Fields:**
- `qid` - Event ID
- `qtime` - Origin time (epoch time)
- `qlat`, `qlon`, `qdep` - Location
- `qmag` - Catalog magnitude
- `qmagest` - Estimated moment magnitude
- `qmomest` - Estimated seismic moment
- `spec` - Source spectrum
- `nspec` - Number of stations used
- `qmagresid` - Magnitude residual

---

## Stage 3: Stress Drop Inversion

**Notebook:** `analysis/DESC_inversion.ipynb`
**MATLAB:**   `analysis/DESC_inversion.m`

**Goal:** Estimate corner frequencies and stress drops using floating EGF optimization.

### Step 3.1: Load Decomposed Spectra

**Input Files:**
- `evspec_XXXXX.mat` - Event spectra from Stage 2
- `P.mat` - Processing parameters
- `maginfo.mat` - Magnitude calibration

**Data Loaded:**
- `evspec` - List of event spectra (source terms)
- `freq` - Frequency array
- `freqlim` - Calibration frequency range

### Step 3.2: Configure Inversion Parameters

**Source Models:**
```python
prefix_all = ['brune', 'boatwright', 'doublecorner']
method_all = [1, 5, 11]
```

**Model Descriptions:**
1. **Brune model** (`imethod=1`): 
   - Single corner frequency
   - `Ω(f) = Ω₀ / (1 + (f/fc)²)`
   
2. **Boatwright model** (`imethod=5`):
   - Single corner frequency
   - `Ω(f) = Ω₀ / √(1 + (f/fc)²)²`
   
3. **Double-corner** (`imethod=11`):
   - Two corner frequencies
   - `Ω(f) = Ω₀ / [√(1 + (f/fc1)²) × √(1 + (f/fc2)²)]`

**Inversion Parameters (change accordingly):**
```python
maglim = [0.9, 4.0]    # Magnitude range for inversion
d_mag = 0.2            # Magnitude bin size
ddsig = [0.01, 100]    # Stress drop search range (MPa)
min_nspec = 10         # Minimum events per magnitude bin
freqfit = [2, 40]      # Frequency band for fitting (Hz)
iwave = 1              # Wave type (1=P, 2=S)
```

### Step 3.3: Stack Spectra by Magnitude

**Process:**
1. Divide events into magnitude bins (width = `d_mag`)
2. For each bin:
   - Compute median spectrum
   - Calculate geometric mean of seismic moments
   - Count number of events

**Key Functions:**
- `specprocess_stackmag()` - Stack spectra by magnitude

**Output:**
- `stack_spec` - Stacked spectra for each magnitude bin
  - `mag` - Bin center magnitude
  - `spec` - Median spectrum
  - `nspec` - Number of events
  - `freq` - Frequency array

### Step 3.4: ECS Inversion

**Method:** Differential Evolution optimization to find stress drops that minimize ECS variability across magnitude bins.

**Theory:**

For each event with moment `M₀` and stress drop `Δσ`:

1. **Corner frequency from stress drop:**
   ```
   fc = (Δσ × fact / M₀)^(1/3)
   ```
   where `fact = (κ × β)³` and `κ` depends on wave type

2. **Predicted source spectrum:**
   ```
   S(f) = Ω₀ × Model(f, fc)
   ```

3. **Empirical Correction Spectrum (ECS):**
   ```
   ECS = log10(Observed_spectrum) - log10(Predicted_source_spectrum)
   ```

**Optimization Goal:**

Find stress drop values for each magnitude bin that minimize:
```
Cost = Σ std(ECS across magnitude bins at each frequency)
```

This ensures ECS is consistent across all events (check Zhang et al., 2024 for further details).

**Process:**

1. **Assign events to magnitude bins**
   - Filter events in magnitude range `[maglim[0], maglim[1]]`
   - Assign to nearest bin

2. **Set up Differential Evolution (change accordingly):**
   - Decision variables: `log10(Δσ)` for each magnitude bin
   - Population size: 25
   - Max iterations: 1000
   - Strategy: DE/rand/1 with per-generation dither
   - Bounds: `[log10(0.01), log10(100)]` MPa

3. **For each DE iteration:**
   - Calculate corner frequencies from stress drops
   - Compute predicted spectra using source model
   - Calculate ECS for all events in diffrent magnitude bins
   - Compute median ECS for each magnitude bin
   - Return median ECS standard deviation across bins as cost

4. **Extract final results:**
   - Optimized stress drop for each bin using DE
   - A uniform ECS = median of all ECSs

**Key Functions:**
- `specprocess_floategf_allEQs()` - Main EGF optimization
- `Rundeopt()` - Differential Evolution setup
- `deopt()` - DE optimization loop
- `objfun()` - Objective function (cost calculation)

**Output:**
- `ECS` - Empirical Calibration Spectra for all magnitude bins
- `log10dsigbins` - Optimized log10(stress drop) for each bin
- Updated `stack_spec` with:
  - `dsig` - Stress drop (MPa)
  - `fc` - Corner frequency (Hz)

### Step 3.5: Compute Individual Event Parameters

**Process:**

For each event:
1. Subtract ECS from observed spectrum
2. Fit source model to corrected spectrum
3. Extract corner frequency `fc`
4. Calculate source parameters:

**Source Parameters:**

```
Corner frequency:    fc (Hz)
Stress drop:         Δσ = fc³ × M₀ / fact / 10⁶ (MPa)
Source radius:       r = κ × β / fc (m)
Seismic moment:      M₀ (N·m)
Radiated energy:     Es (J) - integrated from spectrum
Apparent stress:     σₐ = μ × Es / M₀ / 10⁶ (MPa)
Rupture area:        A = π × r² (m²)
Average slip:        D = M₀ / (μ × A) (m)
Fracture energy:     G = 0.5 × (Δσ - 2σₐ) × D × 10⁶ (J/m²)
```

Physical constants used:
- `β = 3464 m/s` - Shear wave velocity
- `ρ = 2700 kg/m³` - Density
- `μ = ρ × β²` - Shear modulus
- `κ = 0.42` (P-wave) or `0.2766` (S-wave) - Shape factor

**Corner Frequency Uncertainty:**

Test corner frequency uncertainties:
- For each test frequency, re-optimize other parameters
- Find 5% error bounds around minimum
- Report `[fcb1, fcb2]` as uncertainty range

**Key Functions:**
- `specprocess_sourcepara()` - Fit individual event spectrum for all events

**Output:**
- `para_EGF` - List of parameter dictionaries for each event
  - `fc`, `fc2` - Corner frequencies
  - `delsig` - Stress drop
  - `asig` - Apparent stress
  - `r` - Source radius
  - `Es` - Radiated energy
  - `G` - Fracture energy
  - `err` - Fitting error
  - `fcb1`, `fcb2` - Corner frequency bounds

### Step 3.6: Save Results

**Output Directory Structure:**
```
analysis/
├── result-all-[model]-f[flow]-[fhigh]/
│   ├── para_ECS_[model].mat          # All results
│   ├── stacked-spectra-after-EGF.jpg # Diagnostic plot
│   └── *.pdf                          # Individual event fits
└── result-indieq-[model]-f[flow]-[fhigh]/
    └── result-constant-egf.txt        # Parameter table
```

**Saved Files:**

1. **para_ECS_[model].mat** - Main results file containing:
   - `para_EGF` - Individual event parameters
   - `stack_spec` - Magnitude-binned results
   - `ECS` - Empirical Calibration Spectra
   - `evspec`, `stspec`, `distspec` - Decomposed spectra

2. **result-constant-egf.txt** - Human-readable table:
   ```
   i    Mw_est  Mw_cat  Event_ID   Ω₀    fc     r      Δσ    σₐ    G      Es
   1    2.50    2.45    12345      12.3  5.2    200m   2.5   1.2   5.4e5  1.2e8
   ...
   ```

3. **Diagnostic plots:**
   - Stacked spectra with corner frequencies marked
   - Individual event spectral fits (if enabled)

---

## Key Algorithms and Methods

### Multitaper Spectrum Estimation

**Purpose:** Robust spectral estimation with reduced variance

**Parameters:**
- `nw = 2.5` - Time-bandwidth product
- `kspec = 3` - Number of tapers
- `Nfft = 1024` - FFT length

**Process:**
1. Apply orthogonal Slepian tapers to signal
2. Compute FFT for each taper
3. Average across tapers for final spectrum

### Differential Evolution (DE) Optimizer

**Purpose:** Global optimization for stress drop estimation

**Strategy:** DE/rand/1 with per-generation dither

**Algorithm:**
```
1. Initialize population of stress drop vectors
2. For each generation:
   a. For each population member x:
      - Select random members: a, b, c
      - Create mutant: v = c + F × (a - b)
      - Crossover with x to create trial vector u
      - Evaluate cost of u
      - If cost(u) < cost(x), replace x with u
   b. Track best solution
3. Return best stress drop values
```

**Parameters (change if not converging, read Storn and Price, 1997 for best practice):**
- Population: 25 members
- F_weight: 0.6 (mutation scale)
- F_CR: 0.9 (crossover probability)
- Max iterations: 1000

**Cost Function:**
```
Cost = 1 / Σ_f std(ECS across magnitude bins at frequency f)
```

Maximizes consistency of EGF across all events.

---

## Best Practices and Recommendations

### Data Quality

1. **Station Coverage:**
   - Minimum stations per event should be tested for result stability
   - Good azimuthal coverage reduces bias
   - Similar station distribution across magnitude range

2. **Frequency Range:**
   - Low frequency: Choose carefully to avoid low-freq problems e.g. drifting
   - High frequency: Choiise carefully to avoid instrument upper limit
   - Fitting band: Choose based on SNR and corner frequency range

3. **SNR Threshold:**
   - Minimum SNRs for frequency bands should be tested to balance between data quality and amount

### Processing Parameters

1. **Magnitude Binning:**
   - Bin width: 0.2-0.3 magnitude units
   - Minimum events per bin: 10+ for stable median
   - Wider range and more M>1 earthquakes provides better ECS constraint

2. **Stress Drop Range:**
   - Search range: 0.01-100 MPa should be sufficient for most cases
   - Very high/low values may indicate poor fit

3. **Iteration Counts:**
   - Spectral inversion: 20 iterations usually sufficient
   - DE optimization: populations and iterations for convergence should be tested for consistency
   - Monitor convergence in output

### Quality Checks

1. **After Spectral Decomposition:**
   - Check magnitude calibration plots
   - Verify station terms are reasonable (typically -0.5 to +0.5 in log)
   - Examine distance dependence (should be smooth)

2. **After Stress Drop Inversion:**
   - Verify corner frequencies are in fitting band
   - Check if stress drop levels are reasonable
   - Examine spectral fit residuals
   - Compare with scaling relations (e.g., fc vs. M₀)

3. **Physical Constraints:**
   - Stress drops: 0.01-100 MPa (extreme: 0.001-1000)

---

## Troubleshooting (Note: If you don't understand some certain parts, avoid making any changes)

### Common Issues

1. **Poor magnitude calibration:**
   - Check frequency range for moment estimation
   - Verify station terms are stable
   - Ensure sufficient low-frequency data

2. **High stress drop scatter:**
   - Increase `min_nspec` (more events per bin)
   - Tighten quality filters (SNR, station count, etc.)

3. **DE optimization not converging:**
   - Increase `I_NP` (population per iteration)
   - Increase `I_itermax` (max iterations)
   - Adjust stress drop search range
   - Check for data quality issues

4. **Unrealistic corner frequencies:**
   - Verify fitting frequency range
   - Check ECS subtraction is working if you accidentally changed something
   - Examine individual spectral fits

### Performance Notes

The package has been optimized for computational efficiency:
  
- **Vectorized operations** in DE objective function
  - Replaced Python loops with NumPy array operations
  - Major performance gain in innermost optimization loop

- **Parallel processing** in data preprocessing
  - Parallel processing modules are ready for SAC-to-MAT conversion and spectra calculation
  - Number of cores can be specified

For large datasets (>1000 events):
- Consider parallel processing of individual events
- Monitor memory usage with very dense station networks
- Use subset for testing before full processing

---

## Source Model Equations

**Brune (1970):**
```
Ω(f) = Ω₀ / (1 + (f/fc)²)
fc = 0.42 × β × (Δσ / M₀)^(1/3)
```

**Boatwright (1978):**
```
Ω(f) = Ω₀ / √(1 + (f/fc)²)²
fc = 0.2766 × β × (Δσ / M₀)^(1/3)
```

**Stress Drop:**
```
Δσ = (7/16) × M₀ / r³  (Brune)
```

---

## Summary Checklist

### Stage 1: Data Preparation
- [ ] Create event catalog
- [ ] Organize SAC files by event ID
- [ ] Convert SAC to MAT format
- [ ] Compute displacement spectra
- [ ] Check SNR and data quality

### Stage 2: Spectral Decomposition
- [ ] Configure `specinp` parameters
- [ ] Run spectral inversion
- [ ] Filter by station count
- [ ] Calibrate magnitudes
- [ ] Save `evspec_*.mat`, `P.mat`, `maginfo.mat`
- [ ] Verify magnitude calibration plots

### Stage 3: Stress Drop Inversion
- [ ] Load decomposed spectra
- [ ] Configure inversion parameters
- [ ] Stack spectra by magnitude
- [ ] Run floating ECS optimization
- [ ] Compute individual event parameters
- [ ] Save results and plots
- [ ] Check stress drop scaling relations

### Quality Control
- [ ] Magnitude residuals < 0.5
- [ ] Corner frequencies in fitting band
- [ ] Stress drops physically reasonable
- [ ] Good spectral fit quality

---

## Contact and Support

For questions about the DESC method or this package implementation, refer to the original publication or contact me.
