# `This code is a work in progress. It likely has errors.`

# F280 Cooling Load Calculator with AIM-2 Infiltration

This repository contains tools for calculating residential cooling loads per CSA F280-12 standard, with advanced infiltration modeling using the AIM-2 (Air Infiltration Model) based on Bradley's HOT2000 implementation.

## Table of Contents
- [Overview](#overview)
- [AIM-2 Infiltration Model](#aim-2-infiltration-model)
- [F280 Cooling Load Calculator](#f280-cooling-load-calculator)
- [Integration: Using AIM-2 with F280](#integration-using-aim-2-with-f280)
- [Installation](#installation)
- [Examples](#examples)

---

## Overview

### What These Tools Do

**`AIM-2.py`**: Calculates natural air infiltration rates for residential buildings using physics-based modeling that accounts for:
- Stack effect (temperature-driven pressure differences)
- Wind pressure effects
- Leakage distribution (ceiling, floor, walls) based on building characteristics
- Building geometry, orientation, and site conditions

**`calculate_f280_cooling_sizing.py`**: Performs complete CSA F280-12 cooling load calculations including:
- Envelope heat gains (walls, roof, windows)
- Solar heat gains with orientation-specific adjustments
- Internal gains (occupants, appliances, lighting)
- Infiltration loads (using AIM-2 or simple approximation)
- Ventilation loads
- Duct gains
- Equipment sizing with CSA F280-12 compliance checks

---

## AIM-2 Infiltration Model

### What is AIM-2?

AIM-2 (Air Infiltration Model) calculates natural air change rates (ACH) based on building characteristics and environmental conditions. This implementation follows Bradley's HOT2000 version using L. Lew (1993) leakage distribution tables.

### Key Features

- **Physics-based**: Models both stack and wind effects
- **Building-specific**: Uses foundation type, house type, and number of storeys
- **Site-specific**: Accounts for terrain roughness and shelter
- **Temperature units**: Accepts Celsius (converts to Kelvin internally)

### How AIM-2 Works

1. **Leakage Coefficient (C)**: Calculated from blower door test result (ACH@50Pa)
2. **Leakage Distribution**: Assigns leakage fractions to ceiling/floor/walls based on:
   - House type (detached vs semi-detached)
   - Foundation type (crawl, slab, shallow, full)
   - Number of storeys (1, 1.5, 2, 2.5, 3)
3. **Stack Effect**: Calculates pressure difference from indoor/outdoor temperature
4. **Wind Effect**: Calculates pressure from wind speed, terrain, and shelter
5. **Superposition**: Combines stack and wind effects using Bradley's interaction term

### Using AIM-2 Standalone

#### Python Code Example

```python
from AIM_2 import compute_infiltration, Inputs

# Define building and conditions
aim2_inputs = Inputs(
    # Building characteristics
    volume=750.0,              # m³ (e.g., 300 m² × 2.5 m ceiling)
    n=0.67,                    # Flow exponent (typical: 0.65-0.75)
    n50=2.5,                   # ACH @ 50 Pa from blower door test
    
    # Geometry
    eave_height=5.0,           # m (2 storey × 2.5 m)
    
    # Temperatures (Celsius)
    indoor_temp=24.0,          # °C
    outdoor_temp=31.0,         # °C
    temp_unit='C',             # 'C' or 'K'
    
    # Building characteristics
    house_type='detached',     # 'detached' or 'semi-detached'
    foundation='crawl',        # 'crawl', 'slab', 'shallow', 'full'
    storeys='2',               # '1', '1.5', '2', '2.5', '3'
    
    # Flue (chimney)
    flue_diam_mm=0.0,         # mm (0 = no flue)
    
    # Wind and site
    wind_speed_met=12.0,      # km/h at met station
    met_height=10.0,          # m
    terrain_class_met=3,      # Davenport class (1-8)
    terrain_class_site=7,     # Davenport class at building
    shelter_walls=1.0,        # 0-1 (1=exposed, 0=sheltered)
    shelter_flue=1.0,         # 0-1 (for flue if present)
)

# Calculate infiltration
results = compute_infiltration(aim2_inputs)

# Display results
print(f"Natural ACH: {results['Qnat']:.3f}")
print(f"Stack component: {results['Qs']:.3f} ACH")
print(f"Wind component: {results['Qw']:.3f} ACH")
print(f"Leakage distribution - R:{results['R']:.3f}, X:{results['X']:.3f}, Y:{results['Y']:.3f}")
```

#### Command Line Example

```bash
python AIM-2.py \
  --volume 750 \
  --n 0.67 \
  --n50 2.5 \
  --eave_height 5.0 \
  --indoor_temp 24.0 \
  --outdoor_temp 31.0 \
  --temp_unit C \
  --house_type detached \
  --foundation crawl \
  --storeys 2 \
  --wind_speed_met 12.0 \
  --terrain_class_met 3 \
  --terrain_class_site 7
```

**Output:**
```
AIM-2 Infiltration Results (Bradley HOT2000 + L. Lew fractions)
---------------------------------------------------------------
C_base: 0.0123456
C_total: 0.0123456
Qnat: 0.0891
Qs: 0.0654
Qw: 0.0512
...
```

### AIM-2 Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `volume` | float | **required** | Building volume (m³) |
| `n` | float | 0.67 | Flow exponent (0.65-0.75 typical) |
| `n50` | float | **required** | ACH @ 50 Pa (from blower door) |
| `eave_height` | float | **required** | Height to eave (m) |
| `indoor_temp` | float | 22.0 | Indoor temperature |
| `outdoor_temp` | float | 0.0 | Outdoor temperature |
| `temp_unit` | str | 'C' | 'C' (Celsius) or 'K' (Kelvin) |
| `house_type` | str | 'detached' | 'detached' or 'semi-detached' |
| `foundation` | str | 'crawl' | 'crawl', 'slab', 'shallow', 'full' |
| `storeys` | str | '1' | '1', '1.5', '2', '2.5', '3' |
| `flue_diam_mm` | float | 0.0 | Flue diameter (mm), 0=no flue |
| `wind_speed_met` | float | 0.0 | Met station wind speed (km/h) |
| `met_height` | float | 10.0 | Met station height (m) |
| `terrain_class_met` | int | 3 | Davenport class at met (1-8) |
| `terrain_class_site` | int | 7 | Davenport class at site (1-8) |
| `shelter_walls` | float | 1.0 | Wall shelter coefficient (0-1) |
| `shelter_flue` | float | 1.0 | Flue shelter coefficient (0-1) |

**Davenport Terrain Classes:**
1. Open sea
2. Mud flats
3. Open flat terrain (grass)
4. Low crops
5. High crops
6. Parkland/bushes
7. Suburb/forest
8. City centre

---

## F280 Cooling Load Calculator

### What is CSA F280-12?

CSA F280-12 is the Canadian standard for determining heating and cooling loads in residential buildings for proper HVAC equipment sizing.

### Key Features

- **Complete cooling load calculation**: All CSA F280-12 components
- **IDF file parsing**: Extracts geometry from EnergyPlus models
- **Climate integration**: Looks up design conditions from NBC_F280_merge.csv or weather.csv
- **AIM-2 integration**: Automatic infiltration calculation when parameters provided
- **Orientation-specific**: Accounts for solar gains by direction
- **Equipment sizing**: Calculates required capacity with allowable ranges

### Calculation Components

#### Sensible Loads
1. **Opaque assemblies** (HGcop): Walls and roof with solar correction
2. **Transparent assemblies** (HGct): Windows with solar gains + conduction
3. **Internal gains**: Occupants, appliances, lighting
4. **Infiltration** (HGsalb): Via AIM-2 or simple ACH/20 approximation
5. **Ventilation** (HGsvb): Mechanical ventilation loads
6. **Duct gains** (HGdr): Based on location and insulation

#### Equipment Sizing
- **Latent multiplier** (default 1.3): Applied to sensible total per CSA F280-12 to account for latent loads (occupants, infiltration, ventilation moisture)
- Allowable capacity range per Clauses 6.3.2-6.3.5

### Using F280 Calculator

#### Python Code Example

```python
from calculate_f280_cooling_sizing import f280_cooling_load_from_idf

# Calculate cooling load from IDF file
results = f280_cooling_load_from_idf(
    'path/to/building.idf',
    
    # Climate/Design Conditions
    outdoor_cooling_temp=31.0,      # °C
    indoor_cooling_temp=24.0,       # °C
    daily_temp_range=12.0,          # °C
    latitude_deg=51.0,
    longitude_deg=-114.0,
    
    # Envelope Properties (override IDF)
    ag_walls_rsi=3.5,               # m²·K/W
    attic_rsi=10.0,                 # m²·K/W
    windows_u=1.4,                  # W/(m²·K)
    shgc=0.35,                      # Solar Heat Gain Coefficient
    
    # AIM-2 Infiltration Parameters
    n50=2.5,                        # ACH @ 50 Pa
    house_type='detached',
    foundation='full',
    city='Calgary',                 # Loads wind from weather.csv
    
    # HVAC
    hrv_efficiency=0.6,
    duct_location='attic_open',
    duct_rsi=1.4,
)

# Display results
print(f"Infiltration method: {results['infiltration_method']}")
print(f"Natural ACH: {results['natural_ach_used']:.3f}")
print(f"\n--- Load Components (W) ---")
print(f"Opaque envelope: {results['HGcop_total_W']:.0f}")
print(f"Windows: {results['HGct_total_W']:.0f}")
print(f"Solar transmission: {results['q_solar_transmission_W']:.0f}")
print(f"Internal gains: {results['HGsp_people_W'] + results['HGapk_appliances_W']:.0f}")
print(f"Infiltration: {results['HGsalb_infiltration_W']:.0f}")
print(f"Ventilation: {results['HGsvb_ventilation_W']:.0f}")
print(f"Duct gains: {results['HGdr_duct_gain_W']:.0f}")
print(f"\n--- Equipment Sizing ---")
print(f"Sensible total: {results['HGsr_sensible_total_W']:.0f} W")
print(f"Latent multiplier: {results['latent_multiplier_used']:.2f}")
print(f"Nominal capacity: {results['cooling_tons_nominal']:.2f} tons ({results['CSCn_nominal_W']:.0f} W)")
print(f"Selected: {results['cooling_tons_selected']:.2f} tons")
print(f"Capacity: {results['cooling_kw_selected']:.2f} kW")
print(f"Airflow: {results['airflow_cfm']:.0f} CFM")
```

#### Command Line Example

```bash
python calculate_f280_cooling_sizing.py \
  --idf building.idf \
  --latitude 51.0 \
  --longitude -114.0 \
  --toc 31.0 \
  --tic 24.0 \
  --dtr 12.0 \
  --ag-walls-rsi 3.5 \
  --attic-rsi 10.0 \
  --windows-u 1.4 \
  --shgc 0.35 \
  --duct-location attic_open \
  --duct-rsi 1.4 \
  --output results.csv
```

**Output:** JSON results printed to console and optionally saved to CSV.

### F280 Parameters Reference

#### Required
- `idf_path`: Path to EnergyPlus IDF file

#### Envelope (override IDF parsing)
- `ag_walls_rsi`: Above-grade wall RSI (m²·K/W)
- `attic_rsi`: Attic/roof RSI (m²·K/W)
- `windows_u`: Window U-value (W/(m²·K))
- `shgc`: Solar Heat Gain Coefficient (0-1)
- `ach`: ACH @ 50 Pa (default: 2.5)

#### Climate & Design
- `outdoor_cooling_temp`: Outdoor design temperature (°C)
- `indoor_cooling_temp`: Indoor design temperature (°C, default: 24)
- `daily_temp_range`: Daily temperature range for solar corrections
- `latitude_deg`: Site latitude (for solar calculations)
- `longitude_deg`: Site longitude (for climate lookup)
- `city`: City name (loads wind from weather.csv)

#### Infiltration
- `natural_ach`: Direct ACH override (bypasses AIM-2)
- `n50`: ACH @ 50 Pa (for AIM-2)
- `house_type`: 'detached' or 'semi-detached'
- `foundation`: 'crawl', 'slab', 'shallow', 'full'
- `aim2_flow_exponent`: Flow exponent n (default: 0.67)
- `flue_diam_mm`: Flue diameter (mm)
- `wind_speed_met`: Met station wind (km/h)
- `terrain_class_met`: Davenport class at met (1-8)
- `terrain_class_site`: Davenport class at site (1-8)
- `shelter_walls`: Wall shelter coefficient (0-1)
- `shelter_flue`: Flue shelter coefficient (0-1)

#### HVAC
- `PVC_l_s`: Principal ventilation rate (L/s)
- `hrv_efficiency`: Heat/Energy recovery ventilator efficiency (0-1)
- `duct_location`: 'attic_open', 'crawl_enclosed', 'unconditioned_bsmt', 'slab_perimeter'
- `duct_rsi`: Duct insulation RSI

#### Other
- `bedrooms`: Number of bedrooms (for occupant count)
- `latent_multiplier`: Latent load multiplier applied to sensible total (default: 1.3)
- `shading_factor_global`: Global window shading factor (0-1)
- `laundry_gain_W`: Laundry equipment gain (W)

---

## Integration: Using AIM-2 with F280

The F280 calculator automatically uses AIM-2 when infiltration parameters are provided. This gives more accurate results than the simple ACH50/20 approximation.

### Comparison Example

```python
from calculate_f280_cooling_sizing import calculate_cooling_load_f280

params = {
    'footprint': 150.0,
    'storeys': 2,
    'ag_walls_rsi': 3.5,
    'windows_u': 1.6,
    'attic_rsi': 10.0,
    'shgc': 0.35,
    'ach': 2.5,
}

# Method 1: Simple approximation (ACH50/20)
inputs_simple = {
    'outdoor_cooling_temp': 31.0,
    'indoor_cooling_temp': 24.0,
    'natural_ach': 2.5 / 20.0,  # Force simple method
}
results_simple = calculate_cooling_load_f280(params, inputs_simple)

# Method 2: AIM-2 model
inputs_aim2 = {
    'outdoor_cooling_temp': 31.0,
    'indoor_cooling_temp': 24.0,
    'n50': 2.5,
    'house_type': 'detached',
    'foundation': 'full',
    'city': 'Calgary',
}
results_aim2 = calculate_cooling_load_f280(params, inputs_aim2)

# Compare
print("SIMPLE METHOD:")
print(f"  Natural ACH: {results_simple['natural_ach_used']:.4f}")
print(f"  Infiltration: {results_simple['HGsalb_infiltration_W']:.0f} W")

print("\nAIM-2 METHOD:")
print(f"  Natural ACH: {results_aim2['natural_ach_used']:.4f}")
print(f"  Infiltration: {results_aim2['HGsalb_infiltration_W']:.0f} W")
print(f"  Method: {results_aim2['infiltration_method']}")

difference = (results_aim2['natural_ach_used'] / results_simple['natural_ach_used'] - 1) * 100
print(f"\nDifference: {difference:+.1f}%")
```

**Typical Output:**
```
SIMPLE METHOD:
  Natural ACH: 0.1250
  Infiltration: 219 W

AIM-2 METHOD:
  Natural ACH: 0.0913
  Infiltration: 160 W
  Method: AIM-2

Difference: -27.0%
```

### When to Use Each Method

**Use AIM-2 when:**
- You have blower door test results (ACH@50Pa)
- Building characteristics are known (type, foundation, storeys)
- More accurate sizing is required
- Accounting for specific site/wind conditions

**Use Simple method when:**
- Quick estimates needed
- Building details unknown
- Conservative estimates acceptable
- Already have measured/calculated natural ACH

---

## IDF File Requirements

The F280 calculator parses EnergyPlus IDF files to extract building geometry automatically. The parser is flexible and works with standard EnergyPlus objects - **no special naming conventions or modifications are required**.

### What Gets Parsed Automatically

**Building Geometry:**
- `BuildingSurface:Detailed` objects
  - **Walls**: Surfaces with type='Wall' and OutsideBoundaryCondition='Outdoors'
  - **Roofs**: Surfaces with type='Roof', 'RoofCeiling', or 'Ceiling'
  - **Floors**: Surfaces with type='Floor' (used to calculate footprint and storey count)
  - Wall orientations (North/East/South/West) calculated from surface normal vectors

**Windows:**
- `FenestrationSurface:Detailed` objects with type='Window', 'GlassDoor', or 'Skylight'
- `Window` objects (simplified window input)
  - Automatically linked to parent wall surface for orientation
  - Window orientations calculated from parent surface

**Glazing Properties:**
- `WindowMaterial:SimpleGlazingSystem` objects
  - Extracts U-value and SHGC
  - If multiple glazing systems exist, uses most conservative values (lowest U, highest SHGC)

**Occupancy:**
- `Zone` objects with floor area field (field 10)
- `People` objects with 'Area/Person' calculation method
  - Automatically calculates occupant count from zone floor area and area-per-person

### Zone Name Filtering

The parser **automatically excludes** certain zones from above-grade calculations to avoid double-counting:

**Excluded from wall area calculations:**
- Any zone with "basement" in the name (case-insensitive)
- Any zone with "attic" in the name (case-insensitive)

**Excluded from floor area calculations:**
- Any zone with "basement" in the name
- Any zone with "attic" in the name

This ensures that only conditioned living space is included in cooling load calculations.

### What You DON'T Need to Do

❌ **No special naming required** - zones can have any names
❌ **No specific construction names** - parser works with any construction
❌ **No coordinate system alignment** - orientations calculated automatically
❌ **No manual area calculations** - all computed from vertices

### What You CAN Override

Even with a complete IDF, you can override any parsed values:

```python
results = f280_cooling_load_from_idf(
    'building.idf',
    # Override parsed envelope properties
    ag_walls_rsi=5.0,      # Override wall insulation
    windows_u=1.2,         # Override glazing U-value
    shgc=0.3,              # Override solar heat gain coefficient
    # Add missing information
    ach=1.5,               # Add blower door test result
    latitude_deg=51.0,     # Add site location
)
```

### Fallback Behavior

If IDF parsing fails or returns incomplete data, the calculator uses estimation methods:

**Geometry estimation** (if IDF parsing incomplete):
- Requires: `footprint` (m²) and `storeys` parameters
- Estimates wall area from footprint perimeter × storey height
- Distributes windows equally across orientations

**Default envelope values** (if not in IDF or overrides):
- `ag_walls_rsi`: 3.5 m²·K/W
- `attic_rsi`: 8.8 m²·K/W  
- `windows_u`: 1.8 W/(m²·K)
- `shgc`: 0.4
- `ach`: 2.5 ACH@50Pa

### Example IDF Objects

**Typical BuildingSurface:Detailed** (no special requirements):
```
BuildingSurface:Detailed,
  Living_Room_Wall_North,     ! Name (any name is fine)
  Wall,                        ! Surface Type
  ExteriorWall,               ! Construction Name (any name is fine)
  Living_Room,                ! Zone Name (any name is fine)
  Outdoors,                   ! Outside Boundary Condition (must be 'Outdoors' for outdoor walls)
  ,                           ! Outside Boundary Condition Object
  SunExposed,                 ! Sun Exposure
  WindExposed,                ! Wind Exposure
  0.5,                        ! View Factor to Ground
  4,                          ! Number of Vertices
  0.0, 0.0, 2.5,             ! Vertex coordinates (parser calculates area and orientation)
  0.0, 0.0, 0.0,
  10.0, 0.0, 0.0,
  10.0, 0.0, 2.5;
```

**Typical Window** (simplified):
```
Window,
  Living_Room_Window_1,       ! Name (any name is fine)
  SimpleGlazing,              ! Construction Name
  Living_Room_Wall_North,     ! Building Surface Name (links to parent wall)
  ,                           ! Frame and Divider Name
  ,                           ! Multiplier
  2.0,                        ! Starting X Coordinate
  1.0,                        ! Starting Z Coordinate
  1.5,                        ! Length
  1.2;                        ! Height
```

### Tips for Best Results

✅ **Include SimpleGlazingSystem** - Provides accurate U-value and SHGC for windows
✅ **Use descriptive zone names** - Include "basement" or "attic" in names for automatic exclusion
✅ **Define People objects** - Enables automatic occupant count calculation
✅ **Include outdoor walls only** - Parser correctly ignores interior partitions
✅ **Standard EnergyPlus format** - Any valid EnergyPlus IDF will work

---

## Installation

### Requirements

```bash
pip install pandas numpy
```

### Files Needed

**Core:**
- `AIM-2.py` - Infiltration model
- `calculate_f280_cooling_sizing.py` - Cooling load calculator

**Optional Data:**
- `NBC_F280_merge.csv` - Climate design data (includes wind speeds, design temps, etc.)
  - Wind speed column: `JulWind` (July wind in km/h)
  - Design cooling temp: `DCDBT`
  - Daily temp range: `Strange`
- `f280_shading_factors.csv` - Custom shading overrides

---

## Examples

### Example 1: Quick Cooling Load Estimate

```python
from calculate_f280_cooling_sizing import f280_cooling_load_from_idf

results = f280_cooling_load_from_idf(
    'house.idf',
    outdoor_cooling_temp=30.0,
    indoor_cooling_temp=24.0,
)

print(f"Cooling required: {results['cooling_tons_selected']:.1f} tons")
```

### Example 2: Detailed Calculation with AIM-2

```python
results = f280_cooling_load_from_idf(
    'house.idf',
    # Climate
    city='Calgary',
    outdoor_cooling_temp=28.0,
    daily_temp_range=12.0,
    # Envelope
    ag_walls_rsi=5.0,  # High insulation
    attic_rsi=12.0,
    windows_u=1.2,     # Triple glazed
    shgc=0.3,
    # Infiltration (AIM-2)
    n50=1.5,           # Tight house
    house_type='detached',
    foundation='full',
    # HVAC
    hrv_efficiency=0.85,  # High efficiency HRV
)

print(f"Method: {results['infiltration_method']}")
print(f"ACH: {results['natural_ach_used']:.3f}")
print(f"Load: {results['cooling_tons_selected']:.2f} tons")
```

### Example 3: Sensitivity Analysis

```python
# Test different ACH50 values
for ach50 in [1.0, 1.5, 2.0, 2.5, 3.0]:
    results = f280_cooling_load_from_idf(
        'house.idf',
        n50=ach50,
        house_type='detached',
        foundation='full',
        city='Calgary',
    )
    print(f"ACH50={ach50}: Natural ACH={results['natural_ach_used']:.3f}, "
          f"Load={results['cooling_tons_selected']:.2f} tons")
```

### Example 4: Batch Processing

```python
import pandas as pd
from pathlib import Path

results_list = []

for idf_file in Path('buildings').glob('*.idf'):
    results = f280_cooling_load_from_idf(
        str(idf_file),
        n50=2.5,
        house_type='detached',
        foundation='full',
        city='Calgary',
    )
    results_list.append({
        'building': idf_file.stem,
        'cooling_tons': results['cooling_tons_selected'],
        'natural_ach': results['natural_ach_used'],
        'method': results['infiltration_method'],
    })

df = pd.DataFrame(results_list)
df.to_csv('cooling_results.csv', index=False)
print(df)
```

---

## Output Details

### AIM-2 Output Dictionary

```python
{
    'C_base': 0.0123,      # Base flow coefficient
    'C_total': 0.0123,     # Total flow coefficient (with flue)
    'C_flue': 0.0,         # Flue flow coefficient
    'Cc0': 0.0037,         # Ceiling leakage coefficient
    'Cf0': 0.0062,         # Floor leakage coefficient
    'Cw0': 0.0025,         # Wall leakage coefficient
    'R': 0.80,             # Leakage distribution parameter
    'X': -0.20,            # Leakage distribution parameter
    'Y': 0.00,             # Flue fraction
    'Ps': 2.45,            # Stack pressure (Pa)
    'fs': 0.67,            # Stack factor
    'Ue': 3.21,            # Effective wind speed (m/s)
    'Sw': 1.00,            # Shelter coefficient
    'Qs': 0.065,           # Stack flow (ACH)
    'Qw': 0.051,           # Wind flow (ACH)
    'Qnat': 0.089,         # Natural infiltration rate (ACH)
}
```

### F280 Output Dictionary (Selected Fields)

```python
{
    # Load components (W)
    'HGcop_total_W': 1234.5,
    'HGct_total_W': 2345.6,
    'HGsp_people_W': 280.0,
    'HGapk_appliances_W': 1200.0,
    'HGsalb_infiltration_W': 164.7,
    'HGsvb_ventilation_W': 125.0,
    'HGsr_sensible_total_W': 9271.9,
    
    # Equipment sizing
    'CSCn_nominal_W': 12053.5,
    'cooling_tons_selected': 3.5,
    'cooling_kw_selected': 12.31,
    'airflow_cfm': 1400.0,
    
    # Building info
    'floor_area_total_m2': 300.0,
    'volume_m3': 750.0,
    'num_occupants': 4,
    
    # Infiltration details
    'natural_ach_used': 0.0941,
    'infiltration_method': 'AIM-2',
}
```

---

## Tips & Best Practices

### Blower Door Testing
- Conduct test per ASTM E779 or CAN/CGSB-149.10
- Report as ACH@50Pa (not CFM@50Pa)
- Test during construction for tighter envelope goals

### AIM-2 Parameter Selection
- **Flow exponent (n)**: 0.67 is typical, 0.65 for very tight homes, 0.70-0.75 for older/leaky
- **Terrain classes**: Be realistic - most suburban homes are class 6-7
- **Shelter**: 1.0 (exposed) is conservative, reduce for dense neighborhoods

### F280 Sizing
- Review allowable capacity range (min/max)
- Consider oversizing penalties (humidity control, cycling)
- Account for future improvements (envelope upgrades)

### Climate Data
- Use local design conditions from ASHRAE or NBC
- Daily temperature range affects solar corrections significantly
- Wind speed: use seasonal average for cooling (July)

---

## Troubleshooting

### AIM-2 Issues

**Problem:** Results seem too low/high
- Check temperature units (Celsius vs Kelvin)
- Verify n50 is ACH@50Pa (not CFM@50Pa)
- Ensure storeys matches actual building (affects leakage distribution)

**Problem:** Import fails
- AIM-2.py has hyphen in filename - use importlib (see code examples)
- Check file is in same directory or Python path

### F280 Issues

**Problem:** Infiltration method shows 'Simple' instead of 'AIM-2'
- Verify AIM-2 parameters are provided (n50, house_type, foundation)
- Check AIM-2.py import succeeded
- Try providing explicit storeys parameter

**Problem:** Load seems unreasonably high/low
- Verify IDF geometry parsed correctly (check wall_area_gross_m2)
- Check design temperatures are reasonable for climate
- Review solar gains (shading_factor may be needed)
- Ensure RSI values are correct (not R-value in imperial units)

---

## References

- CSA F280-12: Determining the required capacity of residential space heating and cooling appliances
- L. Lew (1993): "Evaluation of AIM-2", Energy Mines and Resources Canada
- Bradley & Saunders: HOT2000 infiltration model implementation
- ASHRAE Handbook of Fundamentals: Infiltration and ventilation
- Davenport (1960): Terrain roughness classification

---

## License

See LICENSE file for details.

## Contact

For issues or questions about this implementation, please open an issue on the repository.
