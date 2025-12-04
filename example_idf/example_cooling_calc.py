"""
Example: Running F280 Cooling Load Calculation with AIM-2 Infiltration

This script demonstrates how to calculate cooling loads for a building using:
1. EnergyPlus IDF file for geometry
2. CSA F280-12 cooling load method
3. AIM-2 infiltration model for accurate air change calculations

The example uses the Montreal IDF file and shows how to:
- Override envelope properties
- Specify AIM-2 infiltration parameters
- Set climate conditions
- Extract and display results
"""

import sys
from pathlib import Path

# Add parent directory to path to import the F280 calculator
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from calculate_f280_cooling_sizing import f280_cooling_load_from_idf

# Path to the IDF file
idf_path = Path(__file__).parent / "CZ6_T1_408_Montreal.idf"

print("=" * 80)
print("F280 Cooling Load Calculation - Example with AIM-2 Infiltration")
print("=" * 80)
print(f"\nIDF File: {idf_path.name}")
print("\n" + "-" * 80)

# Run cooling load calculation with AIM-2 parameters
results = f280_cooling_load_from_idf(
    str(idf_path),
    
    # ========== CLIMATE & DESIGN CONDITIONS ==========
    latitude_deg=45.5,              # Montreal latitude
    longitude_deg=-73.6,            # Montreal longitude (auto-matches climate data)
    outdoor_cooling_temp=30.0,      # °C (design outdoor temperature)
    indoor_cooling_temp=24.0,       # °C (design indoor temperature)
    daily_temp_range=9.0,           # °C (for solar correction factors)
    
    # ========== ENVELOPE OVERRIDES (optional) ==========
    # Uncomment to override values parsed from IDF
    # ag_walls_rsi=3.5,             # m²·K/W - Above-grade wall insulation
    # attic_rsi=10.0,               # m²·K/W - Attic/roof insulation
    # windows_u=1.4,                # W/(m²·K) - Window U-value
    # shgc=0.35,                    # Solar heat gain coefficient
    
    # ========== AIM-2 INFILTRATION PARAMETERS ==========
    # These trigger the AIM-2 model instead of simple ACH/20 approximation
    n50=2.5,                        # ACH @ 50 Pa (from blower door test)
    aim2_flow_exponent=0.67,        # Flow exponent (0.65-0.75 typical)
    house_type='detached',          # 'detached' or 'semi-detached'
    foundation='full',              # 'crawl', 'slab', 'shallow', 'full'
    # storeys='1' is auto-detected from IDF, but can override: '1', '1.5', '2', '2.5', '3'
    
    # Wind conditions (optional - will load from NBC_F280_merge.csv if city provided)
    wind_speed_met=12.0,            # km/h at meteorological station
    terrain_class_met=3,            # Davenport terrain class at met station (1-8)
    terrain_class_site=7,           # Davenport terrain class at building site (1-8)
    shelter_walls=1.0,              # Wall shelter coefficient 0-1 (1=exposed, 0=fully sheltered)
    
    # Flue/chimney (if present)
    flue_diam_mm=0.0,               # mm (0 = no flue)
    shelter_flue=1.0,               # Flue shelter coefficient 0-1
    
    # ========== HVAC PARAMETERS ==========
    hrv_efficiency=0.6,             # Heat/Energy recovery ventilator efficiency (0-1)
    duct_location='conditioned',     # 'attic_open', 'crawl_enclosed', 'unconditioned_bsmt', 'slab_perimeter'
    duct_rsi=1.4,                   # Duct insulation RSI (m²·K/W)
    
    # ========== OTHER PARAMETERS ==========
    bedrooms=2,                     # Number of bedrooms (for occupant count if not in IDF)
    latent_multiplier=1.3,          # Latent load safety factor (CSA F280 default)
    shading_factor_global=1.0,      # Global shading factor 0-1 (1=no shading)
)

# ========== DISPLAY RESULTS ==========

print("\n" + "=" * 80)
print("RESULTS SUMMARY")
print("=" * 80)

# Building characteristics
print("\n--- Building Characteristics ---")
print(f"Floor Area:           {results['floor_area_total_m2']:.1f} m²")
print(f"Volume:               {results['volume_m3']:.1f} m³")
print(f"Wall Area (opaque):   {results['wall_area_opaque_m2']:.1f} m²")
print(f"Window Area:          {results['window_area_m2']:.1f} m²")
print(f"Window-to-Wall:       {results['window_area_m2'] / results['wall_area_gross_m2'] * 100:.1f}%")
print(f"Number of Occupants:  {results['num_occupants']}")

# Climate information
if 'climate_site_name' in results:
    print("\n--- Climate Data (Auto-Matched) ---")
    print(f"Site:                 {results['climate_site_name']}, {results['climate_site_province']}")
    print(f"Coordinates:          {results['climate_site_lat']:.2f}°N, {results['climate_site_lon']:.2f}°W")
    print(f"Design Cooling Temp:  {results['climate_site_dcdbt_C']:.1f}°C")
    print(f"Daily Temp Range:     {results['climate_site_daily_temp_range_C']:.1f}°C")
    print(f"Climate Zone:         {results['climate_zone_code']}")

# Infiltration details
print("\n--- Infiltration Analysis ---")
print(f"Method:               {results['infiltration_method']}")
print(f"Natural ACH:          {results['natural_ach_used']:.4f} h⁻¹")
print(f"Infiltration Load:    {results['HGsalb_infiltration_W']:.1f} W")

# Load components
print("\n--- Cooling Load Components (W) ---")
print(f"Opaque Envelope:      {results['HGcop_total_W']:>8.0f} W")
print(f"Transparent (Windows):{results['HGct_total_W']:>8.0f} W")
print(f"  - Solar Gains:      {results['q_solar_transmission_W']:>8.0f} W")
print(f"  - Conduction:       {results['q_windows_conduction_W']:>8.0f} W")
print(f"Internal Gains:")
print(f"  - People:           {results['HGsp_people_W']:>8.0f} W")
print(f"  - Appliances:       {results['HGapk_appliances_W']:>8.0f} W")
print(f"  - Laundry:          {results['HGapl_laundry_W']:>8.0f} W")
print(f"Infiltration:         {results['HGsalb_infiltration_W']:>8.0f} W")
print(f"Ventilation:          {results['HGsvb_ventilation_W']:>8.0f} W")
print(f"Duct Gains:           {results['HGdr_duct_gain_W']:>8.0f} W")
print(f"{'-' * 40}")
print(f"Total Sensible:       {results['HGsr_sensible_total_W']:>8.0f} W")
print(f"Latent (via multiplier): {(results['CSCn_nominal_W'] - results['HGsr_sensible_total_W']):>8.0f} W")
print(f"{'=' * 40}")
print(f"TOTAL COOLING LOAD:   {results['CSCn_nominal_W']:>8.0f} W")

# Equipment sizing
print("\n--- Equipment Sizing (CSA F280-12) ---")
print(f"Nominal Capacity:     {results['cooling_tons_nominal']:.2f} tons  ({results['CSCn_nominal_W']:.0f} W)")
print(f"Allowable Range:      {results['cooling_tons_min']:.2f} - {results['cooling_tons_max']:.2f} tons")
print(f"\n>>> SELECTED SIZE:    {results['cooling_tons_selected']:.2f} tons  ({results['cooling_kw_selected']:.2f} kW)")
print(f"Required Airflow:     {results['airflow_cfm']:.0f} CFM  ({results['airflow_ls']:.0f} L/s)")

# Solar details by orientation
print("\n--- Solar Gains by Orientation ---")
for orientation, details in results['solar_details'].items():
    print(f"{orientation.upper():>5}: {details['Solar_incident']:>6.0f} W/m² × {details['Area_m2']:>5.2f} m² × SHGC → "
          f"{details['Solar_incident'] * details['Area_m2'] * 0.6:>6.0f} W")

# Performance metrics
print("\n--- Performance Metrics ---")
print(f"Load Intensity:       {results['CSCn_nominal_W'] / results['floor_area_total_m2']:.1f} W/m²")
print(f"Sensible/Total:       {results['HGsr_sensible_total_W'] / results['CSCn_nominal_W'] * 100:.1f}%")
print(f"Latent Multiplier:    {results['latent_multiplier_used']:.2f}")
print(f"Latitude Factor:      {results['LFactor']:.3f}")
print(f"Duct Gain Multiplier: {results['duct_gain_multiplier']:.2f}")

print("\n" + "=" * 80)
print("Calculation complete!")
print("=" * 80)

# Optional: Save results to CSV
save_csv = False  # Set to True to save results
if save_csv:
    import pandas as pd
    output_path = Path(__file__).parent / "cooling_results.csv"
    df = pd.DataFrame([results])
    df.to_csv(output_path, index=False)
    print(f"\nResults saved to: {output_path}")
