"""CSA F280-12 Cooling Equipment Sizing (Single-IDF Mode)

This module provides a single-building cooling load calculation per CSA F280-12 Section 6.
Primary entry point: `f280_cooling_load_from_idf(idf_path, **overrides)` or CLI via `_run_cli()`.

Implemented Components (Sensible):
    - Opaque assemblies conduction + solar correction (Table 3)
    - Transparent assemblies solar + conduction (Solar0 + SHGC)
    - Internal gains (occupants, appliances, laundry)
    - Infiltration & ventilation loads (simplified ACH method / principal ventilation L/s)
    - Duct gains (location + insulation multiplier)
    - Equipment sizing & allowable range (Clauses 6.3.1–6.3.5)
    - Latent loads accounted for via latent_multiplier (default 1.3)

Climate Integration:
    - Optional nearest-site lookup from `NBC_F280_merge.csv` when latitude & longitude supplied.
    - Populates design cooling dry bulb and daily temperature range ("Strange" column) if not overridden.

Geometry Extraction (Lightweight):
    - Parses IDF surfaces for wall, window, roof, and floor areas + orientation distribution.
    - Extracts Simple Glazing U-value & SHGC; estimates occupant count if People objects present.

Notes:
    - Legacy multi-building batch processing has been removed for clarity.
    - Shading factors may be globally overridden (`shading_factor_global`) or via optional CSV `f280_shading_factors.csv`.
    - Output dictionary includes detailed component breakdown and climate metadata when available.
"""

import pandas as pd
import numpy as np
import math
import re
import json
from pathlib import Path
import warnings
import importlib.util
import sys

# Import AIM-2 infiltration model
try:
    # Load AIM-2.py (with hyphen in filename)
    aim2_path = Path(__file__).parent / "AIM-2.py"
    if aim2_path.exists():
        spec = importlib.util.spec_from_file_location("aim2_module", aim2_path)
        aim2_module = importlib.util.module_from_spec(spec)
        sys.modules["aim2_module"] = aim2_module
        spec.loader.exec_module(aim2_module)
        aim2_compute_infiltration = aim2_module.compute_infiltration
        AIM2Inputs = aim2_module.Inputs
    else:
        aim2_compute_infiltration = None
        AIM2Inputs = None
except Exception:
    # Fallback if import fails
    aim2_compute_infiltration = None
    AIM2Inputs = None

# (Legacy multi-building batch globals removed)


# CSA F280-12 Constants
W_PER_TON = 3517.0  # Watts per ton of refrigeration
CFM_PER_TON = 400.0  # Standard airflow: 400 CFM per ton
CFM_TO_M3S = 0.0004719474  # Conversion factor: CFM to m³/s
BTU_TO_W = 0.2930711  # Conversion factor: BTU/h to W

# Cooling design conditions (typical Canadian summer design)
OUTDOOR_TEMP_DESIGN = 31.0  # °C (design outdoor temperature)
INDOOR_TEMP_DESIGN = 24.0  # °C (design indoor temperature)
OUTDOOR_RH = 0.60  # Relative humidity (fraction)
INDOOR_RH = 0.50  # Relative humidity (fraction)

###############################################
# CSA F280-12 TABLE DATA (EXCERPTS IMPLEMENTED)
###############################################

# Clause 6.2.2.1 Estimated solar radiation base (Solar0) W/m²
# Directions grouped per standard (East/West share value; NE/NW; SE/SW; Horizontal for sloped glazing/skylights)
SOLAR0_BASE = {
    'north': 93,
    'south': 160,
    'east': 285,
    'west': 285,
    'northeast': 194,
    'northwest': 194,
    'southeast': 252,
    'southwest': 252,
    'horizontal': 534
}

# Solar correction (SC) °C for opaque assemblies (Table 3) keyed by orientation and daily temperature range class
# We store tuples: (SC for DTR <=14°C, SC for DTR >14°C)
SOLAR_CORRECTION_TABLE = {
    'north': (-3.8, -6.7),
    'northeast': (-0.6, -3.3),
    'northwest': (-0.6, -3.3),
    'east': (1.7, -1.1),
    'west': (1.7, -1.1),
    'southeast': (0.6, -2.2),
    'southwest': (0.6, -2.2),
    'south': (-2.2, -5.0),
    'shaded': (-3.3, -6.1),  # Fully shaded treated as interior partitions/fully shaded
    'roof': (15.0, 12.2),
    'floor_over_uncond': (-3.3, -6.1)
}

# Duct gain multipliers (Table 1) simplified mapping by (location, rsi_category)
DUCT_GAIN_MULTIPLIERS = {
    ('attic_open', 'none'): 0.25,
    ('attic_open', '0.70'): 0.15,
    ('attic_open', '1.4'): 0.10,
    ('attic_open', '2.1'): 0.05,
    ('attic_open', '>=3.5'): 0.0,
    ('unconditioned_bsmt', '2.1'): 0.05,
    ('crawl_enclosed', 'none'): 0.25,
    ('crawl_enclosed', '0.7'): 0.15,
    ('crawl_enclosed', '1.4'): 0.10,
    ('crawl_enclosed', '3.5'): 0.05,
    ('crawl_enclosed', '>=3.5'): 0.0,
    ('slab_perimeter', 'any'): 0.0  # slab-on-grade perimeter gain handled differently; no cooling gain multiplier
}

# Latent multiplier default (Clause 6.3.1 Note)
DEFAULT_LATENT_MULTIPLIER = 1.3

# Minimum / maximum capacity rules (Clauses 6.3.2, 6.3.4, 6.3.5)
MIN_CAPACITY_REDUCTION_W = 1800  # allowed reduction from nominal
MAX_CAPACITY_PERCENT = 1.25
SMALL_LOAD_THRESHOLD_W = 6000
SMALL_LOAD_EXTRA_W = 1750

def latitude_factor(latitude):
    """Clause 6.2.2.1: LFactor = 1 + (Latitude - 40)*0.0375 bounded [0.7, 1.3].
    Applied ONLY to south and southeast/southwest windows.
    """
    lf = 1 + (latitude - 40.0) * 0.0375
    return max(0.7, min(1.3, lf))

def load_shading_factors(csv_path="f280_shading_factors.csv"):
    """Optional shading factors (Table 4) override.
    Expected columns: orientation, SFactor (fraction transmitted).
    If absent, default SFactor=1.0 except user-specified CLI global value.
    """
    p = Path(csv_path)
    if not p.exists():
        return {}
    try:
        df = pd.read_csv(p)
        out = {}
        for _, row in df.iterrows():
            ori = str(row.get('orientation','')).lower().strip()
            sf = row.get('SFactor')
            try:
                out[ori] = float(sf)
            except (TypeError, ValueError):
                continue
        return out
    except Exception:
        return {}

def determine_solar0(direction):
    """Return Solar0 base value for a given simplified direction label."""
    return SOLAR0_BASE.get(direction, SOLAR0_BASE['south'])

def duct_gain_multiplier(location, rsi):
    """Map duct location + insulation RSI to gain multiplier (DGMc)."""
    # Categorize RSI
    if rsi is None:
        cat = 'none'
    elif rsi >= 3.5:
        cat = '>=3.5'
    elif rsi >= 2.1:
        cat = '2.1'
    elif rsi >= 1.4:
        cat = '1.4'
    elif rsi >= 0.7:
        cat = '0.70' if location == 'attic_open' else '0.7'
    else:
        cat = 'none'
    key = (location, cat)
    return DUCT_GAIN_MULTIPLIERS.get(key, 0.0)

# Internal gains (typical residential)
OCCUPANT_SENSIBLE_GAIN = 75  # W per person (sensible)
OCCUPANT_LATENT_GAIN = 55  # W per person (latent)
APPLIANCE_GAIN = 300  # W (average appliance load)
LIGHTING_GAIN = 150  # W (average lighting load)

# Air properties
AIR_DENSITY = 1.2  # kg/m³
AIR_SPECIFIC_HEAT = 1006  # J/(kg·K)


def calculate_cooling_load_f280(params, inputs):
    """
    Calculate cooling load using CSA F280-12 method.
    
    CSA F280-12 Sections 6.1-6.3:
    - Envelope heat gains (transmission)
    - Solar heat gains (glazing)
    - Internal gains (occupants, appliances, lighting)
    - Ventilation and infiltration loads
    - Duct gains (if applicable)
    - Latent loads handled via latent_multiplier applied to sensible total
    """
    
    # Extract building geometry
    footprint = params.get('footprint', 100.0)  # m²
    storeys = params.get('storeys', 2)
    aspect = params.get('aspect', 1.0)
    fdwr = params.get('fdwr', 0.2)  # Fraction (not percentage)
    
    # Calculate dimensions
    width = math.sqrt(footprint / aspect)
    length = aspect * width
    floor_height = 2.5  # m (standard ceiling height)
    
    # Calculate areas - use parsed values if available from IDF, otherwise estimate from footprint
    floor_area_total = footprint * storeys
    perimeter = 2 * (length + width)
    
    # Use parsed geometry if available (more accurate than estimation)
    if params.get('idf_parsed'):
        wall_area_gross = sum(params.get('wall_orientation_areas', {}).values())
        window_area = params.get('parsed_window_area_m2', 0.0)
        wall_area_opaque = wall_area_gross - window_area  # For heat loss calculations
        roof_area = params.get('footprint', footprint)
    else:
        # Fallback: estimate from footprint
        wall_area_gross = perimeter * floor_height * storeys
        window_area = wall_area_gross * fdwr
        wall_area_opaque = wall_area_gross - window_area  # For heat loss calculations
        roof_area = footprint
    
    # Extract envelope properties
    ag_walls_rsi = params.get('ag_walls_rsi', 3.5)  # m²·K/W
    windows_u = params.get('windows_u', 1.8)  # W/(m²·K)
    attic_rsi = params.get('attic_rsi', 8.8)  # m²·K/W
    shgc = params.get('shgc', 0.4)  # Solar Heat Gain Coefficient
    ach = params.get('ach', 2.5)  # Air changes per hour at 50 Pa
    
    # Convert RSI to U-value (W/(m²·K))
    wall_u = 1.0 / ag_walls_rsi if ag_walls_rsi > 0 else 0.5
    roof_u = 1.0 / attic_rsi if attic_rsi > 0 else 0.2
    
    # Design temperature difference
    # Clause 6.1.1 design temperature difference
    Toc = inputs.get('outdoor_cooling_temp', OUTDOOR_TEMP_DESIGN)
    Tic = inputs.get('indoor_cooling_temp', INDOOR_TEMP_DESIGN)
    delta_t = Toc - Tic
    daily_temp_range = inputs.get('daily_temp_range', 14.0)  # For SC selection
    dtr_class_index = 0 if daily_temp_range <= 14.0 else 1
    
    # ========== SENSIBLE COOLING LOADS ==========
    
    # 1. Opaque assemblies conductive gain with solar correction (HGcop)
    orientations = ['north','east','south','west']
    wall_orientation_areas = params.get('wall_orientation_areas')
    hgcop_components = []
    if wall_orientation_areas:
        for ori in orientations:
            area_o = wall_orientation_areas.get(ori, 0.0)
            if area_o <= 0:
                continue
            sc = SOLAR_CORRECTION_TABLE[ori][dtr_class_index]
            hgcop_components.append((area_o / ag_walls_rsi) * (delta_t + sc))
    else:
        wall_area_each = wall_area_opaque / 4.0
        for ori in orientations:
            sc = SOLAR_CORRECTION_TABLE[ori][dtr_class_index]
            hgcop_components.append((wall_area_each / ag_walls_rsi) * (delta_t + sc))
    # Roof solar correction
    sc_roof = SOLAR_CORRECTION_TABLE['roof'][dtr_class_index]
    hgcop_roof = (roof_area / attic_rsi) * (delta_t + sc_roof)
    # Windows conduction portion handled separately (DTDc/RSI term in HGct)
    q_windows_conduction = window_area * windows_u * delta_t  # For reporting only (not SC adjusted)
    HGcop_total = sum(hgcop_components) + hgcop_roof
    
    # 2. Transparent assemblies heat gain HGct = A * (SHGC*Solar + DTDc/RSI)
    latitude = inputs.get('latitude_deg', 49.0)
    lfactor = latitude_factor(latitude)
    sfactor_global = inputs.get('shading_factor_global', 1.0)
    shading_factors = load_shading_factors()
    window_orientation_areas = params.get('window_orientation_areas')
    window_area_each = window_area / 4.0 if not window_orientation_areas else None
    window_rsi = 1.0 / windows_u if windows_u > 0 else 0.5
    solar_details = {}
    hgct_components = []
    for ori in orientations:
        area_o = window_orientation_areas.get(ori, 0.0) if window_orientation_areas else window_area_each
        if area_o is None or area_o <= 0:
            continue
        if ori == 'south':
            solar0 = determine_solar0('south') * lfactor
        elif ori in ['east','west']:
            solar0 = determine_solar0('east')
        elif ori == 'north':
            solar0 = determine_solar0('north')
        else:
            solar0 = determine_solar0(ori)
        sfactor = shading_factors.get(ori, sfactor_global)
        solar_incident = solar0 * sfactor
        hgct = area_o * (shgc * solar_incident + delta_t / window_rsi)
        hgct_components.append(hgct)
        solar_details[ori] = {'Solar0': solar0,'SFactor': sfactor,'Solar_incident': solar_incident,'Area_m2': area_o}
    HGct_total = sum(hgct_components)
    q_solar_total = sum(d['Solar_incident'] * d['Area_m2'] * shgc for d in solar_details.values())
    
    # 3. Internal gains (CSA F280 Section 6.1.3)
    # Assume occupancy based on floor area (1 person per 25 m²)
    bedrooms = inputs.get('bedrooms')
    if params.get('occupant_count'):
        num_occupants = params['occupant_count']
    elif bedrooms is not None:
        num_occupants = bedrooms + 1
    else:
        num_occupants = max(2, int(round(floor_area_total / 40.0)))
    # Clause 6.2.4 occupants sensible 70 W
    HGsp = num_occupants * 70.0
    # Clause 6.2.5 appliances & plug loads: 4 W/m2, not lower than 800 W
    HGapk_raw = 4.0 * floor_area_total
    HGapk = max(800.0, HGapk_raw)
    # Laundry equipment (optional) - allow override else assume 0
    HGapl = inputs.get('laundry_gain_W', 0.0)
    q_internal = HGsp + HGapk + HGapl
    
    # 4. Infiltration load (CSA F280 Section 6.1.4)
    # Use AIM-2 model to calculate natural ACH if available and configured
    volume = floor_area_total * floor_height
    # Clause 6.2.6: HGsalb = LFair * Vb/3.6 * DTDc * 1.2  (LFair in air changes per hour)
    lfair = inputs.get('natural_ach')
    
    # Try to use AIM-2 if parameters provided and module available
    if lfair is None and aim2_compute_infiltration is not None:
        aim2_params = {
            'n': inputs.get('aim2_flow_exponent', 0.67),
            'n50': inputs.get('n50', ach),
            'house_type': inputs.get('house_type', 'detached'),
            'foundation': inputs.get('foundation', 'crawl'),
            'storeys': str(storeys) if storeys in [1, 2, 3] else '2',
            'flue_diam_mm': inputs.get('flue_diam_mm', 0.0),
            'shelter_walls': inputs.get('shelter_walls', 1.0),
            'shelter_flue': inputs.get('shelter_flue', 1.0),
            'terrain_class_met': inputs.get('terrain_class_met', 3),
            'terrain_class_site': inputs.get('terrain_class_site', 7),
        }
        
        # Get wind speed from climate data or inputs
        wind_speed_met = inputs.get('wind_speed_met')
        if wind_speed_met is None:
            city_name = inputs.get('city')
            if city_name:
                # Load from NBC_F280_merge.csv
                climate_db = _load_climate_db()
                if climate_db is not None:
                    city_row = climate_db[climate_db['City'].str.lower() == city_name.lower()]
                    if not city_row.empty:
                        # Use JulWind column for cooling season
                        wind_speed_met = city_row.iloc[0].get('JulWind', 10.0)
            
            if wind_speed_met is None:
                wind_speed_met = 10.0  # Default 10 km/h
        
        try:
            aim2_inputs = AIM2Inputs(
                volume=volume,
                n=aim2_params['n'],
                n50=aim2_params['n50'],
                eave_height=storeys * floor_height,
                indoor_temp=Tic,
                outdoor_temp=Toc,
                temp_unit='C',
                flue_diam_mm=aim2_params['flue_diam_mm'],
                wind_speed_met=wind_speed_met,
                met_height=10.0,
                terrain_class_met=aim2_params['terrain_class_met'],
                terrain_class_site=aim2_params['terrain_class_site'],
                shelter_walls=aim2_params['shelter_walls'],
                shelter_flue=aim2_params['shelter_flue'],
                house_type=aim2_params['house_type'],
                foundation=aim2_params['foundation'],
                storeys=aim2_params['storeys']
            )
            aim2_results = aim2_compute_infiltration(aim2_inputs)
            lfair = aim2_results['Qnat']  # Natural ACH from AIM-2
        except Exception:
            # Fall back to simple approximation if AIM-2 fails
            lfair = None
    
    if lfair is None:
        # Fallback: approximate from ACH50 if provided
        ach50 = ach
        lfair = ach50 / inputs.get('ach50_to_natural_divisor', 20.0)
    
    HGsalb = lfair * volume / 3.6 * delta_t * 1.2
    q_infiltration_sensible = HGsalb  # building-level
    infiltration_m3s = (lfair * volume) / 3600.0  # Convert ACH to m³/s
    
    # 5. Ventilation load (CSA F280 Section 6.1.5)
    # Assume 0.3 ACH for mechanical ventilation
    # Clause 6.2.7 ventilation
    PVC_ls = inputs.get('PVC_l_s')  # principal ventilation rate L/s
    if PVC_ls is None:
        # fallback simple estimate 0.3 ACH in L/s
        ventilation_m3s = (0.3 * volume) / 3600.0
        PVC_ls = ventilation_m3s * 1000.0
    hrv_efficiency = params.get('hrv_efficiency', inputs.get('hrv_efficiency', 0.0))
    HGsvb = PVC_ls * delta_t * 1.2 * (1 - hrv_efficiency)
    q_ventilation_sensible = HGsvb
    
    # Total sensible cooling load
    # Conductive total HGcr = opaque (HGcop_total) + transparent (HGct_total)
    HGcr = HGcop_total + HGct_total
    # Distribute infiltration & ventilation proportionally to conductive (room-based). Single-room approximation.
    HGsalr = HGsalb  # whole building single aggregation
    HGsvr = HGsvb
    # Duct gain Clause 6.2.8
    dgm_location = inputs.get('duct_location', 'attic_open')
    duct_rsi = inputs.get('duct_rsi')
    DGMc = duct_gain_multiplier(dgm_location, duct_rsi)
    HGdr = DGMc * (HGcr + HGsp + HGapk + HGapl + HGsalr + HGsvr)
    # Clause 6.2.9 total sensible
    HGsr = HGcr + HGsp + HGapk + HGapl + HGsalr + HGdr + HGsvr
    q_sensible_total = HGsr
    
    # ========== EQUIPMENT SIZING (CSA F280 Section 6.3) ==========
    
    # Apply latent multiplier per Clause 6.3.1
    # This accounts for latent loads (occupants, infiltration, ventilation moisture)
    latent_multiplier = inputs.get('latent_multiplier', DEFAULT_LATENT_MULTIPLIER)
    CSCn_nominal = latent_multiplier * HGsr  # Clause 6.3.1
    # Allowable installed capacity range
    min_capacity = max(CSCn_nominal - MIN_CAPACITY_REDUCTION_W, 0.8 * CSCn_nominal)
    max_capacity = MAX_CAPACITY_PERCENT * CSCn_nominal
    if CSCn_nominal < SMALL_LOAD_THRESHOLD_W:
        max_capacity = CSCn_nominal + SMALL_LOAD_EXTRA_W
    # Convert to tons (nominal)
    cooling_tons_nominal = CSCn_nominal / W_PER_TON
    cooling_tons_min = min_capacity / W_PER_TON
    cooling_tons_max = max_capacity / W_PER_TON
    cooling_tons_selected = 0.25 * math.ceil(cooling_tons_nominal / 0.25)
    cooling_w_selected = cooling_tons_selected * W_PER_TON
    
    # Calculate required airflow (400 CFM per ton)
    airflow_cfm = CFM_PER_TON * cooling_tons_selected
    airflow_m3s = airflow_cfm * CFM_TO_M3S
    
    # Return detailed results
    return {
        # Load components (W)
        'HGcop_total_W': HGcop_total,
        'HGct_total_W': HGct_total,
        'HGcr_conductive_total_W': HGcr,
        'HGsp_people_W': HGsp,
        'HGapk_appliances_W': HGapk,
        'HGapl_laundry_W': HGapl,
        'HGsalb_infiltration_W': HGsalb,
        'HGsvb_ventilation_W': HGsvb,
        'HGdr_duct_gain_W': HGdr,
        'HGsr_sensible_total_W': HGsr,
        'q_solar_transmission_W': q_solar_total,
        'q_windows_conduction_W': q_windows_conduction,
        'CSCn_nominal_W': CSCn_nominal,
        'CSCn_min_allowable_W': min_capacity,
        'CSCn_max_allowable_W': max_capacity,
        'cooling_tons_nominal': cooling_tons_nominal,
        'cooling_tons_min': cooling_tons_min,
        'cooling_tons_max': cooling_tons_max,
        'cooling_tons_selected': cooling_tons_selected,
        'cooling_w_selected': cooling_w_selected,
        
        # Equipment sizing
        'cooling_kw_selected': cooling_w_selected / 1000.0,
        'cooling_kbtu_h_selected': cooling_w_selected / BTU_TO_W / 1000.0,
        'airflow_cfm': airflow_cfm,
        'airflow_m3s': airflow_m3s,
        'airflow_ls': airflow_m3s * 1000.0,
        
        # Building geometry (for reference)
        'floor_area_total_m2': floor_area_total,
        'wall_area_gross_m2': wall_area_gross,
        'wall_area_opaque_m2': wall_area_opaque,
        'window_area_m2': window_area,
        'roof_area_m2': roof_area,
        'volume_m3': volume,
        'num_occupants': num_occupants,
        'latitude_deg': latitude,
        'LFactor': lfactor,
        'daily_temp_range_C': daily_temp_range,
        'latent_multiplier_used': latent_multiplier,
        'duct_gain_multiplier': DGMc,
        'wall_orientation_areas_m2': wall_orientation_areas if wall_orientation_areas else {},
        'window_orientation_areas_m2': window_orientation_areas if window_orientation_areas else {},
        'solar_details': solar_details,
        'natural_ach_used': lfair,
        'infiltration_method': 'AIM-2' if (inputs.get('natural_ach') is None and aim2_compute_infiltration is not None) else 'Simple',
    }


def parse_idf_parameters(sim_dir: Path):
    """Parse IDF for simplified geometry, glazing, and occupant estimation."""
    try:
        idf_files = list(sim_dir.glob('*.idf'))
        if not idf_files:
            return {}
        idf_path = idf_files[0]
        text = idf_path.read_text(encoding='utf-8', errors='ignore')
        raw_objects = text.split(';')
        wall_area_gross = 0.0
        window_area = 0.0
        roof_area = 0.0
        floor_areas = []
        wall_orientation_areas = {'north':0.0,'east':0.0,'south':0.0,'west':0.0}
        window_orientation_areas = {'north':0.0,'east':0.0,'south':0.0,'west':0.0}
        glazing_u = None
        glazing_shgc = None
        occupant_count = 0.0
        zone_floor_area = {}
        wall_orientation_by_surface = {}

        def clean_token(token):
            token = token.strip()
            if '!' in token:
                token = token.split('!')[0].strip()
            return token

        def polygon_area(vertices):
            if len(vertices) < 3:
                return 0.0
            v0 = vertices[0]
            translated = [(x - v0[0], y - v0[1], z - v0[2]) for (x,y,z) in vertices]
            cx = cy = cz = 0.0
            for i in range(len(translated)):
                x1,y1,z1 = translated[i]
                x2,y2,z2 = translated[(i+1) % len(translated)]
                cx += y1*z2 - z1*y2
                cy += z1*x2 - x1*z2
                cz += x1*y2 - y1*x2
            return 0.5 * math.sqrt(cx*cx + cy*cy + cz*cz)

        for obj in raw_objects:
            lines = [l for l in obj.split('\n') if l.strip()]
            if not lines:
                continue
            first = lines[0].strip().lower()
            if first.startswith('buildingsurface:detailed'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                if len(tokens) < 10:
                    continue
                surface_name = tokens[0]
                surface_type = tokens[1].lower()
                construction_name = tokens[2] if len(tokens) > 2 else ''
                zone_name = tokens[3] if len(tokens) > 3 else ''
                outside_bc = tokens[4].lower() if len(tokens) > 4 else ''
                num_vertices = None
                coords_start_index = None
                
                # Try to find explicit vertex count in the tokens
                for idx, tk in enumerate(tokens):
                    if num_vertices is None and re.fullmatch(r'[0-9]+', tk):
                        num_vertices = int(tk)
                        coords_start_index = idx + 1
                        break
                
                # If no explicit vertex count found, look for coordinate data
                # Starting from a reasonable position (after standard fields)
                if not num_vertices:
                    # Standard fields: Name, Type, Construction, Zone, Space, OutsideBoundary, 
                    # BoundaryObject, SunExposure, WindExposure, ViewFactor, NumVertices (often blank)
                    # So coordinates typically start around index 10-11
                    search_start = min(10, len(tokens))
                    coord_tokens = []
                    for idx in range(search_start, len(tokens)):
                        try:
                            float(tokens[idx])
                            coord_tokens.append(tokens[idx])
                        except ValueError:
                            # Stop when we hit a non-numeric token
                            break
                    
                    # Calculate number of vertices from coordinate count
                    if len(coord_tokens) >= 9 and len(coord_tokens) % 3 == 0:
                        num_vertices = len(coord_tokens) // 3
                        coords_start_index = search_start
                
                if not num_vertices or coords_start_index is None:
                    continue
                    
                coord_tokens = tokens[coords_start_index:coords_start_index + num_vertices*3]
                if len(coord_tokens) < num_vertices*3:
                    continue
                verts = []
                try:
                    for i in range(0, num_vertices*3, 3):
                        verts.append((float(coord_tokens[i]), float(coord_tokens[i+1]), float(coord_tokens[i+2])))
                except ValueError:
                    continue
                area = polygon_area(verts)
                if outside_bc == 'outdoors':
                    if surface_type == 'wall':
                        # Exclude basement and attic walls from above-grade wall calculation
                        if 'basement' not in zone_name.lower() and 'attic' not in zone_name.lower():
                            wall_area_gross += area
                            if len(verts) >= 3:
                                v1,v2,v3 = verts[0],verts[1],verts[2]
                                nx = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v2[2]-v1[2])*(v3[1]-v1[1])
                                ny = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v2[0]-v1[0])*(v3[2]-v1[2])
                                nz = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v2[1]-v1[1])*(v3[0]-v1[0])
                                if abs(nz) < 0.2:
                                    az = (math.degrees(math.atan2(nx, ny)) + 360.0) % 360.0
                                    if az >= 315 or az < 45:
                                        wall_orientation_areas['north'] += area
                                        wall_orientation_by_surface[surface_name] = 'north'
                                    elif az < 135:
                                        wall_orientation_areas['east'] += area
                                        wall_orientation_by_surface[surface_name] = 'east'
                                    elif az < 225:
                                        wall_orientation_areas['south'] += area
                                        wall_orientation_by_surface[surface_name] = 'south'
                                    else:
                                        wall_orientation_areas['west'] += area
                                        wall_orientation_by_surface[surface_name] = 'west'
                    elif surface_type in ['roof','roofceiling','ceiling']:
                        roof_area += area
                # Collect floor areas (check all floors, not just outdoors)
                if surface_type == 'floor':
                    # Only count above-grade conditioned floors (exclude basement and attic)
                    if 'basement' not in zone_name.lower() and 'attic' not in zone_name.lower():
                        floor_areas.append(area)
            elif first.startswith('fenestrationsurface:detailed'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                if len(tokens) < 10:
                    continue
                surf_type = tokens[1].lower()
                num_vertices = None
                coords_start_index = None
                
                # Try to find explicit vertex count
                for idx, tk in enumerate(tokens):
                    if num_vertices is None and re.fullmatch(r'[0-9]+', tk):
                        num_vertices = int(tk)
                        coords_start_index = idx + 1
                        break
                
                # If no explicit vertex count, try to count coordinate triplets
                if not num_vertices:
                    # FenestrationSurface fields: Name, Type, Construction, Surface, OutsideBoundary,
                    # ViewFactor, FrameDivider, Multiplier, NumVertices (often blank)
                    search_start = min(8, len(tokens))
                    coord_tokens = []
                    for idx in range(search_start, len(tokens)):
                        try:
                            float(tokens[idx])
                            coord_tokens.append(tokens[idx])
                        except ValueError:
                            break
                    
                    if len(coord_tokens) >= 9 and len(coord_tokens) % 3 == 0:
                        num_vertices = len(coord_tokens) // 3
                        coords_start_index = search_start
                
                if not num_vertices or coords_start_index is None:
                    continue
                    
                coord_tokens = tokens[coords_start_index:coords_start_index + num_vertices*3]
                if len(coord_tokens) < num_vertices*3:
                    continue
                verts = []
                try:
                    for i in range(0, num_vertices*3, 3):
                        verts.append((float(coord_tokens[i]), float(coord_tokens[i+1]), float(coord_tokens[i+2])))
                except ValueError:
                    continue
                area = polygon_area(verts)
                if surf_type in ['window','glassdoor','skylight']:
                    window_area += area
                    if len(verts) >= 3:
                        v1,v2,v3 = verts[0],verts[1],verts[2]
                        nx = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v2[2]-v1[2])*(v3[1]-v1[1])
                        ny = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v2[0]-v1[0])*(v3[2]-v1[2])
                        nz = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v2[1]-v1[1])*(v3[0]-v1[0])
                        if abs(nz) < 0.2:
                            az = (math.degrees(math.atan2(nx, ny)) + 360.0) % 360.0
                            if az >= 315 or az < 45:
                                window_orientation_areas['north'] += area
                            elif az < 135:
                                window_orientation_areas['east'] += area
                            elif az < 225:
                                window_orientation_areas['south'] += area
                            else:
                                window_orientation_areas['west'] += area
            elif first.startswith('window'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                # Format: Name, Construction, ParentSurface, StartX, StartZ, LengthX, HeightZ
                # (blank fields for Frame/Divider and Multiplier are filtered out by clean_token)
                if len(tokens) < 7:
                    continue
                parent_surface = tokens[2]
                try:
                    start_x = float(tokens[3]); start_z = float(tokens[4])  # not used, retained for potential future placement logic
                    length_x = float(tokens[5]); height_z = float(tokens[6])
                except ValueError:
                    continue
                area = length_x * height_z
                if area <= 0:
                    continue
                window_area += area
                orientation = wall_orientation_by_surface.get(parent_surface)
                if orientation in window_orientation_areas:
                    window_orientation_areas[orientation] += area
            elif first.startswith('windowmaterial:simpleglazingsystem'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                if len(tokens) >= 4:
                    try:
                        uval = float(tokens[1]); shgc = float(tokens[2])
                        glazing_u = uval if glazing_u is None else min(glazing_u, uval)
                        glazing_shgc = shgc if glazing_shgc is None else max(glazing_shgc, shgc)
                    except ValueError:
                        pass
            elif first.startswith('zone'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                if len(tokens) >= 10:
                    try:
                        floor_a = float(tokens[9])
                        if floor_a > 0:
                            zone_floor_area[tokens[0]] = floor_a
                    except Exception:
                        pass
            elif first.startswith('people'):
                tokens = [clean_token(t) for t in ','.join(lines[1:]).split(',') if clean_token(t)]
                if len(tokens) >= 7:
                    zone_name = tokens[1]
                    method = tokens[3].lower()
                    if method == 'area/person':
                        try:
                            area_per_person = float(tokens[6])
                            zone_area = zone_floor_area.get(zone_name)
                            if zone_area and area_per_person > 0:
                                occupant_count += zone_area / area_per_person
                        except Exception:
                            pass

        footprint = max(floor_areas) if floor_areas else (roof_area if roof_area > 0 else None)
        floor_area_total = sum(floor_areas) if floor_areas else (footprint if footprint else None)
        storeys = None
        if footprint and floor_area_total:
            est = floor_area_total / footprint
            if est >= 0.9:
                storeys = int(round(est))
        params = {}
        if footprint:
            params['footprint'] = footprint
        if storeys:
            params['storeys'] = storeys
        if wall_area_gross > 0:
            params['aspect'] = 1.0
            # Store parsed wall area for direct use
            params['parsed_wall_area_gross_m2'] = wall_area_gross
        if glazing_u is not None:
            params['windows_u'] = glazing_u
        if glazing_shgc is not None:
            params['shgc'] = glazing_shgc
        if roof_area > 0 and 'footprint' not in params:
            params['footprint'] = roof_area
        if occupant_count > 0:
            params['occupant_count'] = int(round(occupant_count))
        params['wall_orientation_areas'] = wall_orientation_areas
        params['window_orientation_areas'] = window_orientation_areas
        params['parsed_window_area_m2'] = window_area
        params['idf_parsed'] = True
        return params
    except Exception:
        return {}

def f280_cooling_load_from_idf(idf_path: str, **overrides):
    """Compute CSA F280 cooling load from a single IDF file.
    idf_path: path to an EnergyPlus IDF file.
    overrides: optional keyword parameters supplying data not present in IDF
        ag_walls_rsi, attic_rsi, windows_u, shgc, ach,
        outdoor_cooling_temp, indoor_cooling_temp, daily_temp_range,
        latitude_deg, longitude_deg, natural_ach, PVC_l_s, bedrooms,
        duct_location, duct_rsi, shading_factor_global, latent_multiplier.
    Returns: dict of load components and sizing values.
    """
    idf_path_obj = Path(idf_path)
    sim_dir = idf_path_obj.parent
    params = parse_idf_parameters(sim_dir)
    # Allow direct overrides of parsed envelope values
    for k in ['ag_walls_rsi','attic_rsi','windows_u','shgc','ach']:
        if k in overrides and overrides[k] is not None:
            params[k] = overrides[k]
    # Integrate nearest climate data if latitude & longitude provided
    lat = overrides.get('latitude_deg')
    lon = overrides.get('longitude_deg')
    climate_meta = {}
    if lat is not None and lon is not None:
        try:
            climate_row = _nearest_climate_row(lat, lon)
            if climate_row is not None:
                # Populate outdoor cooling design temp if not explicitly overridden
                if 'outdoor_cooling_temp' not in overrides and 'DCDBT' in climate_row:
                    try:
                        overrides['outdoor_cooling_temp'] = float(climate_row['DCDBT'])
                    except (TypeError, ValueError):
                        pass
                # Populate daily temperature range (summer) if not explicitly overridden
                if 'daily_temp_range' not in overrides and 'Strange' in climate_row:
                    try:
                        overrides['daily_temp_range'] = float(climate_row['Strange'])
                    except (TypeError, ValueError):
                        pass
                # Collect metadata
                climate_meta.update({
                    'climate_site_name': climate_row.get('Location') or climate_row.get('City'),
                    'climate_site_city': climate_row.get('City'),
                    'climate_site_province': climate_row.get('Prov'),
                    'climate_site_lat': _safe_float(climate_row.get('Latitude')),
                    'climate_site_lon': _safe_float(climate_row.get('Longitude')),
                    'climate_site_dcdbt_C': _safe_float(climate_row.get('DCDBT')),
                    'climate_site_dhdbt_C': _safe_float(climate_row.get('DHDBT')),
                    'climate_site_jul_2p5_db_C': _safe_float(climate_row.get('July 2.5% DB (°C)')),
                    'climate_site_jul_2p5_wb_C': _safe_float(climate_row.get('July 2.5% WB (°C)')),
                    'climate_site_daily_temp_range_C': _safe_float(climate_row.get('Strange')),
                    'climate_zone_code': climate_row.get('Climate_Zone'),
                })
        except Exception:
            pass
    results = calculate_cooling_load_f280(params, overrides)
    if climate_meta:
        results.update(climate_meta)
    return results

# ------------------ Climate Data Helpers ------------------
_CLIMATE_DB = None

def _safe_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None

def _load_climate_db(path='NBC_F280_merge.csv'):
    global _CLIMATE_DB
    if _CLIMATE_DB is not None:
        return _CLIMATE_DB
    csv_path = Path(path)
    if not csv_path.exists():
        _CLIMATE_DB = None
        return None
    try:
        # First line is banner, second line headers
        # Use engine='c' for faster parsing, low_memory=False to avoid dtype warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df = pd.read_csv(csv_path, header=1, engine='c', low_memory=False)
        # Convert lat/lon to numeric upfront and filter (only keep needed columns for speed)
        df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
        df['Longitude'] = pd.to_numeric(df['Longitude'], errors='coerce')
        df = df[df['Latitude'].notnull() & df['Longitude'].notnull()].copy()
        _CLIMATE_DB = df.reset_index(drop=True)
        return _CLIMATE_DB
    except Exception as e:
        warnings.warn(f"Could not load climate database: {e}")
        _CLIMATE_DB = None
        return None

def _nearest_climate_row(lat, lon, path='NBC_F280_merge.csv'):
    df = _load_climate_db(path)
    if df is None or df.empty:
        return None
    # Vectorized haversine (using numpy from module imports)
    lat_arr = df['Latitude'].values
    lon_arr = df['Longitude'].values
    lat1 = math.radians(lat)
    lon1 = math.radians(lon)
    lat2 = np.radians(lat_arr)
    lon2 = np.radians(lon_arr)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + math.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distances = 6371.0 * c
    idx = int(distances.argmin())
    return df.iloc[idx].to_dict()

    # Legacy multi-building batch function removed; use single-IDF interface instead.


def _run_cli():
    import argparse
    parser = argparse.ArgumentParser(description='CSA F280-12 cooling load from a single IDF',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Provide envelope overrides if not derivable from IDF surfaces.')
    parser.add_argument('--idf', required=True, help='Path to IDF file')
    parser.add_argument('-o', '--output', type=str, help='Optional CSV output path for single-row results')
    # Envelope overrides
    parser.add_argument('--ag-walls-rsi', type=float, help='Above-grade wall RSI (m2*K/W)')
    parser.add_argument('--attic-rsi', type=float, help='Attic/roof insulation RSI (m2*K/W)')
    parser.add_argument('--windows-u', type=float, help='Window U-value (W/m2*K) override')
    parser.add_argument('--shgc', type=float, help='Window SHGC override')
    # Design & site
    parser.add_argument('--latitude', type=float, help='Site latitude (deg) override')
    parser.add_argument('--longitude', type=float, help='Site longitude (deg) override (for nearest climate match)')
    parser.add_argument('--dtr', type=float, help='Daily temperature range (°C) for SC table selection')
    parser.add_argument('--tic', type=float, help='Indoor cooling design temperature (°C)')
    parser.add_argument('--toc', type=float, help='Outdoor cooling design temperature (°C)')
    parser.add_argument('--natural-ach', type=float, help='Natural air change rate (ACH) override')
    parser.add_argument('--ventilation-lps', type=float, help='Principal ventilation rate (L/s) override')
    parser.add_argument('--duct-location', type=str, help='Duct location (attic_open, crawl_enclosed, unconditioned_bsmt, slab_perimeter)')
    parser.add_argument('--duct-rsi', type=float, help='Duct effective insulation RSI value')
    parser.add_argument('--bedrooms', type=int, help='Number of bedrooms (for occupant count)')
    parser.add_argument('--latent-multiplier', type=float, help='Latent multiplier for CSCn (default 1.3)')
    parser.add_argument('--shading-factor', type=float, help='Global shading factor (0-1)')
    args = parser.parse_args()

    overrides = {}
    for arg_key, override_key in [
        ('ag_walls_rsi','ag_walls_rsi'),('attic_rsi','attic_rsi'),('windows_u','windows_u'),('shgc','shgc'),
        ('latitude','latitude_deg'),('longitude','longitude_deg'),('dtr','daily_temp_range'),('tic','indoor_cooling_temp'),('toc','outdoor_cooling_temp'),
        ('natural_ach','natural_ach'),('ventilation_lps','PVC_l_s'),('duct_location','duct_location'),('duct_rsi','duct_rsi'),
        ('bedrooms','bedrooms'),('latent_multiplier','latent_multiplier'),('shading_factor','shading_factor_global')]:
        val = getattr(args, arg_key.replace('-','_'))
        if val is not None:
            overrides[override_key] = val

    results = f280_cooling_load_from_idf(args.idf, **overrides)
    import json as _json
    print(_json.dumps(results, indent=2))
    if args.output:
        import pandas as _pd
        _pd.DataFrame([results]).to_csv(args.output, index=False)

if __name__ == "__main__":
    _run_cli()
