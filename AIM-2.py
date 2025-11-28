
#!/usr/bin/env python3
"""
AIM-2 Infiltration Model (Bradley HOT2000 Implementation)
=========================================================

Implements Bradley's HOT2000 version of AIM-2:
- Leakage distribution (R, X, Y) using L. Lew (1993) Table for fractions
  (ceiling, floor, walls) selected by house type, foundation, and storeys.
- Flue leakage coefficient (from diameter)
- Stack and wind effect calculations
- Shelter coefficient (Bradley version)
- Wind factor (no flue; Bradley eqn)
- Superposition of stack and wind flows
- ACH conversions

Key leakage parameters:
-----------------------
- n   : Flow exponent (dimensionless). Typical range: 0.65–0.75 for houses.
        Appears in Q = C * (ΔP)^n.

- N50 : Air change rate at 50 Pa (ACH @ 50 Pa).
        Used to compute C:
            Q50 = N50 * Volume / 3600
            C   = Q50 / (50)^n

- C   : Flow coefficient (m³/s/Pa^n). Computed internally from N50 and n.

Leakage distribution source:
----------------------------
L. Lew, "Evaluation of AIM-2", EMR, Apr 23/93.
Fortran table mapping (ceiling, floor, walls) fractions by:
  i = house type (1 detached, 2 semi-detached)
  j = foundation (1 crawl, 2 slab-on-grade, 3 shallow, 4 full)
  k = storeys (1=1, 2=1.5, 3=2, 4=2.5, 5=3)

"""

import math
import argparse
from dataclasses import dataclass
from typing import Tuple, Dict

# ------------------------
# Physical constants
# ------------------------
RHO_AIR = 1.204097  # kg/m³ at ~20°C
G = 9.80665         # m/s²
AP_REF = 4.0        # Pa reference pressure for flue coefficient
B_INTERACTION = -1.0 / 3.0  # Bradley empirical interaction term

# ------------------------
# Terrain roughness (Davenport)
# ------------------------
DAVENPORT_ZO: Dict[int, float] = {
    1: 0.0002,  # Open sea
    2: 0.005,   # Mud flats
    3: 0.03,    # Open flat terrain (grass)
    4: 0.10,    # Low crops
    5: 0.25,    # High crops
    6: 0.50,    # Parkland/bushes
    7: 1.00,    # Suburb/forest
    8: 2.00     # City centre (assumed)
}

# ------------------------
# Leakage distribution tables (L. Lew, 1993)
# Each tuple = (ceiling, floor, walls)
# ------------------------

# Detached (i = 1), foundations j=1..4, storeys k=1..5
DETACHED_ROWS = [
    # j=1 Crawl Space
    [(0.20, 0.50, 0.30), (0.20, 0.50, 0.30), (0.15, 0.60, 0.25),
     (0.15, 0.60, 0.25), (0.10, 0.70, 0.20)],
    # j=2 Slab-on-grade
    [(0.30, 0.50, 0.20), (0.30, 0.50, 0.20), (0.20, 0.60, 0.20),
     (0.20, 0.65, 0.15), (0.15, 0.70, 0.15)],
    # j=3 Shallow
    [(0.30, 0.50, 0.20), (0.30, 0.50, 0.20), (0.20, 0.60, 0.20),
     (0.15, 0.70, 0.15), (0.10, 0.80, 0.10)],
    # j=4 Full
    [(0.30, 0.50, 0.20), (0.20, 0.60, 0.20), (0.20, 0.65, 0.15),
     (0.20, 0.70, 0.10), (0.10, 0.80, 0.10)]
]

# Semi-detached (i = 2), foundations j=1..4, storeys k=1..5
SEMI_ROWS = [
    # j=1 Crawl Space
    [(0.30, 0.40, 0.30), (0.30, 0.40, 0.30), (0.20, 0.50, 0.30),
     (0.20, 0.60, 0.20), (0.20, 0.60, 0.20)],
    # j=2 Slab-on-grade
    [(0.30, 0.40, 0.30), (0.30, 0.50, 0.20), (0.20, 0.60, 0.20),
     (0.20, 0.60, 0.20), (0.20, 0.65, 0.15)],
    # j=3 Shallow
    [(0.30, 0.50, 0.20), (0.30, 0.55, 0.15), (0.25, 0.60, 0.15),
     (0.20, 0.70, 0.10), (0.20, 0.70, 0.10)],
    # j=4 Full
    [(0.30, 0.50, 0.20), (0.30, 0.50, 0.20), (0.20, 0.60, 0.20),
     (0.20, 0.70, 0.10), (0.20, 0.70, 0.10)]
]

HOUSE_TYPE_MAP = {
    "detached": 1,
    "semi-detached": 2,
    "semi": 2
}

FOUNDATION_MAP = {
    "crawl": 1,
    "slab": 2,        # slab-on-grade
    "shallow": 3,
    "full": 4
}

STOREYS_MAP = {
    "1": 1, "1.0": 1,
    "1.5": 2,
    "2": 3, "2.0": 3,
    "2.5": 4,
    "3": 5, "3.0": 5
}

@dataclass
class Inputs:
    # Building & leakage
    volume: float
    n: float
    n50: float
    # Geometry & temps
    eave_height: float
    indoor_temp: float = 22.0  # Celsius (converted to Kelvin internally)
    outdoor_temp: float = 0.0  # Celsius (converted to Kelvin internally)
    temp_unit: str = 'C'  # 'C' for Celsius or 'K' for Kelvin
    # Flue
    flue_diam_mm: float = 0.0
    flue_height: float = None
    # Wind & terrain
    wind_speed: float = None
    wind_speed_met: float = 0.0
    met_height: float = 10.0
    terrain_class_met: int = 5
    terrain_class_site: int = 7
    shelter_walls: float = 1.0
    shelter_flue: float = 1.0
    # Leakage distribution selection
    house_type: str = "detached"   # 'detached' or 'semi-detached'
    foundation: str = "crawl"      # 'crawl'|'slab'|'shallow'|'full'
    storeys: str = "1"             # '1','1.5','2','2.5','3'

# ------------------------
# Helper functions
# ------------------------

def convert_to_kelvin(temp: float, unit: str = 'C') -> float:
    """Convert temperature to Kelvin.
    unit: 'C' for Celsius, 'K' for Kelvin
    """
    if unit.upper() == 'C':
        return temp + 273.15
    return temp

def C_from_n50(volume_m3: float, n: float, n50: float) -> float:
    """
    Compute flow coefficient C from N50 and n.
    N50 = ACH @ 50 Pa
    Q50 = N50 * Volume / 3600
    C   = Q50 / (50)^n
    """
    q50 = n50 * volume_m3 / 3600.0
    return q50 / (50.0 ** n)

def effective_wind_speed(inputs: Inputs) -> float:
    """
    Compute site windspeed:
    - If explicit wind_speed provided, use it.
    - Else adjust met tower wind to site using log-profile approximation (Davenport classes).
    """
    if inputs.wind_speed is not None:
        return inputs.wind_speed
    Zo_met = DAVENPORT_ZO.get(inputs.terrain_class_met, 0.25)
    Zo_site = DAVENPORT_ZO.get(inputs.terrain_class_site, 0.25)
    H = max(inputs.eave_height, 1.0)
    h_met = max(inputs.met_height, 1.0)
    num = math.log(H / Zo_site)*math.log(60/Zo_met)
    den = math.log(h_met / Zo_met)*math.log(60/Zo_site)
    if den <= 0.0:
        den = 1.0
    return inputs.wind_speed_met * (num / den) / 3.6  # Convert km/h to m/s

def select_fractions(house_type: str, foundation: str, storeys: str) -> Tuple[float, float, float]:
    """
    Select (ceiling, floor, walls) fractions from L. Lew table.

    house_type: 'detached' | 'semi-detached' (or 'semi')
    foundation: 'crawl' | 'slab' | 'shallow' | 'full'
    storeys: '1'|'1.5'|'2'|'2.5'|'3'
    """
    i = HOUSE_TYPE_MAP.get(house_type.lower())
    j = FOUNDATION_MAP.get(foundation.lower())
    k = STOREYS_MAP.get(storeys.strip().lower())

    if i is None:
        raise ValueError("Invalid house_type. Use 'detached' or 'semi-detached'.")
    if j is None:
        raise ValueError("Invalid foundation. Use 'crawl', 'slab', 'shallow', or 'full'.")
    if k is None:
        raise ValueError("Invalid storeys. Use one of '1','1.5','2','2.5','3'.")

    rows = DETACHED_ROWS if i == 1 else SEMI_ROWS
    frac = rows[j - 1][k - 1]  # tuples sum to 1.0
    return frac  # (a_c, a_f, a_w)

def leakage_coefficients(C_base: float, n: float, flue_diam_mm: float,
                         fractions: Tuple[float, float, float]) -> Tuple[float, float, float, float, float]:
    """
    Split base C across ceiling/floor/walls, add flue coefficient if flue present,
    and return total coefficients.
    """
    a_c, a_w, a_f = fractions
    # Normalize (safety) and distribute base C
    s = a_c + a_f + a_w
    if s <= 0.0:
        raise ValueError("Leakage fractions must sum > 0")
    a_c, a_f, a_w = a_c / s, a_f / s, a_w / s

    Cc0 = a_c * C_base
    Cf0 = a_f * C_base
    Cw0 = a_w * C_base

    # Flue coefficient
    Cflue = 0.0
    if flue_diam_mm and flue_diam_mm > 0.0:
        d_m = flue_diam_mm / 1000.0
        area = math.pi * (d_m ** 2) / 4.0
        # Bradley: C_flue = 0.5 * A_flue * (rho / (2 * ΔP_ref))^(n - 0.5)
        Cflue = 0.5 * area * ((RHO_AIR / (2.0 * AP_REF)) ** (n - 0.5))

    C_total = Cc0 + Cf0 + Cw0 + Cflue
    return C_total, Cc0, Cf0, Cw0, Cflue

def leakage_parameters(C_total: float, Cc0: float, Cf0: float, Cw0: float, Cflue: float, fractions: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Compute R, X, Y from coefficients:
      Y = C_flue / C_total
      R = (a_c + a_f) * (1 - Y)
      X = (a_c - a_f) * (1 - Y)
    where a_c = Cc0 / (C_total - C_flue), a_f similarly.
    """
    a_c, a_w, a_f = fractions
    Y = (Cflue / C_total) if C_total > 0.0 else 0.0
    R = (a_c + a_f) * (1.0 - Y)
    X = (a_c - a_f) * (1.0 - Y)
    return R, X, Y

def stack_pressure(eave_height: float, Ti: float, To: float) -> float:
    """
    Bradley: Ps = rho * g * H * (Ti - To) / Ti
    """
    return RHO_AIR * G * eave_height * (abs(Ti - To) / Ti)

def wind_pressure(Ue: float, Sw: float) -> float:
    """
    Bradley: Pw = 0.5 * rho * (Sw * Ue)^2
    """
    return 0.5 * RHO_AIR * (Sw * Ue) ** 2

def shelter_coefficient(Y: float, Swo: float, Swflue: float) -> float:
    """
    Bradley (HOT2000): Sw = Swo*(1 + Y) + Swflue*(1.5*Y)
    """
    return Swo * (1.0 - Y) + Swflue * (1.5 * Y) #This equation is shown as Swo * (1.0 + Y) + Swflue * (1.5 * Y) in Bradley but this is a typo

def fw_no_flue(n: float, R: float, X: float) -> float:
    """
    Bradley: fw = 0.19 * (2 - n) * (1 - (X + R)^(3/2))
    """
    return 0.19 * (2.0 - n) * (1.0 - ((X + R) / 2) ** 1.5)

def wind_flow(C_total: float, n: float, Ue: float, Sw: float, R: float, X: float, volume_m3: float) -> float:
    C_converted = C_total * 3600/volume_m3
    Pw = wind_pressure(Ue, Sw)
    fw = fw_no_flue(n, R, X)
    return C_converted * fw * (Pw ** n)

def stack_flow(C_total: float, n: float, Ps: float, fs: float, volume_m3: float) -> float:
    C_converted = C_total * 3600/volume_m3
    return C_converted * fs * (Ps ** n)

def fs_no_flue(n: float, R: float, X: float) -> float:
    """
    Bradley approximation for stack factor (no flue).
    """
    return ((1 + n * R) / (n + 1)) * ((0.5 - 0.5 * (((X ** 2)/(2 - R)) ** (5 / 4)))) ** (n + 1)

def superpose(Qs: float, Qw: float, n: float) -> float:
    """
    Bradley superposition:
      Q_nat = [ Qs^n + Qw^n + B * (Qs * Qw)^n ]^(1/n),  B ≈ -1/3
    """

    return ((Qs ** (1/n)) + (Qw ** (1/n)) + B_INTERACTION * ((Qs * Qw) ** (0.5 / n))) ** n

def compute_infiltration(inputs: Inputs) -> Dict[str, float]:
    # Convert temperatures to Kelvin if needed
    Ti = convert_to_kelvin(inputs.indoor_temp, inputs.temp_unit)
    To = convert_to_kelvin(inputs.outdoor_temp, inputs.temp_unit)
    
    # Derive C from N50 and n
    C_base = C_from_n50(inputs.volume, inputs.n, inputs.n50)

    # Fractions from L. Lew table
    fractions = select_fractions(inputs.house_type, inputs.foundation, inputs.storeys)

    # Coefficients
    C_total, Cc0, Cf0, Cw0, Cflue = leakage_coefficients(
        C_base=C_base,
        n=inputs.n,
        flue_diam_mm=inputs.flue_diam_mm,
        fractions=fractions
    )

    # Distribution parameters
    R, X, Y = leakage_parameters(C_total, Cc0, Cf0, Cw0, Cflue, fractions)

    # Stack flow
    Ps = stack_pressure(inputs.eave_height, Ti, To)
    fs = fs_no_flue(inputs.n, R, X)
    Qs = stack_flow(C_total, inputs.n, Ps, fs, inputs.volume)

    # Wind flow
    Ue = effective_wind_speed(inputs) #converted to m/s
    Sw = shelter_coefficient(Y, inputs.shelter_walls, inputs.shelter_flue)
    Qw = wind_flow(C_total, inputs.n, Ue, Sw, R, X, inputs.volume)

    # Superposition
    Qnat = superpose(Qs, Qw, inputs.n)

    return {
        "C_base": C_base, "C_total": C_total, "C_flue": Cflue,
        "Cc0": Cc0, "Cf0": Cf0, "Cw0": Cw0,
        "R": R, "X": X, "Y": Y,
        "Ps": Ps, "fs": fs,
        "Ue": Ue, "Sw": Sw,
        "Qs": Qs, "Qw": Qw,
        "Qnat": Qnat
    }

# ------------------------
# CLI
# ------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="AIM-2 infiltration calculator (Bradley HOT2000) with L. Lew fractions")
    # Required for C computation
    p.add_argument("--volume", type=float, required=True, help="House volume (m³)")
    p.add_argument("--n", type=float, default=0.67, help="Flow exponent (typical 0.65–0.75)")
    p.add_argument("--n50", type=float, required=True, help="Air change rate at 50 Pa (ACH @ 50 Pa)")
    # Geometry & temps
    p.add_argument("--eave_height", type=float, required=True, help="Eave height (m)")
    p.add_argument("--indoor_temp", type=float, default=24.0, help="Indoor temperature (°C)")
    p.add_argument("--outdoor_temp", type=float, default=28.0, help="Outdoor temperature (°C)")
    p.add_argument("--temp_unit", type=str, default='C', choices=['C', 'K'], help="Temperature unit (C or K)")
    # Flue
    p.add_argument("--flue_diam_mm", type=float, default=0.0, help="Flue diameter (mm)")
    p.add_argument("--flue_height", type=float, help="Flue height (m)")
    # Wind & terrain
    p.add_argument("--wind_speed", type=float, help="Explicit site windspeed (km/h). If omitted, model uses met tower + terrain correction.")
    p.add_argument("--wind_speed_met", type=float, default=0.0, help="Met tower windspeed (m/s)")
    p.add_argument("--met_height", type=float, default=10.0, help="Met tower height (m)")
    p.add_argument("--terrain_class_met", type=int, default=3, help="Davenport class @ met")
    p.add_argument("--terrain_class_site", type=int, default=7, help="Davenport class @ site")
    p.add_argument("--shelter_walls", type=float, default=1.0, help="Wall shelter coefficient Swo")
    p.add_argument("--shelter_flue", type=float, default=1.0, help="Flue shelter coefficient Swflue")
    # Leakage distribution selectors
    p.add_argument("--house_type", type=str, required=True, choices=["detached", "semi-detached", "semi"],
                   help="House type")
    p.add_argument("--foundation", type=str, required=True, choices=["crawl", "slab", "shallow", "full"],
                   help="Foundation type")
    p.add_argument("--storeys", type=str, required=True, choices=["1", "1.5", "2", "2.5", "3"],
                   help="Number of storeys")
    return p

def main():
    parser = build_parser()
    args = parser.parse_args()

    inputs = Inputs(
        volume=args.volume,
        n=args.n,
        n50=args.n50,
        eave_height=args.eave_height,
        indoor_temp=args.indoor_temp,
        outdoor_temp=args.outdoor_temp,
        temp_unit=args.temp_unit,
        flue_diam_mm=args.flue_diam_mm,
        flue_height=args.flue_height,
        wind_speed=args.wind_speed,
        wind_speed_met=args.wind_speed_met,
        met_height=args.met_height,
        terrain_class_met=args.terrain_class_met,
        terrain_class_site=args.terrain_class_site,
        shelter_walls=args.shelter_walls,
        shelter_flue=args.shelter_flue,
        house_type=args.house_type,
        foundation=args.foundation,
        storeys=args.storeys
    )

    res = compute_infiltration(inputs)

    print("AIM-2 Infiltration Results (Bradley HOT2000 + L. Lew fractions)")
    print("---------------------------------------------------------------")
    for k, v in res.items():
        print(f"{k}: {v:.6g}")

if __name__ == "__main__":
    main()
