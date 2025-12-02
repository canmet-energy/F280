from typing import Dict
import math
class aim2:
    def __init__(self,volume,n,n50,eave_height,house_type="detached",foundation="crawl",storeys="1",shelter_walls=1.0,shelter_flue=1.0,flue_diam_mm=0.0,flue_height=None,met_height=10.0,terrain_class_met=5,terrain_class_site=7):
        """
        AIM-2 Infiltration Model (Bradley HOT2000 Implementation)
        =========================================================

        Implements Bradley's HOT2000 version of AIM-2:
        - Leakage distribution (R, X, Y) using L. Lew (1993) Table for fractions
          (ceiling, floor, walls) selected by house type, foundation, and storeys.
        - Flue leakage coefficient (from diameter)
        - Stack and wind effect calculations
        - Shelter coefficient (Bradley version)
        - Superposition of stack and wind flows
        - ACH conversions

        Inputs:
        -----------------------
        - volume  : volume of the conditioned space (m³)
        - n   : Flow exponent (dimensionless). Typical range: 0.65–0.75 for houses
        - n50 : Air change rate at 50 Pa (ACH @ 50 Pa)
        - eave_height : Height from ground to eaves (m)
        - house_type : "detached", "semi-detached" or "semi"
        - foundation : "crawl", "slab", "shallow", or "full"
        - storeys : "1", "1.5", "2", "2.5", or "3"
        - shelter_walls : Wall shelter coefficient. Less than or equal to 1
        - shelter_flue : Flue shelter coefficient. 1 is flue rises above surrounding obstructions, overwise same as shelter_walls
        - flue_diam_mm : Flue diameter (mm). If no flue, set to 0
        - flue_height : Height of the flue top relative to ground (m). If no flue, set to None
        - met_height : Height at which wind observations were made (m)
        - terrain_class_met : Terrain class at the meteorlogical observation site (1 to 8)
        - terrain_class_site : Terrain class at the building site (1 to 8)

        """

        # Store attributes and type-check
        self.__volume: float = volume
        self.__n: float = n
        self.__n50: float = n50
        self.__eave_height: float = eave_height
        self.__flue_diam_mm: float = flue_diam_mm
        self.__flue_height: float = flue_height
        self.__met_height: float = met_height
        self.__terrain_class_met: int = terrain_class_met
        self.__terrain_class_site: int = terrain_class_site
        self.__shelter_walls: float = shelter_walls
        self.__shelter_flue: float = shelter_flue
        self.__house_type: str = house_type
        self.__foundation: str = foundation
        self.__storeys: str = storeys
        
        self.__setGlobals() # Set the global constants
        self.__checkInputs() # Check the inputs
        
        # Determine base flow coefficient
        self.__C_base = (self.__n50 * self.__volume / 3600.0) / (50.0 ** self.__n)
        # Determine the flow fractions
        self.__setFractions()
        # Determine leakage coefficients
        self.__setLeakageCoefficients()
        # Determine leakage fractions
        self.__setLeakageFractions()
        # Shelter coefficient
        self.__setShelterCoefficient()
        
        if self.__flue_diam_mm > 0.0:
            #self.__hasFlue = True
            self.__Bf = self.__flue_height/self.__eave_height
            self.__fs_has_flue()
            self.__fw_has_flue()
        else:
            #self.__hasflue = False
            self.__fs_no_flue() # Stack factor (no-flue)
            self.__fw_no_flue()
    
    def getInfiltration(self,indoor_temp=22.0,outdoor_temp=0.0,wind_speed=0,temp_unit='C',is_met_wind_speed=True):
        """
        Inputs
        =========================================================
        - indoor_temp : Indoor drybulb temperature (C or K)
        - outdoor_temp : Outdoor drybulb temperature (C or K)
        - wind_speed : Wind speed (km/h)
        - temp_unit : Temperature units, either 'C' or 'K'
        - is_met_wind_speed : Boolean. True if entered windspeed was measured at meteorlogical site

        """
        if temp_unit != 'C' and temp_unit != 'K':
            raise Exception(f"Temperature unit entered is {temp_unit}. Must 'C' or 'K'.")
        if temp_unit == 'C':
            Ti = indoor_temp+273.15
            To = outdoor_temp+273.15
        else:
            Ti = indoor_temp
            To = outdoor_temp
        
        # Adjust temperatures (if required)
        Ps = self.__stack_pressure(Ti, To)
        Qs = self.__stack_flow(Ps)

        # Wind flow
        Ue = self.__effective_wind_speed(wind_speed,is_met_wind_speed) #converted to m/s
        Qw = self.__wind_flow(Ue)
        
        # Superposition
        Qnat = self.__superpose(Qs, Qw)
        
        return {
            "C_base": self.__C_base, "C_total": self.__C_total, "C_flue": self.__Cflue,
            "Cc0": self.__Cc0, "Cf0": self.__Cf0, "Cw0": self.__Cw0,
            "R": self.__R, "X": self.__X, "Y": self.__Y,
            "Ps": Ps, "fs": self.__fs,
            "Ue": Ue, "Sw":self.__Sw,
            "Qs": Qs, "Qw": Qw,
            "Qnat": Qnat
        }
    
    #################################################################
    #                       PRIVATE METHODS
    #################################################################
    
    # INITIALIZATION METHODS
    def __checkInputs(self):
        if self.__volume <= 0.0:
            raise Exception(f"Volume entered is {self.__volume}. Must be greater than zero")
        if self.__n < 0.45 or self.__n > 1.0:
            raise Exception(f"Flow exponent entered is {self.__n}. Must be greater than 0.45 or less than or equal to 1.")
        if self.__n50 <= 0.0:
            raise Exception(f"ACH50 entered is {self.__n50}. Must be greater than zero")
        # TODO: More error checking
    
    def __setGlobals(self):
        self.__RHO_AIR = 1.204097  # kg/m³ at ~20°C
        self.__G = 9.80665         # m/s²
        self.__AP_REF = 4.0        # Pa reference pressure for flue coefficient
        self.__B_INTERACTION = -1.0 / 3.0  # Bradley empirical interaction term
        self.__DAVENPORT_ZO: Dict[int, float] = { # Terrain roughness (Davenport)
            1: 0.0002,  # Open sea
            2: 0.005,   # Mud flats
            3: 0.03,    # Open flat terrain (grass)
            4: 0.10,    # Low crops
            5: 0.25,    # High crops
            6: 0.50,    # Parkland/bushes
            7: 1.00,    # Suburb/forest
            8: 2.00     # City centre (assumed)
        }
        # Leakage distribution tables (L. Lew, 1993)
        self.__DETACHED_ROWS = [
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
        self.__SEMI_ROWS = [
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

        self.__HOUSE_TYPE_MAP = {
            "detached": 1,
            "semi-detached": 2,
            "semi": 2
        }

        self.__FOUNDATION_MAP = {
            "crawl": 1,
            "slab": 2,        # slab-on-grade
            "shallow": 3,
            "full": 4
        }

        self.__STOREYS_MAP = {
            "1": 1, "1.0": 1,
            "1.5": 2,
            "2": 3, "2.0": 3,
            "2.5": 4,
            "3": 5, "3.0": 5
        }
    
    def __setFractions(self):
        """"
        Leakage distribution source:
        ----------------------------
        L. Lew, "Evaluation of AIM-2", EMR, Apr 23/93.
        Fortran table mapping (ceiling, floor, walls) fractions by:
          i = house type (1 detached, 2 semi-detached)
          j = foundation (1 crawl, 2 slab-on-grade, 3 shallow, 4 full)
          k = storeys (1=1, 2=1.5, 3=2, 4=2.5, 5=3)
        """
        i = self.__HOUSE_TYPE_MAP.get(self.__house_type.lower())
        j = self.__FOUNDATION_MAP.get(self.__foundation.lower())
        k = self.__STOREYS_MAP.get(self.__storeys.strip().lower())

        if i is None:
            raise ValueError("Invalid house_type. Use 'detached' or 'semi-detached'.")
        if j is None:
            raise ValueError("Invalid foundation. Use 'crawl', 'slab', 'shallow', or 'full'.")
        if k is None:
            raise ValueError("Invalid storeys. Use one of '1','1.5','2','2.5','3'.")

        rows = self.__DETACHED_ROWS if i == 1 else self.__SEMI_ROWS
        self.__fractions = rows[j - 1][k - 1]  # tuples sum to 1.0
    
    def __setLeakageCoefficients(self):
        a_c, a_w, a_f = self.__fractions
        # Normalize (safety) and distribute base C
        s = a_c + a_f + a_w
        if s <= 0.0:
            raise ValueError("Leakage fractions must sum > 0")
        a_c, a_f, a_w = a_c / s, a_f / s, a_w / s

        self.__Cc0 = a_c * self.__C_base
        self.__Cf0 = a_f * self.__C_base
        self.__Cw0 = a_w * self.__C_base

        # Flue coefficient
        self.__Cflue = 0.0
        if self.__flue_diam_mm and self.__flue_diam_mm > 0.0:
            d_m = self.__flue_diam_mm / 1000.0
            area = math.pi * (d_m ** 2) / 4.0
            # Bradley: C_flue = 0.5 * A_flue * (rho / (2 * ΔP_ref))^(n - 0.5)
            self.__Cflue = 0.5 * area * ((self.__RHO_AIR / (2.0 * self.__AP_REF)) ** (self.__n - 0.5))

        self.__C_total = self.__Cc0 + self.__Cf0 + self.__Cw0 + self.__Cflue
    
    def __setLeakageFractions(self):
        a_c, a_w, a_f = self.__fractions
        self.__Y = (self.__Cflue / self.__C_total) if self.__C_total > 0.0 else 0.0
        self.__R = (a_c + a_f) * (1.0 - self.__Y)
        self.__X = (a_c - a_f) * (1.0 - self.__Y)
    
    def __setShelterCoefficient(self):
        #This equation is shown as Swo * (1.0 + Y) + Swflue * (1.5 * Y) in Bradley but this is a typo
        self.__Sw=self.__shelter_walls * (1.0 - self.__Y) + self.__shelter_flue * (1.5 * self.__Y)
    
    def __fs_no_flue(self):
        """
        Bradley approximation for stack factor (no flue).
        """
        self.__fs = ((1 + self.__n * self.__R) / (self.__n + 1)) * ((0.5 - 0.5 * (((self.__X ** 2)/(2 - self.__R)) ** (5 / 4)))) ** (self.__n + 1)
    
    def __fs_has_flue(self):
        """
        Bradley approximation for stack factor (with flue).
        """
        M=((self.__X+((2.0*self.__n)+1.0)*self.__Y)**2.0)/(2.0-self.__R)
        if M > 1.0: # Limit as defined in Walker and Wilson Equation 23
            M = 1.0
        Xc = self.__R + ((2.0*(1.0-self.__R-self.Y))/(self.__n+1.0))-(2.0*self.__Y*((self.__Bf-1.0)**self.__n))
        
        exp = ((3.0*self.__n)-1.0)/3.0
        num = 3.0*((Xc-self.__X)**2.0)*(self.__R**(1.0-self.__n))
        denom = 2.0*(self.__Bf+1.0)
        F = self.__n*self.__Y*((self.__Bf-1.0)**exp)*(1.0-(num/denom))
        self.__fs = ((1.0+(self.__n*self.__R))/(self.__n+1.0))*((0.5-(0.5*(M**(5.0/4.0))))**(self.__n+1))+F

    def __fw_no_flue(self):
        self.__fw = 0.19 * (2.0 - self.__n) * (1.0 - ((self.__X + self.__R) / 2) ** 1.5)
    
    def __fw_has_flue(self):
        S = (self.__X+self.__R+(2.0*self.__Y))/2.0
        comp1 = 1.0-(((self.__X+self.__R)/2.0)**((3.0/2.0)-self.__Y))
        self.__fw = (0.19*(2.0-self.__n)*comp1)-((self.__Y/4.0)*(S-(2.0*self.__Y*(S**4.0))))
        
    # CALCULATIONS METHODS
    def __stack_pressure(self,Ti,To):
        """
        Bradley: Ps = rho * g * H * (Ti - To) / Ti
        """
        return self.__RHO_AIR * self.__G * self.__eave_height * (abs(Ti - To) / Ti)
    
    def __stack_flow(self,Ps):
        C_converted = self.__C_total * 3600/self.__volume
        return C_converted * self.__fs * (Ps ** self.__n)
    
    def __effective_wind_speed(self,wind_speed,bIsMet):
        """
        Compute site windspeed:
        - If explicit wind_speed provided, use it.
        - Else adjust met tower wind to site using log-profile approximation (Davenport classes).
        """
        if not bIsMet:
            return wind_speed / 3.6  # Convert km/h to m/s
        Zo_met = self.__DAVENPORT_ZO.get(self.__terrain_class_met, 0.25)
        Zo_site = self.__DAVENPORT_ZO.get(self.__terrain_class_site, 0.25)
        H = max(self.__eave_height, 1.0)
        h_met = max(self.__met_height, 1.0)
        num = math.log(H / Zo_site)*math.log(60/Zo_met)
        den = math.log(h_met / Zo_met)*math.log(60/Zo_site)
        if den <= 0.0:
            den = 1.0
        return wind_speed * (num / den) / 3.6  # Convert km/h to m/s
    
    def __wind_flow(self,Ue):
        C_converted = self.__C_total * 3600/self.__volume
        # Bradley: Pw = 0.5 * rho * (Sw * Ue)^2
        Pw = 0.5 * self.__RHO_AIR * (self.__Sw * Ue) ** 2  
        return C_converted * self.__fw * (Pw ** self.__n)
    
    def __superpose(self,Qs,Qw):
        """
        Bradley superposition:
        Q_nat = [ Qs^n + Qw^n + B * (Qs * Qw)^n ]^(1/n),  B ≈ -1/3
        """
        return ((Qs ** (1/self.__n)) + (Qw ** (1/self.__n)) + self.__B_INTERACTION * ((Qs * Qw) ** (0.5 / self.__n))) ** self.__n