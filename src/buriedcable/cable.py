"""Cable physical parameters and internal thermal resistances (IEC 60287-1-1)."""
import math
from dataclasses import dataclass


@dataclass
class CableParams:
    """
    Single-core XLPE power cable.
    Defaults: 240 mm² Cu, 110 kV, XLPE insulation, PE outer sheath.
    """
    # --- Conductor ---
    conductor_area_mm2: float = 240.0
    conductor_material: str = "Cu"      # "Cu" or "Al"
    conductor_diameter_mm: float = 18.0

    # --- Insulation (XLPE) ---
    insulation_thickness_mm: float = 15.0
    insulation_resistivity_Km_W: float = 3.5   # K·m/W

    # --- Outer sheath (PE) ---
    sheath_thickness_mm: float = 3.5
    sheath_resistivity_Km_W: float = 3.5       # K·m/W

    # --- Operating limits ---
    voltage_kV: float = 110.0
    max_conductor_temp_C: float = 90.0          # XLPE continuous rating

    # --- Installation ---
    burial_depth_m: float = 1.0                 # depth to cable centre (m)

    # --- Temperature coefficient of resistance ---
    _ALPHA = {"Cu": 3.93e-3, "Al": 4.03e-3}    # 1/K at 20°C
    _RHO20 = {"Cu": 1.7241e-8, "Al": 2.8264e-8}  # Ω·m at 20°C

    # ------------------------------------------------------------------ #
    # Derived geometry                                                     #
    # ------------------------------------------------------------------ #

    @property
    def external_diameter_mm(self) -> float:
        return (self.conductor_diameter_mm
                + 2 * self.insulation_thickness_mm
                + 2 * self.sheath_thickness_mm)

    # ------------------------------------------------------------------ #
    # Electrical resistance                                                #
    # ------------------------------------------------------------------ #

    def dc_resistance_ohm_per_m(self, temp_C: float = 20.0) -> float:
        """IEC 60228 DC resistance at given temperature (Ω/m)."""
        rho20 = self._RHO20[self.conductor_material]
        alpha = self._ALPHA[self.conductor_material]
        area_m2 = self.conductor_area_mm2 * 1e-6
        return (rho20 / area_m2) * (1.0 + alpha * (temp_C - 20.0))

    def ac_resistance_ohm_per_m(self, temp_C: float = 90.0) -> float:
        """
        AC resistance at operating temperature.
        Skin/proximity factors neglected (< 2% for 240 mm² at 50 Hz).
        """
        return self.dc_resistance_ohm_per_m(temp_C)

    # ------------------------------------------------------------------ #
    # Internal thermal resistances (IEC 60287-2-1)                        #
    # ------------------------------------------------------------------ #

    def t1_insulation(self) -> float:
        """T1: insulation thermal resistance (K·m/W)."""
        d_c = self.conductor_diameter_mm
        t1 = self.insulation_thickness_mm
        rho = self.insulation_resistivity_Km_W
        return (rho / (2.0 * math.pi)) * math.log(1.0 + 2.0 * t1 / d_c)

    def t3_sheath(self) -> float:
        """T3: sheath thermal resistance (K·m/W)."""
        d_i = self.conductor_diameter_mm + 2.0 * self.insulation_thickness_mm
        t3 = self.sheath_thickness_mm
        rho = self.sheath_resistivity_Km_W
        return (rho / (2.0 * math.pi)) * math.log(1.0 + 2.0 * t3 / d_i)

    def t_cable(self) -> float:
        """Total internal resistance T_cable = T1 + T3 (K·m/W)."""
        return self.t1_insulation() + self.t3_sheath()
