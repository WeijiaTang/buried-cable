"""
Soil thermal model.

Johansen (1975) formula for thermal conductivity vs. saturation degree,
IEC 60287-2-1 buried cable external thermal resistance (T4),
and soil drying-out threshold.
"""
import math
from dataclasses import dataclass
from enum import Enum


class SoilType(Enum):
    SAND    = "sand"
    SILT    = "silt"
    CLAY    = "clay"
    LOAM    = "loam"


# Johansen parameters per soil type
# (lambda_dry, lambda_sat_unfrozen, n_porosity, coarse_grained)
# coarse_grained: True  → Ke = log10(Sr) + 1        (sand)
#                 False → Ke = 0.7·log10(Sr) + 1     (silt, clay, loam)
_JOHANSEN_PARAMS: dict[SoilType, tuple[float, float, float, bool]] = {
    #                   λ_dry  λ_sat   n      coarse
    SoilType.SAND:  (0.30,  2.20,  0.40,  True),
    SoilType.SILT:  (0.25,  1.80,  0.45,  False),
    SoilType.CLAY:  (0.20,  1.50,  0.50,  False),
    SoilType.LOAM:  (0.25,  1.70,  0.43,  False),
}

# IEC 60287-2-1 §2.2.3 drying-out critical heat flux (W/m).
# Derived from the cable-surface temperature reaching the soil's critical
# drying temperature T_crit (≈50 °C for moist soil): q_c = (T_crit − T_amb)/T4.
# With T_amb ≈ 15 °C and T4_moist ≈ 0.65 K·m/W this gives q_c ≈ 54 W/m for
# sand; fine-grained soils are assigned a lower threshold because they crack
# and lose capillary continuity at a lower surface temperature. Below q_c the
# soil stays moist; above it sustained moisture migration dries the cable
# vicinity and the dry-zone resistivity takes over.
_DRYOUT_FLUX_W_PER_M: dict[SoilType, float] = {
    SoilType.SAND:  55.0,
    SoilType.SILT:  45.0,
    SoilType.CLAY:  35.0,
    SoilType.LOAM:  50.0,
}


@dataclass
class SoilParams:
    soil_type: SoilType = SoilType.SAND
    burial_depth_m: float = 1.0           # depth to cable centre (m)
    cable_ext_diameter_m: float = 0.075   # outer diameter of installed cable (m)

    # Two-region model: moist_rho and dry_rho in K·m/W (used as override)
    # When set to None the Johansen model is used dynamically.
    static_rho_moist: float | None = None   # K·m/W, IEC 60287 "T4 moist zone"
    static_rho_dry:   float | None = None   # K·m/W, IEC 60287 "T4 dry zone"

    def _johansen_lambda(self, moisture_vol: float) -> float:
        """
        Johansen (1975): soil thermal conductivity (W/m·K) from
        volumetric moisture content (m³/m³).

        Kersten number:
          coarse-grained (sand):      Ke = log10(Sr) + 1
          fine-grained (silt/clay):   Ke = 0.7·log10(Sr) + 1
        """
        lam_dry, lam_sat, n, coarse = _JOHANSEN_PARAMS[self.soil_type]
        Sr = moisture_vol / n
        Sr = max(0.0, min(Sr, 1.0))
        if Sr <= 0.0:
            return lam_dry
        Ke = math.log10(Sr) + 1.0 if coarse else 0.7 * math.log10(Sr) + 1.0
        Ke = max(0.0, Ke)
        return Ke * (lam_sat - lam_dry) + lam_dry

    def thermal_resistivity_Km_W(self, moisture_vol: float) -> float:
        """
        Effective soil thermal resistivity (K·m/W) at given volumetric
        moisture content.  Returns static_rho_moist if override is set.
        """
        if self.static_rho_moist is not None:
            return self.static_rho_moist
        lam = self._johansen_lambda(moisture_vol)
        return 1.0 / lam

    def is_dryout(self, heat_flux_W_per_m: float) -> bool:
        """True if cable linear heat flux exceeds drying-out threshold."""
        return heat_flux_W_per_m > _DRYOUT_FLUX_W_PER_M[self.soil_type]

    def effective_resistivity(
        self,
        moisture_vol: float,
        heat_flux_W_per_m: float,
    ) -> float:
        """
        Return effective soil thermal resistivity accounting for drying-out.
        If drying-out occurs, use dry-zone resistivity (2× moist as fallback).
        """
        rho_moist = self.thermal_resistivity_Km_W(moisture_vol)
        if self.is_dryout(heat_flux_W_per_m):
            if self.static_rho_dry is not None:
                return self.static_rho_dry
            # Dry resistivity ≈ lambda_dry^-1
            lam_dry = _JOHANSEN_PARAMS[self.soil_type][0]
            return 1.0 / lam_dry
        return rho_moist

    # ------------------------------------------------------------------ #
    # IEC 60287-2-1 external thermal resistance T4                        #
    # ------------------------------------------------------------------ #

    def t4_external(self, rho_soil_Km_W: float) -> float:
        """
        T4: soil external thermal resistance for a single buried cable
        using IEC 60287-2-1 image method (K·m/W).

            T4 = (rho / 2π) · ln(L/r + sqrt((L/r)² - 1))
                ≈ (rho / 2π) · ln(2L / r)    when L >> r
        """
        L = self.burial_depth_m
        r = self.cable_ext_diameter_m / 2.0
        u = L / r
        return (rho_soil_Km_W / (2.0 * math.pi)) * math.log(u + math.sqrt(u**2 - 1.0))

    def t4_dynamic(
        self,
        moisture_vol: float,
        heat_flux_W_per_m: float,
    ) -> float:
        """T4 computed from dynamic Johansen resistivity."""
        rho = self.effective_resistivity(moisture_vol, heat_flux_W_per_m)
        return self.t4_external(rho)
