"""
IEC 60287 thermal circuit model.

Provides both static (standard) and dynamic ampacity calculation.
The thermal circuit per-unit-length (K·m/W) is:

  T_total = T_cable + T4
  I_max   = sqrt[(T_conductor_max - T_ambient) / (R_ac · T_total)]
"""
import math
from dataclasses import dataclass

from .cable import CableParams
from .soil import SoilParams


@dataclass
class ThermalCircuit:
    cable: CableParams
    soil: SoilParams

    # ------------------------------------------------------------------ #
    # Static IEC 60287 ampacity                                           #
    # ------------------------------------------------------------------ #

    def static_ampacity(
        self,
        T_ambient_C: float,
        rho_soil_Km_W: float,
    ) -> float:
        """
        IEC 60287 steady-state ampacity (A).
        rho_soil_Km_W : soil thermal resistivity (constant, as per standard)
        """
        delta_T = self.cable.max_conductor_temp_C - T_ambient_C
        R_ac    = self.cable.ac_resistance_ohm_per_m(self.cable.max_conductor_temp_C)
        T_cable = self.cable.t_cable()
        T4      = self.soil.t4_external(rho_soil_Km_W)
        T_total = T_cable + T4

        if delta_T <= 0.0 or T_total <= 0.0:
            return 0.0
        return math.sqrt(delta_T / (R_ac * T_total))

    # ------------------------------------------------------------------ #
    # Dynamic ampacity                                                     #
    # ------------------------------------------------------------------ #

    def dynamic_ampacity(
        self,
        T_ambient_C: float,
        moisture_vol: float,
        current_estimate_A: float | None = None,
        load_factor: float = 0.65,
    ) -> tuple[float, float, bool]:
        """
        Dynamic ampacity using Johansen soil thermal resistivity.

        Drying-out is evaluated against the *sustained operating* heat flux
        ``q = I_op² · R_ac`` (W/m), where the operating current ``I_op`` is
        ``load_factor`` times the IEC 60287 reference ampacity (rho = 1.0
        K·m/W), or ``current_estimate_A`` if supplied. This follows IEC
        60287-2-1 §2.2.3: drying is driven by the actual thermal load on the
        soil, not by the computed ampacity limit. When ``q`` exceeds the
        soil's critical flux the dry-zone resistivity governs; otherwise the
        Johansen moist resistivity at the current moisture content is used.

        Parameters
        ----------
        T_ambient_C        : ambient temperature (°C)
        moisture_vol       : volumetric soil moisture (m³/m³)
        current_estimate_A : explicit operating current for drying-out check.
                             If None, estimated as load_factor × IEC ampacity.
        load_factor        : fraction of IEC reference ampacity used as the
                             operating current (default 0.65).

        Returns
        -------
        (I_max, rho_eff, dryout_flag)
        """
        R_ac    = self.cable.ac_resistance_ohm_per_m(self.cable.max_conductor_temp_C)
        T_cable = self.cable.t_cable()

        if current_estimate_A is not None:
            I_op = current_estimate_A
        else:
            I_op = load_factor * self.static_ampacity(T_ambient_C, 1.0)

        heat_flux = I_op**2 * R_ac
        dryout    = self.soil.is_dryout(heat_flux)
        rho_eff   = self.soil.effective_resistivity(moisture_vol, heat_flux)

        T4      = self.soil.t4_external(rho_eff)
        T_total = T_cable + T4
        delta_T = self.cable.max_conductor_temp_C - T_ambient_C

        if delta_T <= 0.0 or T_total <= 0.0:
            return 0.0, rho_eff, dryout

        I_max = math.sqrt(delta_T / (R_ac * T_total))
        return I_max, rho_eff, dryout

    # ------------------------------------------------------------------ #
    # Conductor temperature at given current                              #
    # ------------------------------------------------------------------ #

    def conductor_temperature(
        self,
        current_A: float,
        T_ambient_C: float,
        rho_soil_Km_W: float,
    ) -> float:
        """Steady-state conductor temperature at given current (°C)."""
        R_ac    = self.cable.ac_resistance_ohm_per_m()
        T_cable = self.cable.t_cable()
        T4      = self.soil.t4_external(rho_soil_Km_W)
        return T_ambient_C + current_A**2 * R_ac * (T_cable + T4)
