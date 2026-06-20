"""
Online Recursive Least Squares (RLS) identification of soil thermal resistance.

Two observable temperatures are supported:

  (a) Conductor temperature (from transient thermal model or FEM):
        T_cond = T_amb + I² · R_ac · (T_cable + T4)
      Observation: y = T_cond - T_amb,  φ = I² · R_ac,  unknown θ = T_cable + T4
      → T4 = θ_est - T_cable

  (b) Cable outer-surface temperature (from DTS or thermistor):
        T_surf ≈ T_amb + I² · R_ac · T4          (sheath drop already excluded)
      Observation: y = T_surf - T_amb,  φ = I² · R_ac,  unknown θ = T4

Default mode is (b) because DTS fibre-optic is the practical sensor.
"""
from dataclasses import dataclass, field


@dataclass
class RLSIdentifier:
    """
    Scalar recursive least squares for T4 estimation.

    Parameters
    ----------
    init_t4        : initial guess for T4 (K·m/W)
    forgetting     : forgetting factor λ ∈ (0, 1].  1.0 = no forgetting.
    conductor_mode : if True, the temperature input is conductor temperature
                     and T_cable is subtracted from the estimate.
                     If False (default), input is outer-surface temperature
                     and T_cable is not needed.
    """
    init_t4:        float = 0.5
    forgetting:     float = 0.98
    conductor_mode: bool  = False

    _theta_est: float = field(init=False)  # estimated (T_cable+T4) or T4
    _P:         float = field(init=False)

    def __post_init__(self) -> None:
        self._theta_est = self.init_t4
        self._P         = 1.0

    def update(
        self,
        T_obs_C:    float,
        T_ambient_C: float,
        current_A:   float,
        R_ac:        float,
        T_cable:     float = 0.0,
    ) -> float:
        """
        Ingest one measurement and return updated T4 estimate (K·m/W).

        Parameters
        ----------
        T_obs_C     : observed temperature (°C).
                      conductor temperature if conductor_mode=True,
                      cable outer-surface temperature if False.
        T_ambient_C : ambient temperature (°C)
        current_A   : conductor current (A)
        R_ac        : AC resistance per unit length (Ω/m)
        T_cable     : internal cable thermal resistance (K·m/W).
                      Only used when conductor_mode=True.
        """
        phi = current_A**2 * R_ac
        y   = T_obs_C - T_ambient_C

        lam = self.forgetting
        K   = (self._P * phi) / (lam + phi * self._P * phi)
        self._theta_est += K * (y - phi * self._theta_est)
        self._P = (1.0 / lam) * (1.0 - K * phi) * self._P

        return self.t4(T_cable)

    def t4(self, T_cable: float = 0.0) -> float:
        """Current T4 estimate (K·m/W)."""
        if self.conductor_mode:
            return max(0.0, self._theta_est - T_cable)
        return max(0.0, self._theta_est)

    def reset(self, init_t4: float | None = None) -> None:
        self._theta_est = init_t4 if init_t4 is not None else self.init_t4
        self._P         = 1.0
