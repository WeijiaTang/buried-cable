# Buried Cable Dynamic Ampacity Model

Layered Python implementation of a dynamic ampacity model for buried power
cables: IEC 60287 thermal circuit with a Johansen moisture-dependent soil
thermal resistivity, a heat-flux drying-out criterion, and a recursive
least-squares (RLS) identifier for online estimation of the external thermal
resistance from cable-surface temperature.

## Contents

```
src/buriedcable/        core model package
  cable.py              conductor/insulation/sheath geometry, R_ac, T1+T3
  soil.py               Johansen lambda(Sr), drying-out threshold, T4
  thermal.py            static_ampacity / dynamic_ampacity
  climate.py            climate CSV parser, depth interpolation
  identification.py     RLS online T4 estimation
data/raw/               hourly climate records (Open-Meteo, 8760 rows each)
data/processed/         derived seasonal correction factors
```

## Environment

Python 3.13, managed with [`uv`](https://docs.astral.sh/uv/). Dependencies:
`numpy`, `scipy`, `matplotlib`.

```bash
uv sync
uv run python -c "from buriedcable import ThermalCircuit, CableParams, SoilParams, SoilType; \
tc=ThermalCircuit(CableParams(), SoilParams(SoilType.SAND, 1.0, 0.055)); \
print(tc.dynamic_ampacity(15.0, 0.20))"
```

The package is not installed; import via `sys.path.insert(0, .../src)`.

## Data source

Hourly soil temperature and moisture from the
[Open-Meteo](https://open-meteo.com/) historical archive (free, no API key).

## License

MIT License. Climate data subject to the Open-Meteo terms of use.
