"""Load and preprocess climate CSV files produced by scripts/fetch_climate.py."""
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class HourlyRecord:
    datetime: str
    temperature_2m: float
    soil_temperature_0_7cm: float
    soil_moisture_0_7cm: float
    soil_moisture_7_28cm: float

    @property
    def ambient_temp_C(self) -> float:
        return self.temperature_2m

    def soil_moisture_at_depth(self, depth_m: float) -> float:
        """
        Linear interpolation of volumetric moisture between the two layers.
        Layer 0–0.07 m → soil_moisture_0_7cm
        Layer 0.07–0.28 m → soil_moisture_7_28cm
        Beyond 0.28 m extrapolates with the deeper value.
        """
        if depth_m <= 0.07:
            return self.soil_moisture_0_7cm
        if depth_m >= 0.28:
            return self.soil_moisture_7_28cm
        t = (depth_m - 0.07) / (0.28 - 0.07)
        return (1.0 - t) * self.soil_moisture_0_7cm + t * self.soil_moisture_7_28cm


def load_climate(path: Path) -> list[HourlyRecord]:
    """Parse a climate CSV and return a list of HourlyRecord."""
    records: list[HourlyRecord] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                records.append(HourlyRecord(
                    datetime=row["datetime"],
                    temperature_2m=float(row["temperature_2m"]),
                    soil_temperature_0_7cm=float(row["soil_temperature_0_to_7cm"]),
                    soil_moisture_0_7cm=float(row["soil_moisture_0_to_7cm"]),
                    soil_moisture_7_28cm=float(row["soil_moisture_7_to_28cm"]),
                ))
            except (ValueError, KeyError):
                # Skip rows with missing sensor data (rare in Open-Meteo)
                continue
    return records
