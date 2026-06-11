# rescore.py
import json
from pathlib import Path
from calibrate_meltpool import peak_metrics, load_config, write_summary, RUNS_DIR, SUMMARY_CSV

config = load_config(Path("calibration_config.json"))

rows = []
for case_dir in sorted(RUNS_DIR.iterdir()):
    if not case_dir.is_dir():
        continue
    try:
        metrics = peak_metrics(case_dir, config)
        radius_um = float(case_dir.name.split("um_")[0].lstrip("r"))
        er = float(case_dir.name.split("er")[1].replace("e-6", "")) * 1e-6
        row = {
            "case_id": case_dir.name,
            "laserRadius_m": radius_um * 1e-6,
            "laserRadius_um": radius_um,
            "elec_resistivity_ohm_m": er,
            "target_width_um": config["target_width_um"],
            "target_depth_um": config["target_depth_um"],
            **metrics,
            "case_dir": str(case_dir),
        }
        rows.append(row)
        print(f"{case_dir.name}: width={metrics['peak_width_um']:.1f} depth={metrics['peak_depth_um']:.1f} score={metrics['score']:.3f}")
    except Exception as e:
        print(f"SKIP {case_dir.name}: {e}")

write_summary(rows, config)
print(f"\nDone. Updated {SUMMARY_CSV}")
