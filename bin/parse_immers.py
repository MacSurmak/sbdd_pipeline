#!/usr/bin/env python3
"""
bin/parse_immers.py
Parses PPM 3.0 datapar1 output file into JSON and text thickness deliverables.
"""
import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout)

def main():
    parser = argparse.ArgumentParser(description="Parse PPM 3.0 immersion datapar1 output")
    parser.add_argument("--datapar", default="datapar1", help="Path to datapar1 file")
    parser.add_argument("--out_json", required=True, help="Output JSON path")
    parser.add_argument("--out_txt", required=True, help="Output text path (thickness only)")
    args = parser.parse_args()

    data = {"thickness": 28.0, "tilt": 0.0, "dg": 0.0}
    datapar_path = Path(args.datapar)

    if datapar_path.exists():
        try:
            with open(datapar_path, "r") as f:
                parts = [p.strip() for p in f.readline().split(";")]
                data["thickness"] = float(parts[1])
                data["tilt"] = float(parts[3])
                data["dg"] = float(parts[5])
            logging.info(f"Parsed PPM 3.0: thickness={data['thickness']}, tilt={data['tilt']}, dG={data['dg']}")
        except Exception as e:
            logging.warning(f"Failed to parse datapar1 ({e}); using defaults: {data}")
    else:
        logging.warning(f"{datapar_path} not found; using defaults: {data}")

    with open(args.out_json, "w") as f_json:
        json.dump(data, f_json, indent=2)

    with open(args.out_txt, "w") as f_txt:
        f_txt.write(str(data["thickness"]))

if __name__ == "__main__":
    main()

