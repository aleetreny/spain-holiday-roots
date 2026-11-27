"""Build figures in order with the coursework virtual environment and Python -I."""

import argparse
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))

SCRIPTS = {
    "1": "fig1_design", "2": "fig2_calendar", "3": "fig3_who_pulls",
    "4": "fig4_receiving_end", "S1": "figS1_network", "S2": "figS2_parents",
    "S3": "figS3_robustness", "S4": "figS4_method_bridge", "S5": "figS5_zero_cells",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figures", nargs="*", metavar="FIGURE", help="defaults to all nine figures")
    args = parser.parse_args()
    for number in args.figures or SCRIPTS:
        if number not in SCRIPTS:
            parser.error(f"unknown figure {number}")
        print(f"Building Fig. {number}", flush=True)
        runpy.run_path(str(HERE / f"{SCRIPTS[number]}.py"), run_name="__main__")


if __name__ == "__main__":
    main()
