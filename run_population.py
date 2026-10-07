import argparse
from contextlib import chdir
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_DIR / "config" / "sample_population.ini"
OUTPUT_DIR = PROJECT_DIR / "outputs"


def main():
    parser = argparse.ArgumentParser(description="Run a POSYDON population simulation.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    config_path = args.config.resolve()

    from posydon.popsyn.synthetic_population import PopulationRunner

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # POSYDON writes batch files and the merged population into the working directory.
    with chdir(OUTPUT_DIR):
        poprun = PopulationRunner(str(config_path), verbose=True)
        print("Number of binary populations:", len(poprun.binary_populations))
        print("Metallicity:", poprun.binary_populations[0].metallicity)
        print("Number of binaries:", poprun.binary_populations[0].number_of_binaries)
        poprun.evolve()


if __name__ == "__main__":
    main()
