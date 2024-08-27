import argparse
import os

import danton


def main():
    """Entry point for the CLI."""

    parser = argparse.ArgumentParser(
        prog = "python3 -m danton",
        description = "Command-line utility for Danton.",
        epilog = "Copyright (C) Université Clermont Auvergne, CNRS/IN2P3, LPCA"
    )
    parser.add_argument("-c", "--cache",
        help = "Danton default cache location.",
        action = "store_true"
    )
    parser.add_argument("-p", "--prefix",
        help = "Danton installation prefix.",
        action = "store_true"
    )
    parser.add_argument("-v", "--version",
        help = "Danton version.",
        action = "store_true"
    )

    args = parser.parse_args()

    result = []
    if args.cache:
        result.append(danton.DEFAULT_CACHE)
    if args.prefix:
        result.append(os.path.dirname(__file__))
    if args.version:
        result.append(danton.VERSION)

    if result:
        print(" ".join(result))


if __name__ == "__main__":
    main()
