#! /usr/bin/env python3
import argparse
import gzip
import numpy as np
from pathlib import Path
import pickle
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x: x


def merge(args):
    """Merge simulation results."""

    def get_tag(path):
        return path.name.rsplit("_", 1)[0]

    tag = get_tag(args.files[0])

    data = None
    for path in tqdm(args.files):
        t = get_tag(path)
        if t != tag:
            raise ValueError(f"bad tag (expected '{tag}', found '{t}')")

        with gzip.open(path, "rb") as f:
            d = pickle.load(f)
            n = d["secondaries"].size

        if data is None:
            data = d
            data["seed"] = n * [data["seed"]]
        else:
            data["events"] += d["events"]
            for k in ("secondaries", "random_index"):
                data[k] = np.append(data[k], d[k])
            data["seed"] += n * [d["seed"]]
            data["cpu"] += d["cpu"]

    outfile = args.files[0].parent / f"{tag}.pkl.gz"
    with gzip.open(outfile, "wb") as f:
        pickle.dump(data, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge simulation results")
    parser.add_argument("files",
        help = "Input file(s)",
        nargs = "+",
        type = lambda x: Path(x)
    )

    args = parser.parse_args()
    merge(args)
