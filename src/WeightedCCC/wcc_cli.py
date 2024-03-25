# -*- coding: utf-8 -*-
"""CLI implementation for the weighted-ccc package."""

import argparse

from .graph import GraphClosure
from .calc_lig import (
    cal_node_path_dependent_error,
    cal_node_path_independent_error,
    calcMolEnes,
    getMolEnergyDataFrame,
    set_node_map,
)


# if weights contained bennett_std, pleas display bennett_std before other weights
def parse_arguments(fakeArgs=None):
    parser = argparse.ArgumentParser(description="Process energy data for weighted CCC calculation.")

    parser.add_argument("-f", "--file", dest="file", help="Input file containing pairwise energy data", required=True)
    parser.add_argument(
        "-r",
        "--ref",
        dest="ref",
        default="",
        help=(
            "Reference molecule for calculating the energy for other molecules. "
            "Default: The first molecule in the data file"
        ),
    )
    parser.add_argument(
        "-e",
        "--ref_ene",
        dest="ref_ene",
        type=float,
        default=0.00,
        help="Energy for the reference molecule. Default: 0.00",
    )
    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        choices=["no", "yes"],
        default="no",
        help=(
            "Print option: no(Only print molecule energy), "
            "yes(print pair-wise energy and molecule energy). Default: no"
        ),
    )

    if fakeArgs:
        return parser.parse_args(fakeArgs)
    else:
        return parser.parse_args()


def main():
    args = parse_arguments()
    # fakeArgs = "-f bace_run1_0_with_w -r 3A -e -8.83 -p yes"  # only keep this for test purpose
    # opts = optParser(fakeArgs.strip().split())  # only keep this for test purpose
    args.print = args.print.lower()
    if not args.file:
        raise Exception("No input energy data!")
    g = GraphClosure(filename=args.file)
    g.getAllCyles()
    if len(g.cycles) == 0:
        print("No cycle in this graph.")
        exit()
    g.iterateCycleClosure(minimum_cycles=2)
    if args.print == "yes":
        g.getEnergyPairsDataFrame(verbose=True)
    node_map = set_node_map(g)
    if not args.ref.strip():
        args.ref = g.V[0]
    try:
        ref_node = g.V.index(args.ref)
    except ValueError as err:
        print(f"Check your args. Ref {args.ref} isn't in your input file! Error:\n{err}")
        raise ValueError from err
    path_independent_error = cal_node_path_independent_error(g.V, node_map)
    path_dependent_error, path = cal_node_path_dependent_error(ref_node, g.V, node_map)
    mol_ene = calcMolEnes(args.ref_ene, g, path)
    getMolEnergyDataFrame(g.V, mol_ene, path_dependent_error, path_independent_error, verbose=True)


if __name__ == "__main__":
    main()
