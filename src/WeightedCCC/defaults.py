# -*- coding: utf-8 -*-
"""module with DefaultGraph class, the base class for the GraphClosure class."""

from collections import defaultdict
import decimal


class DefaultGraph:  # DefaultGraph -> GraphClosure
    """Base class for the GraphClosure class. This class implementation contains operations
    that involve reading and preparing the input data. For the cycle closure correction algorithm,
    refer to the .graph.GraphClosure class
    """

    def __init__(self, from_lig=None, to_lig=None, b_ddG=None, weights=None, filename=None):
        # Initialize graph attributes
        self.initialize_attributes()

        if all([from_lig is not None, to_lig is not None, b_ddG is not None]):
            self.load_from_values(from_lig, to_lig, b_ddG, weights=weights)

        if filename is not None:
            self.load_from_file(filename)

        # Update vertices list
        self.V = list(self.graph.keys())

    def load_from_values(self, from_lig, to_lig, b_ddG, weights=None):
        results = self.initialize_data_structures()
        n_edges = 0
        if weights is None:
            weights = [None] * len(from_lig)
        for _from, _to, bddG, weight in zip(from_lig, to_lig, b_ddG, weights):
            weight_num = self.addEdge(_from, _to, bddG, weight, results)
            n_edges += 1
        results.update({"num_e": n_edges, "weight_num": weight_num})
        self.apply_results(results)

    def load_from_file(self, filename):
        results = DefaultGraph.from_file(filename)
        # Update this Graph instance with the loaded data
        self.apply_results(results)

    def initialize_attributes(self):
        self.num_e = 0
        self.graph = defaultdict(list)
        self.cycles = []
        self.nodelist = []
        self.paths = defaultdict(list)
        self.start = ""
        self.ddG = defaultdict(decimal.Decimal)
        self.ddG_cc = defaultdict(list)
        self.ddG_save = defaultdict(decimal.Decimal)
        self.err = defaultdict(decimal.Decimal)
        self.weight = defaultdict(list)
        self.print_e = defaultdict(decimal.Decimal)
        self.cycleset = set()

    def apply_results(self, results):
        # Assuming results is a dictionary of all relevant data structures
        # Update each relevant attribute of the Graph instance
        self.graph = results["graph"]
        self.ddG = results["ddG"]
        self.ddG_cc = results["ddG_cc"]
        self.ddG_save = results["ddG_save"]
        self.err = results["err"]
        self.weight = results["weight"]
        self.print_e = results["print_e"]
        self.num_e = sum(len(edges) for edges in self.graph.values()) // 2
        self.weight = results["weight"]
        self.weight_num = results["weight_num"]
        self.V = list(self.graph.keys())

    @classmethod
    def initialize_data_structures(cls):
        return {
            "graph": defaultdict(list),
            "ddG": defaultdict(decimal.Decimal),
            "ddG_cc": defaultdict(list),
            "ddG_save": defaultdict(decimal.Decimal),
            "err": defaultdict(decimal.Decimal),
            "weight": defaultdict(list),
            "print_e": defaultdict(decimal.Decimal),
            "weight_num": int(0),
        }

    @staticmethod
    def parse_line(line):
        data = line.strip().split()
        if len(data) < 3:
            return None
        from_lig, to_lig, b_ddG = data[:3]
        weights = data[3:] if len(data) > 3 else None
        return from_lig, to_lig, b_ddG, weights

    @classmethod
    def from_file(cls, filename):
        if not filename:
            raise ValueError("Filename must not be None")

        results = cls.initialize_data_structures()
        n_edges = 0
        try:
            with open(filename, "r") as file:
                for line in file:
                    parsed_line = cls.parse_line(line)
                    if not parsed_line:
                        continue  # Skip lines that don't match expected format
                    weight_num = cls.addEdge(*parsed_line, results)
                    n_edges += 1
        except FileNotFoundError as err:
            raise Exception(f"File {filename} not found.") from err
        except Exception as err:
            raise Exception(f"Error processing file {filename}: {str(err)}") from err
        results.update({"num_e": n_edges, "weight_num": weight_num})
        return results

    @classmethod
    def addEdge(cls, from_lig, to_lig, b_ddG, weights, results):
        # Unpack results data structures for easier access
        graph, ddG, ddG_cc, ddG_save, err, weight, print_e, weight_num = (
            results["graph"],
            results["ddG"],
            results["ddG_cc"],
            results["ddG_save"],
            results["err"],
            results["weight"],
            results["print_e"],
            results["weight_num"],
        )
        graph[from_lig].append(to_lig)
        graph[to_lig].append(from_lig)

        # Update delta delta G (ddG) values
        ddG[(from_lig, to_lig)] = decimal.Decimal(b_ddG)
        ddG[(to_lig, from_lig)] = -decimal.Decimal(b_ddG)

        # Initialize print_e and ddG_save with default values
        print_e[(from_lig, to_lig)] = decimal.Decimal(0)
        ddG_save[(from_lig, to_lig)] = decimal.Decimal(0)
        ddG_save[(to_lig, from_lig)] = decimal.Decimal(0)

        # Default error values
        err[(from_lig, to_lig)] = decimal.Decimal(-1)
        err[(to_lig, from_lig)] = decimal.Decimal(-1)

        # Default weight values
        weight[(from_lig, to_lig)].append(decimal.Decimal(1))
        weight[(to_lig, from_lig)].append(decimal.Decimal(1))

        # Update cumulative delta delta G (ddG_cc)
        ddG_cc[(from_lig, to_lig)].append(decimal.Decimal(b_ddG))
        ddG_cc[(to_lig, from_lig)].append(-decimal.Decimal(b_ddG))

        # Handle weights if provided
        if weights is not None:
            for w in weights:
                weight_value = decimal.Decimal(w)
                weight[(from_lig, to_lig)].append(weight_value**2)
                weight[(to_lig, from_lig)].append(weight_value**2)
                ddG_cc[(from_lig, to_lig)].append(decimal.Decimal(b_ddG))
                ddG_cc[(to_lig, from_lig)].append(-decimal.Decimal(b_ddG))
                err[(from_lig, to_lig)] = weight_value
                err[(to_lig, from_lig)] = weight_value

        # Update weight_num for both directions
        weight_num = len(weight[(from_lig, to_lig)])
        return weight_num
