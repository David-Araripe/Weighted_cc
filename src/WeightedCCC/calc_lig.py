# -*- coding: utf-8 -*-
"""functions to calculate the ligand energies for a given set of nodes in a graph."""

import sys
import decimal
import numpy as np
import pandas as pd

decimal.getcontext().rounding = "ROUND_HALF_UP"


def set_node_map(graph):
    node_map = [
        [decimal.Decimal(sys.maxsize) if (val1 != val2) else decimal.Decimal(0.0) for val2 in range(len(graph.V))]
        for val1 in range(len(graph.V))
    ]
    for x, y, val in graph.nodelist:
        node_map[graph.V.index(x)][graph.V.index(y)] = node_map[graph.V.index(y)][graph.V.index(x)] = decimal.Decimal(
            val**2
        )
    return node_map


def cal_node_path_independent_error(node, node_map):
    node_rel = [[decimal.Decimal(-1) if (p == sys.maxsize) else p for p in n] for n in node_map]
    no_ref_error = [np.max(n) for n in node_rel]
    return no_ref_error


def cal_node_path_dependent_error(source, node, node_map):
    vertex_num = len(node)
    visited = [False for i in range(vertex_num)]
    visited[source] = True
    dist = [0 for i in range(vertex_num)]
    pre = [0 for i in range(vertex_num)]
    for i in range(len(node)):
        dist[i] = node_map[source][i]
        pre[i] = source
    for i in range(1, vertex_num):
        min_cost = sys.maxsize
        for j in range(vertex_num):
            if visited[j] is False and dist[j] < min_cost:
                min_cost_index = j
                min_cost = dist[j]
        visited[min_cost_index] = True
        for j in range(vertex_num):
            if (
                visited[j] is False
                and node_map[min_cost_index][j] != sys.maxsize
                and min_cost + node_map[min_cost_index][j] < dist[j]
            ):
                dist[j] = min_cost + node_map[min_cost_index][j]
                pre[j] = min_cost_index

    path_all = []
    for i in range(0, vertex_num):
        path = []
        path.append(i)
        if i != source:
            t = pre[i]
            path.append(t)
            while t != source:
                t = pre[t]
                path.append(t)
        path_all.append(path)

    return dist, path_all


def calcMolEnes(ref_ene, graph, path):
    mol_ene = [[ref_ene for i in range(len(graph.V))] for j in range(0, graph.weight_num)]
    for k in range(0, graph.weight_num):
        for m in range(len(graph.V)):
            route_mol_list = path[m]
            for i in range(len(route_mol_list) - 1):
                curr_path = (graph.V[route_mol_list[i]], graph.V[route_mol_list[i + 1]])
                curr_ene = graph.ddG_cc[curr_path][k]
                mol_ene[k][m] -= float(curr_ene)
    return mol_ene


def getMolEnergyDataFrame(nodes, mol_ene, path_dependent_error, path_independent_error, verbose=False):
    value_cols = [f"dG_wcc{k}" for k in range(0, len(mol_ene))] + ["path_dependent_error", "path_independent_error"]
    columns = ["Node", "dG_cc"] + value_cols
    data = []

    for i in range(len(nodes)):  # rows of values for each molecular pair
        row = [nodes[i], mol_ene[0][i]]
        row += [mol_ene[k][i] for k in range(0, len(mol_ene))]
        row.append(path_dependent_error[i].sqrt().quantize(decimal.Decimal("0.00")))
        row.append(path_independent_error[i])
        data.append(row)

    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    for col in value_cols:
        df[col] = df[col].astype(float)
    if verbose:
        print(df)
    return df
