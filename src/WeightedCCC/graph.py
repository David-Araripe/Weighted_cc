# -*- coding: utf-8 -*-
""""module with GraphClosure class, used to apply the weighted CCC algorithm"""

import math
import copy
import decimal
import pandas as pd

from .defaults import DefaultGraph

decimal.getcontext().rounding = "ROUND_HALF_UP"


class GraphClosure(DefaultGraph):
    def __init__(self, from_lig=None, to_lig=None, b_ddG=None, weights=None, filename=None):
        super().__init__(from_lig, to_lig, b_ddG, weights, filename)

    def getAllPathsUtil(self, u, d, visited, path, iscycle=False):
        self.num_iteration += 1
        if iscycle:
            if self.num_iteration > 1:
                visited[u] = True
        else:
            visited[u] = True
            if self.num_iteration == 1:
                self.start = u
        path.append(u)
        if self.num_iteration > 1 and u == d:
            if iscycle:
                if len(path) > 3:
                    path_s = copy.deepcopy(path)
                    path_s.sort()
                    path_str = {"".join(path_s)}
                    if not path_str < self.cycleset:
                        self.cycles.append(copy.deepcopy(path))
                        self.cycleset.add("".join(path_s))
            else:
                self.paths[(self.start, d)].append(copy.deepcopy(path))

        else:
            for i in self.graph[u]:
                if visited[i] is False:
                    self.getAllPathsUtil(i, d, visited, path, iscycle=iscycle)
        path.pop()
        visited[u] = False

    def getAllCyles(self):
        visit = dict(zip(self.graph.keys(), [False] * (len(self.graph))))
        for mol in self.V:
            visited = visit
            path = []
            self.num_iteration = 0
            self.getAllPathsUtil(mol, mol, visited, path, iscycle=True)
            visit[mol] = True
        return

    def getDelta(self, n, cycle_list):
        delta = decimal.Decimal(0.0)
        edges = 0
        std = decimal.Decimal(0)
        for i in range(len(cycle_list) - 1):
            mol1, mol2 = cycle_list[i : i + 2]
            delta += self.ddG_cc[mol1, mol2][n]
            edges += 1
            std += self.weight[mol1, mol2][n]
        return delta, edges, std

    def CycleClosure(self, n, edge_error):
        # cycle closure for all the cycles of a single molecule
        for k in range(len(self.cycles)):
            cycle_list = self.cycles[k]
            delta, edges, std_sum = self.getDelta(n, cycle_list)
            single_err = abs(delta / decimal.Decimal(math.sqrt(edges)))
            for i in range(len(cycle_list) - 1):
                mol1, mol2 = cycle_list[i : i + 2]
                if edge_error is True:
                    if edges > 6:  # ignore cycles more than 6edges
                        continue
                    if single_err > self.err[mol1, mol2]:
                        self.err[mol1, mol2] = single_err
                        self.err[mol2, mol1] = single_err
                    continue
                scale = self.weight[mol1, mol2][n] / std_sum
                ene = self.ddG_cc[mol1, mol2][n]
                newene = ene - scale * delta
                self.ddG_cc[mol1, mol2][n] = newene
                self.ddG_cc[mol2, mol1][n] = -newene

    def chk_continue(self, n, tol=0.001):
        for curr_key in self.ddG.keys():
            if abs(self.ddG_save[curr_key] - self.ddG_cc[curr_key][n]) > tol:
                return True
        return False

    def iterateCycleClosure(self, minimum_cycles=2):
        for n in range(0, self.weight_num):
            i = 0
            while i < minimum_cycles or self.chk_continue(n, 0.001):
                cal_error = True if (i == 0) else False  # if the first iteration, calculate pair error
                for curr_key in self.ddG.keys():  # save the current energy value for the next step
                    self.ddG_save[curr_key] = self.ddG_cc[curr_key][n]
                self.CycleClosure(n, cal_error)
                i += 1
        for molpair in self.print_e:
            self.nodelist.append([molpair[0], molpair[1], self.err[molpair]])

    def getEnergyPairsDataFrame(self, verbose=False):
        value_cols = [f"ddG_wcc{k}" for k in range(0, self.weight_num)] + ["pair_error"]
        columns = ["Pair"] + value_cols
        data = []

        for molpair in self.print_e:  # rows of values for each molecular pair
            row = [f"{molpair[0]}-{molpair[1]}"]
            row += [self.ddG_cc[molpair][k] for k in range(0, self.weight_num)]
            row.append(decimal.Decimal(self.err[molpair]).quantize(decimal.Decimal("0.00")))
            data.append(row)

        df = pd.DataFrame(data, columns=columns)
        for col in value_cols:
            df[col] = df[col].astype(float)  # convert back to float dtypes
        if verbose:
            print(df)
        return df
