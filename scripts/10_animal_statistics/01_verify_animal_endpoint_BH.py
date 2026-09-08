"""Verify public mouse-level two-group statistics and family-wise BH FDR.

Run from the repository root:
    python scripts/10_animal_statistics/01_verify_animal_endpoint_BH.py

This verifier uses only the Python standard library.
"""

import csv
import math
from collections import defaultdict
from pathlib import Path


INPUT = Path("data_public/figure_source_data/Animal_numerical_data/animal_endpoint_multiplicity_BH.csv")


def parse_values(text):
    return [float(x.strip()) for x in str(text).split(";") if x.strip()]


def beta_continued_fraction(a, b, x):
    max_iter, eps, fpmin = 200, 3.0e-14, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            return h
    raise ArithmeticError("Incomplete-beta continued fraction did not converge")


def regularized_incomplete_beta(x, a, b):
    if not 0.0 <= x <= 1.0:
        raise ValueError("x must be in [0,1]")
    if x in (0.0, 1.0):
        return x
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * beta_continued_fraction(a, b, x) / a
    return 1.0 - bt * beta_continued_fraction(b, a, 1.0 - x) / b


def equal_variance_t_pvalue(x, y):
    nx, ny = len(x), len(y)
    mx, my = sum(x) / nx, sum(y) / ny
    vx = sum((v - mx) ** 2 for v in x) / (nx - 1)
    vy = sum((v - my) ** 2 for v in y) / (ny - 1)
    df = nx + ny - 2
    pooled = ((nx - 1) * vx + (ny - 1) * vy) / df
    t = abs(mx - my) / math.sqrt(pooled * (1.0 / nx + 1.0 / ny))
    return regularized_incomplete_beta(df / (df + t * t), df / 2.0, 0.5)


def bh_adjust(p_values):
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    out = [0.0] * len(p_values)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * len(order) / rank)
        out[original_index] = min(running, 1.0)
    return out


def main():
    with INPUT.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 31:
        raise AssertionError(f"Expected 31 endpoint rows, found {len(rows)}")

    by_family = defaultdict(list)
    for index, row in enumerate(rows):
        hc = parse_values(row["Control_mouse_values"])
        dss = parse_values(row["DSS_mouse_values"])
        if len(hc) != int(row["n_HC"]) or len(dss) != int(row["n_DSS"]):
            raise AssertionError(f"Sample-size mismatch for {row['Endpoint']}")
        row["recomputed_raw_P"] = equal_variance_t_pvalue(hc, dss)
        by_family[row["Family"]].append(index)

    for indices in by_family.values():
        adjusted = bh_adjust([rows[i]["recomputed_raw_P"] for i in indices])
        for i, value in zip(indices, adjusted):
            rows[i]["recomputed_BH_P"] = value

    for row in rows:
        if not math.isclose(float(row["Raw_P"]), row["recomputed_raw_P"], rel_tol=2e-7, abs_tol=1e-12):
            raise AssertionError(f"Raw P mismatch for {row['Endpoint']}")
        if not math.isclose(float(row["BH_adjusted_P"]), row["recomputed_BH_P"], rel_tol=2e-7, abs_tol=1e-12):
            raise AssertionError(f"BH P mismatch for {row['Endpoint']}")
        if row["Significance_changed"].strip().lower() == "yes":
            raise AssertionError("Expected zero significance-status changes after BH correction")

    print("PASS: 31 mouse endpoints reproduced; BH changed significance for 0 endpoints.")


if __name__ == "__main__":
    main()
