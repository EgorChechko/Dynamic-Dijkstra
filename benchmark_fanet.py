import argparse
import csv
import math
import time
from pathlib import Path

from dynamic_graph import DynamicGraph, full_dijkstra

ROOT = Path(__file__).resolve().parent
# DEFAULT_REPORT = ROOT / "the-one" / "reports" / "fanet_split_LinkBudgetReport.txt"
# DEFAULT_CSV = ROOT / "results_fanet_split.csv"

# DEFAULT_REPORT = ROOT / "the-one" / "reports" / "fanet_rotate_LinkBudgetReport.txt"
# DEFAULT_CSV = ROOT / "results_fanet_rotate.csv"

# Эти две строки — журнал и csv, которые берутся, если путь не передан
# аргументом. Для второго отчёта путь передаётся командой, см. ЗАПУСК.md.
# DEFAULT_REPORT = ROOT / "the-one" / "reports" / "fanet_halve_LinkQualityReport.txt"
# DEFAULT_CSV = ROOT / "results_fanet_halve_quality.csv"

DEFAULT_REPORT = ROOT / "the-one" / "reports" / "fanet_halve_LinkBudgetReport.txt"
DEFAULT_CSV = ROOT / "results_fanet_halve.csv"


def parse_report(path):
    events = []
    n_hosts = None
    hosts = None

    with open(path, encoding="utf-8-sig") as f:
        for line_no, raw in enumerate(f, 1):
            line = raw.strip()

            if not line:
                continue

            if line.startswith("#"):
                parts = line.split()

                if len(parts) >= 3 and parts[1].lower() == "nrofhosts":
                    n_hosts = int(parts[2])

                elif len(parts) >= 3 and parts[1].lower() == "hosts":
                    hosts = [int(x) for x in parts[2:]]

                continue

            p = line.split()

            t = float(p[0])
            op = p[1].upper()
            u, v = int(p[2]), int(p[3])

            if op == "W":
                w = float(p[4])
                if not math.isfinite(w) or w < 0:
                    raise ValueError(f"Invalid weight at line {line_no}")
            elif op == "DEL":
                w = math.inf
            else:
                raise ValueError(f"Unknown operation at line {line_no}: {op}")

            events.append((t, op, u, v, w, line_no))

    if not events:
        raise ValueError("No events found in report")

    if hosts is None:
        if n_hosts is not None:
            # The ONE host IDs in this report are 0..24.
            hosts = list(range(n_hosts))
        else:
            hosts = sorted({
                x for e in events for x in (e[2], e[3])
            })

    observed = {x for e in events for x in (e[2], e[3])}
    if not observed.issubset(set(hosts)):
        raise ValueError("Some event IDs are missing from the host list")

    return hosts, events


def add_edge(graph, weights, u, v, w):
    if v not in graph[u]:
        graph[u].append(v)

    if u not in graph[v]:
        graph[v].append(u)

    weights[(u, v)] = w
    weights[(v, u)] = w


def add_to_instance_graph(instance, u, v):
    # recompute() rebuilds adj from graph, so graph must be synchronized.
    if v not in instance.graph[u]:
        instance.graph[u].append(v)

    if u not in instance.graph[v]:
        instance.graph[v].append(u)


def compare_distances(a, b, atol=1e-8, rtol=1e-8):
    max_error = 0.0

    for x, y in zip(a[1:], b[1:]):
        if math.isinf(x) and math.isinf(y):
            continue

        if not math.isfinite(x) or not math.isfinite(y):
            return False, math.inf

        error = abs(x - y)
        max_error = max(max_error, error)

        if error > atol + rtol * max(abs(x), abs(y)):
            return False, max_error

    return True, max_error


def scenario_name(report):
    stem = report.stem
    for suffix in ("_LinkBudgetReport", "_LinkQualityReport"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "report",
        nargs="?",
        type=Path,
        default=DEFAULT_REPORT,
        help="LinkBudgetReport or LinkQualityReport from The ONE",
    )
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    if not args.report.is_file():
        raise SystemExit(f"missing report: {args.report}")

    print(f"Report: {args.report.resolve()}")
    hosts, events = parse_report(args.report)

    n = len(hosts)

    # External The ONE ID -> internal DynamicGraph ID (1..N).
    id_map = {
        external_id: i + 1
        for i, external_id in enumerate(hosts)
    }

    # First timestamp is the initial graph snapshot.
    first_time = min(e[0] for e in events)

    initial = [e for e in events if e[0] == first_time]
    updates = [e for e in events if e[0] != first_time]

    graph = {v: [] for v in range(1, n + 1)}
    weights = {}

    for t, op, ext_u, ext_v, w, line_no in initial:
        u = id_map[ext_u]
        v = id_map[ext_v]
        add_edge(graph, weights, u, v, w)

    # One independent DynamicGraph per source.
    dynamic = []

    for source in range(1, n + 1):
        g = {v: list(neighbors) for v, neighbors in graph.items()}
        W = dict(weights)

        instance = DynamicGraph()
        instance.init(n, g, W, source)

        dynamic.append(instance)

    # Independent graph for the full Dijkstra baseline.
    baseline_graph = {
        v: list(neighbors)
        for v, neighbors in graph.items()
    }
    baseline_weights = dict(weights)

    columns = [
        "event",
        "time",
        "operation",
        "u",
        "v",
        "weight",
        "dynamic_ms",
        "dijkstra_ms",
        "speedup",
        "correct",
        "max_error",
    ]

    rows = []
    failed = False

    for index, event in enumerate(updates, 1):
        t, op, ext_u, ext_v, w, line_no = event

        u = id_map[ext_u]
        v = id_map[ext_v]

        # Dynamic algorithm.
        started = time.perf_counter()

        for instance in dynamic:
            add_to_instance_graph(instance, u, v)
            instance.changeEdgeWeight(u, v, w)

        dynamic_ms = (time.perf_counter() - started) * 1000

        # Full recomputation.
        add_edge(baseline_graph, baseline_weights, u, v, w)

        started = time.perf_counter()

        reference = [
            full_dijkstra(
                n,
                baseline_graph,
                baseline_weights,
                source,
            )[0]
            for source in range(1, n + 1)
        ]

        dijkstra_ms = (time.perf_counter() - started) * 1000

        # Verify every source.
        correct = True
        max_error = 0.0

        for instance, expected in zip(dynamic, reference):
            ok, error = compare_distances(instance.d, expected)

            correct = correct and ok
            max_error = max(max_error, error)

        failed |= not correct

        speedup = (
            dijkstra_ms / dynamic_ms
            if dynamic_ms > 0
            else math.inf
        )

        rows.append({
            "event": index,
            "time": t,
            "operation": op,
            "u": ext_u,
            "v": ext_v,
            "weight": "" if math.isinf(w) else w,
            "dynamic_ms": dynamic_ms,
            "dijkstra_ms": dijkstra_ms,
            "speedup": speedup,
            "correct": correct,
            "max_error": max_error,
        })

        print(
            f"[{index}/{len(updates)}] "
            f"t={t:.3f} {op} {ext_u}-{ext_v} | "
            f"Dynamic: {dynamic_ms:.3f} ms | "
            f"Dijkstra: {dijkstra_ms:.3f} ms | "
            f"Speedup: {speedup:.2f}x | "
            f"{'OK' if correct else 'ERROR'}"
        )

        if not correct:
            print("  Maximum distance error:", max_error)

    args.csv.parent.mkdir(parents=True, exist_ok=True)

    with open(args.csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    print("\n========== RESULTS ==========")
    print("Hosts:", n)
    print("Initial events:", len(initial))
    print("Update events:", len(updates))
    print("CSV:", args.csv.resolve())

    if rows:
        dynamic_total = sum(r["dynamic_ms"] for r in rows)
        dijkstra_total = sum(r["dijkstra_ms"] for r in rows)

        print(f"Dynamic total: {dynamic_total:.3f} ms")
        print(f"Dijkstra total: {dijkstra_total:.3f} ms")

        if dynamic_total > 0:
            print(
                f"Overall speedup: "
                f"{dijkstra_total / dynamic_total:.2f}x"
            )

    if failed:
        print("\nERROR: distance verification failed")
        status = 1
    else:
        print("\nAll distance checks passed")
        status = 0

    print("Scenario:", scenario_name(args.report))
    return status


if __name__ == "__main__":
    raise SystemExit(main())