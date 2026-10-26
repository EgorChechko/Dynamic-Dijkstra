"""Сравнение DynamicGraph с полным пересчётом Дейкстры на USA-road-d.NY.

1. На каждом шаге меняется источник (changeStart) и вес одного ребра.
2. Источник фиксирован, меняются только рёбра: где обновление быстрее полного пересчёта и во сколько раз.
"""
import random
import time
from pathlib import Path
from statistics import median

import datasets
from dynamic_graph import DynamicGraph, full_dijkstra

GR = Path(__file__).resolve().parent / "USA-road-d.NY.gr"
SRC = 1


def load():
    t0 = time.perf_counter()
    N, graph, W = datasets.load_dimacs_gr(GR)
    dt = time.perf_counter() - t0
    edges = sum(len(v) for v in graph.values()) // 2
    print(f"USA-road-d.NY: вершин {N}, рёбер {edges}, загрузка {dt:.1f} с")
    return N, graph, W


def check(g, ref, what):
    if g.d == ref:
        return
    bad = [v for v in range(1, g.N + 1) if g.d[v] != ref[v]][:8]
    print("РАСХОЖДЕНИЕ", what)
    for v in bad:
        print(f"  v={v} алгоритм={g.d[v]} дейкстра={ref[v]} родитель={g.P[v]}")
    raise SystemExit(1)


def timed_dijkstra(g):
    t0 = time.perf_counter()
    ref = full_dijkstra(g.N, g.graph, g.W, g.start)
    return ref, time.perf_counter() - t0


def scenario_change_start(N, graph, W, trials, rnd):
    g = DynamicGraph()
    g.init(N, graph, dict(W), SRC)
    edges = [e for e in g.W if e[0] < e[1]]
    print()
    print("1. Каждый шаг: changeStart и смена веса ребра")
    print(f"{'#':>3} {'источник':>10} {'ребро':>16} {'changeStart, мс':>16} {'ребро, мс':>10} {'дейкстра, мс':>13}")
    for i in range(1, trials + 1):
        src = rnd.randint(1, N)
        t0 = time.perf_counter()
        g.changeStart(src)
        t_start = time.perf_counter() - t0
        a, b = edges[rnd.randrange(len(edges))]
        old = g.W[(a, b)]
        new = old * 5 if rnd.randrange(2) else max(1, old // 2)
        if new == old:
            new = old * 5
        t0 = time.perf_counter()
        g.changeEdgeWeight(a, b, new)
        t_edge = time.perf_counter() - t0
        ref, t_full = timed_dijkstra(g)
        check(g, ref, f"changeStart {src}, {a}-{b}: {old}->{new}")
        edge = f"{a}-{b} {old}->{new}"
        print(f"{i:3} {src:10} {edge:>22} {t_start*1e3:16.1f} {t_edge*1e3:10.1f} {t_full*1e3:13.1f}")
    print("Расстояния совпали с полным пересчётом на каждом шаге.")
    print("changeStart сам заново запускает Дейкстру, поэтому вместе со сменой ребра шаг не быстрее одного полного пересчёта.")


def subtree_sizes(g):
    size = [0] * (g.N + 1)
    order = [v for v in range(1, g.N + 1) if g.d[v] < float("inf")]
    for v in order:
        size[v] = 1
    order.sort(key=lambda v: g.d[v], reverse=True)
    for v in order:
        p = g.P[v]
        if p != -1:
            size[p] += size[v]
    return size


def new_graph(N, graph, W):
    g = DynamicGraph()
    g.init(N, graph, dict(W), SRC)
    return g


def apply_edge(g, a, b, new, kind, rows):
    old = g.W[(a, b)]
    before = g.d[:]
    t0 = time.perf_counter()
    g.changeEdgeWeight(a, b, new)
    t_alg = time.perf_counter() - t0
    ref, t_full = timed_dijkstra(g)
    check(g, ref, f"{kind} {a}-{b}: {old}->{new}")
    changed = sum(x != y for x, y in zip(before[1:], g.d[1:]))
    speed = t_full / t_alg if t_alg > 0 else float("inf")
    rows.append({"kind": kind, "changed": changed, "speed": speed, "t_alg": t_alg, "t_full": t_full})
    print(f"{kind:28} {a:>7}-{b:<7} {changed:11} {t_alg*1e3:12.1f} {t_full*1e3:12.1f} {speed:10.1f}x")


def tree_edges_in_band(g, size, lo, hi):
    found = []
    for v in range(1, g.N + 1):
        p = g.P[v]
        if p == -1 or not (lo <= size[v] <= hi):
            continue
        if g.W.get((p, v), float("inf")) < float("inf"):
            found.append((p, v))
    return found


def random_new_weight(old, rnd):
    new = old * 5 if rnd.randrange(2) else max(1, old // 2)
    if new == old:
        new = old * 5
    return new


def scenario_edges(N, graph, W, rnd):
    lo, hi = int(0.10 * N), int(0.40 * N)
    g = new_graph(N, graph, W)
    print()
    print(f"Случайные рёбра, источник {SRC}, fallback = {g.fallback}")
    print(f"{'что меняли':28} {'ребро':>16} {'изменено':>11} {'алгоритм, мс':>13} {'дейкстра, мс':>13} {'выигрыш':>11}")
    edges = [e for e in g.W if e[0] < e[1] and g.W[e] < float("inf")]
    rows = []
    for _ in range(12):
        a, b = edges[rnd.randrange(len(edges))]
        old = g.W[(a, b)]
        if old == float("inf"):
            continue
        apply_edge(g, a, b, random_new_weight(old, rnd), "случайное", rows)

    probe = new_graph(N, graph, W)
    band = tree_edges_in_band(probe, subtree_sizes(probe), lo, hi)
    rnd.shuffle(band)
    print(f"Древесных рёбер с поддеревом {lo}..{hi}: {len(band)}")
    band_rows = []
    for a, b in band[:6]:
        h = new_graph(N, graph, W)
        old = h.W[(a, b)]
        apply_edge(h, a, b, old * 5, "случайное 10-40%", band_rows)
        rows.append(band_rows[-1])
    print("Расстояния совпали с полным пересчётом на каждом шаге.")
    print()
    print("Выигрыш — во сколько раз смена ребра быстрее полного пересчёта Дейкстры.")
    in_band = [r for r in rows if lo <= r["changed"] <= hi]
    if in_band:
        speeds = [r["speed"] for r in in_band]
        times = [r["t_alg"] for r in in_band]
        print(
            f"Изменилось от {lo} до {hi} расстояний: {len(in_band)} шагов, "
            f"медиана {median(r['changed'] for r in in_band):.0f} вершин, "
            f"медиана времени {median(times)*1e3:.0f} мс против полной Дейкстры "
            f"{median(r['t_full'] for r in in_band)*1e3:.0f} мс, "
            f"медиана выигрыша {median(speeds):.1f}x "
            f"(от {min(speeds):.1f}x до {max(speeds):.1f}x)."
        )
    else:
        print(f"Ни один шаг не изменил от {lo} до {hi} расстояний.")


def main():
    rnd = random.Random(1)
    N, graph, W = load()
    scenario_edges(N, graph, W, rnd)


if __name__ == "__main__":
    main()
