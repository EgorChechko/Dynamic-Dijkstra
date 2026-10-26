"""Бесконечная сверка DynamicGraph со полной Дейкстрой.

Запуск из каталога ramangalim:

    python stress_dijkstra.py

Каждый шаг меняет вес, добавляет ребро или источник. Расстояния и дерево
сравниваются с full_dijkstra. Первое расхождение печатается и процесс
завершается с кодом 1. Пока расхождений нет, цикл не останавливается.
"""

import random
import sys
from collections import defaultdict
from math import inf
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "ramangalim"))
from dynamic_graph import DynamicGraph, full_dijkstra


def make_graph(rnd, n, extra, wmax):
    graph = defaultdict(list)
    W = {}
    edges = []

    def add(a, b, w):
        if a == b or (a, b) in W:
            return False
        graph[a].append(b)
        graph[b].append(a)
        W[(a, b)] = W[(b, a)] = w
        edges.append((a, b) if a < b else (b, a))
        return True

    for v in range(2, n + 1):
        add(rnd.randint(1, v - 1), v, rnd.randint(1, wmax))
    guard = 0
    while len(edges) < n - 1 + extra and guard < extra * 30 + 20:
        guard += 1
        add(rnd.randint(1, n), rnd.randint(1, n), rnd.randint(1, wmax))
    return graph, W, edges


def fail(text):
    print(text, flush=True)
    raise SystemExit(1)


def check(g, head):
    ref = full_dijkstra(g.N, g.graph, g.W, g.start)
    bad = []
    for v in range(1, g.N + 1):
        if g.d[v] != ref[v]:
            bad.append(v)
    if bad:
        lines = [head, "расхождения расстояний:"]
        for v in bad[:12]:
            lines.append(
                f"  v={v} алгоритм={g.d[v]} дейкстра={ref[v]} родитель={g.P[v]}"
            )
        if len(bad) > 12:
            lines.append(f"  ... ещё {len(bad) - 12}")
        fail("\n".join(lines))

    parent_count = [0] * (g.N + 1)
    for v in range(1, g.N + 1):
        if v == g.start:
            if g.d[v] != 0 or g.P[v] != -1:
                fail(f"{head}\nисточник: d={g.d[v]} P={g.P[v]}")
            continue
        if g.d[v] == inf:
            if g.P[v] != -1:
                fail(f"{head}\nнедостижимая {v} хранит родителя {g.P[v]}")
            continue
        p = g.P[v]
        w = g.W.get((p, v))
        if p < 1 or w is None or g.d[v] != g.d[p] + w:
            fail(f"{head}\nребро дерева не тугое: {p}->{v} d={g.d[v]} d[p]={g.d[p]} w={w}")
        if v not in g.childrens[p]:
            fail(f"{head}\n{v} нет в childrens[{p}]")
        parent_count[p] += 1
        x = v
        for _ in range(g.N + 1):
            x = g.P[x]
            if x == g.start:
                break
            if x == -1:
                fail(f"{head}\nпуть из {v} оборвался")
        else:
            fail(f"{head}\nцикл в дереве из {v}")

    for p in range(1, g.N + 1):
        kids = g.childrens[p]
        if len(kids) != parent_count[p] or len(kids) != len(set(kids)):
            fail(f"{head}\nchildrens[{p}]={kids} не совпадает с родителями")
        for c in kids:
            if g.P[c] != p:
                fail(f"{head}\nchildrens[{p}] содержит {c}, родитель {g.P[c]}")


def one_update(rnd, g, edges, wmax):
    kind = rnd.randrange(10)
    if kind == 0:
        src = rnd.randint(1, g.N)
        g.changeStart(src)
        return f"changeStart({src})"
    if kind == 1 and len(edges) < g.N * (g.N - 1) // 2:
        a = rnd.randint(1, g.N)
        b = rnd.randint(1, g.N)
        if a != b and (a, b) not in g.W:
            w = rnd.randint(1, wmax*10)
            g.graph[a].append(b)
            g.graph[b].append(a)
            g.changeEdgeWeight(a, b, w)
            edges.append((a, b) if a < b else (b, a))
            return f"новое {a}-{b} вес {w}"
    if not edges:
        return "пусто"
    a, b = edges[rnd.randrange(len(edges))]
    old = g.W.get((a, b), inf)
    if rnd.randrange(8) == 0:
        new = inf
    else:
        new = rnd.randint(1, wmax * 3)
    g.changeEdgeWeight(a, b, new)
    return f"{a}-{b}: {old} -> {new}"


def main():
    master_seed = random.randrange(1 << 30)
    master = random.Random(master_seed)
    print(f"мастер-seed {master_seed}", flush=True)
    checks = 0
    graph_i = 0
    while True:
        graph_i += 1
        seed = master.randrange(1 << 30)
        rnd = random.Random(seed)
        n = rnd.randint(4, 36)
        extra = rnd.randint(0, n * 2)
        wmax = rnd.randint(2, 15)
        src = rnd.randint(1, n)
        graph, W, edges = make_graph(rnd, n, extra, wmax)
        g = DynamicGraph()
        g.fallback = 1
        g.init(n, graph, W, src)
        head = f"seed={seed} граф={graph_i} N={n} источник={src}"
        check(g, head + " после init")
        checks += 1
        ops_left = rnd.randint(40, 180)
        for step in range(1, ops_left + 1):
            op = one_update(rnd, g, edges, wmax)
            check(g, f"{head} шаг={step} {op}")
            checks += 1
        print(
            f"ok  графов={graph_i}  проверок={checks}  последний seed={seed} N={n}",
            flush=True,
        )


if __name__ == "__main__":
    main()
