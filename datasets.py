import random
from collections import defaultdict

def load_dimacs_gr(path, max_n=None):
    """DIMACS 9th Challenge (.gr): строки 'a u v w'. Дороги направленные, но почти все
    рёбра парные; берём как неориентированный граф с min весом."""
    return _load(path, lambda parts: (int(parts[1]), int(parts[2]), int(parts[3])) if parts[0] == 'a' else None, max_n)

def load_snap_edgelist(path, weight_range=(1, 10), seed=0, max_n=None):
    """SNAP roadNet-*.txt: 'u v' без весов, '#' — комментарии. Веса генерируем случайно."""
    rnd = random.Random(seed)
    return _load(path, lambda parts: (int(parts[0]), int(parts[1]), rnd.randint(*weight_range)), max_n)

def _load(path, parse, max_n):
    ids, W = {}, {}
    graph = defaultdict(list)
    def vid(x):
        if x not in ids: ids[x] = len(ids) + 1       # перенумерация в 1..N
        return ids[x]
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] in '#cp': continue
            r = parse(line.split())
            if r is None: continue
            a, b, w = vid(r[0]), vid(r[1]), r[2]
            if a == b: continue
            if (a, b) in W: W[(a, b)] = W[(b, a)] = min(W[(a, b)], w); continue
            graph[a].append(b); graph[b].append(a); W[(a, b)] = W[(b, a)] = w
    return len(ids), graph, W

def grid(rows, cols, wmax=10, seed=0, diag=0.0):
    """Дорожно-подобный граф: решётка + случайные веса (+ немного диагоналей)."""
    rnd = random.Random(seed); graph = defaultdict(list); W = {}
    vid = lambda r, c: r * cols + c + 1
    def add(a, b):
        w = rnd.randint(1, wmax); graph[a].append(b); graph[b].append(a); W[(a, b)] = W[(b, a)] = w
    for r in range(rows):
        for c in range(cols):
            if c + 1 < cols: add(vid(r, c), vid(r, c + 1))
            if r + 1 < rows: add(vid(r, c), vid(r + 1, c))
            if diag and r + 1 < rows and c + 1 < cols and rnd.random() < diag: add(vid(r, c), vid(r + 1, c + 1))
    return rows * cols, graph, W

def random_sparse(n, avg_deg=3, wmax=10, seed=0):
    """Случайный граф (экспандер) — «плохой» случай: деревья кратчайших путей мелкие и широкие."""
    rnd = random.Random(seed); graph = defaultdict(list); W = {}
    for v in range(2, n + 1):                     # остов для связности
        u = rnd.randint(1, v - 1); w = rnd.randint(1, wmax)
        graph[u].append(v); graph[v].append(u); W[(u, v)] = W[(v, u)] = w
    for _ in range(int(n * (avg_deg - 2) / 2)):
        a, b = rnd.randint(1, n), rnd.randint(1, n)
        if a != b and (a, b) not in W:
            w = rnd.randint(1, wmax); graph[a].append(b); graph[b].append(a); W[(a, b)] = W[(b, a)] = w
    return n, graph, W
