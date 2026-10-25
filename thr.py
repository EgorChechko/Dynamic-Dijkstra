import time, random, statistics as st
from dynamic_graph import DynamicGraph, full_dijkstra
import datasets

def run(name, N, G, W, src, factor=5, rand=120):
    rnd = random.Random(3)
    base = DynamicGraph()
    base.init(N, G, dict(W), src)
    ups = list(base.childrens[src])
    ups += [c for k in list(base.childrens[src]) for c in base.childrens[k]]
    ups += [c for k in ups[:] for c in base.childrens[k]][:40]       # глубина 3
    ups += [rnd.randint(2, N) for _ in range(rand)]
    ups = [v for v in ups if base.P[v] != -1]
    tf = []
    t0 = time.perf_counter(); full_dijkstra(N, G, W, src); tfull = time.perf_counter() - t0
    print(f"\n{name}: N={N}, обновлений={len(ups)}, полный пересчёт ≈ {tfull*1e3:.0f} мс")
    print(f"{'порог':>7} | {'итого, с':>8} | {'худший (dec/full)':>18} | {'выигрыш на всей серии':>22}")
    for fb in (1.0, 0.5, 0.3, 0.15, 0.05):
        g = DynamicGraph(); g.init(N, G, dict(W), src); g.fallback = fb; Wc = dict(W)
        tot = 0; worst = 0
        for v in ups:
            a, b = g.P[v], v
            if a == -1: continue
            old = Wc[(a, b)]; new = old * factor
            Wc[(a, b)] = Wc[(b, a)] = new
            t0 = time.perf_counter(); g._decrementalDijkstra(a, b, new); dt = time.perf_counter() - t0
            tot += dt; worst = max(worst, dt / tfull)
            Wc[(a, b)] = Wc[(b, a)] = old; g._incrementalDijkstra(a, b, old)
        label = "выкл." if fb == 1.0 else f"{fb:g}"
        print(f"{label:>7} | {tot:8.2f} | {worst:>17.2f}x | {len(ups)*tfull/tot:>21.1f}x")

N, G, W = datasets.grid(150, 150, wmax=10, seed=1)
run("решётка 150x150, источник в центре", N, G, W, 75*150+76)
N, G, W = datasets.random_sparse(20000, avg_deg=3, seed=3)
run("случайный граф 20000", N, G, W, 1)
