import heapq
from math import inf


def _relink(P, ch, u, new_parent):
    old = P[u]
    if old == new_parent:
        return
    if old != -1:
        kids = ch[old]
        for i, c in enumerate(kids):
            if c == u:
                last = kids.pop()
                if i < len(kids):
                    kids[i] = last
                break
    P[u] = new_parent
    if new_parent != -1:
        ch[new_parent].append(u)


class DynamicGraph:
    def __init__(self):
        self.W = None
        self.N = None
        self.start = None
        self.graph = None
        self.adj = None
        self.P = None
        self.childrens = None
        self.d = None
        self.fallback = 1

    def init(self, N, graph, W, start):
        self.W = W
        self.N = N
        self.start = start
        self.graph = graph
        self.adj = [[] for _ in range(N + 1)]
        for v in graph:
            for u in graph[v]:
                self.adj[v].append([u, W[(v, u)]])
        self.P = [-1] * (N + 1)
        self.childrens = [[] for _ in range(N + 1)]
        self.d = [inf] * (N + 1)
        self.d[start] = 0
        self._Dijkstra()
        self._buildChildrensByP()

    def _set_weight(self, x, y, newW):
        self.W[(x, y)] = newW
        self.W[(y, x)] = newW
        self._assign(x, y, newW)
        self._assign(y, x, newW)

    def _assign(self, x, y, newW):
        for edge in self.adj[x]:
            if edge[0] == y:
                edge[1] = newW
                return
        self.adj[x].append([y, newW])

    def changeStart(self, start):
        if self.graph is None:
            raise RuntimeError("call init() first")
        self.init(self.N, self.graph, self.W, start)

    def _Dijkstra(self):
        d = self.d
        P = self.P
        adj = self.adj
        heappush = heapq.heappush
        heappop = heapq.heappop
        seen = [0] * (self.N + 1)
        h = [(0, self.start)]
        while h:
            _, v = heappop(h)
            if seen[v]:
                continue
            seen[v] = 1
            dv = d[v]
            for u, w in adj[v]:
                nd = dv + w
                if d[u] > nd:
                    d[u] = nd
                    P[u] = v
                    heappush(h, (nd, u))

    def _buildChildrensByP(self):
        P = self.P
        ch = self.childrens
        for v in range(1, self.N + 1):
            p = P[v]
            if p != -1:
                ch[p].append(v)

    def changeEdgeWeight(self, x, y, newW):
        if self.W is None:
            raise RuntimeError("call init() first")
        old = self.W.get((x, y), self.W.get((y, x), inf))
        if newW < old:
            self._incrementalDijkstra(x, y, newW)
        elif newW > old:
            self._decrementalDijkstra(x, y, newW)

    def _incrementalDijkstra(self, x, y, newW):
        self._set_weight(x, y, newW)
        d = self.d
        P = self.P
        ch = self.childrens
        adj = self.adj
        if d[x] > d[y]:
            x, y = y, x
        if d[y] <= d[x] + newW:
            return
        d[y] = d[x] + newW
        _relink(P, ch, y, x)
        heappush = heapq.heappush
        heappop = heapq.heappop
        h = [(d[y], y)]
        while h:
            dist_v, v = heappop(h)
            if d[v] != dist_v:
                continue
            dv = d[v]
            for u, w in adj[v]:
                nd = dv + w
                if d[u] > nd:
                    d[u] = nd
                    _relink(P, ch, u, v)
                    heappush(h, (nd, u))

    def _decrementalDijkstra(self, x, y, newW):
        # Alshammari, Rezgui, AKCE 2020, Algorithm 2.
        self._set_weight(x, y, newW)
        P = self.P
        ch = self.childrens
        d = self.d
        adj = self.adj
        N = self.N
        start = self.start
        if P[y] != x:
            if P[x] == y:
                x, y = y, x
            else:
                return

        red = bytearray(N + 1)
        affected = 0
        stack = [y]
        while stack:
            w = stack.pop()
            if red[w]:
                continue
            red[w] = 1
            affected += 1
            stack.extend(ch[w])
        if affected > self.fallback * N:
            self.recompute()
            return

        ver = [0] * (N + 1)
        bestkey = [0.0] * (N + 1)
        seen_gen = [0] * (N + 1)
        seq = 0
        gen = 0
        h = []
        heappush = heapq.heappush
        heappop = heapq.heappop

        def insert(v, key):
            nonlocal seq
            if v == start:
                return
            cur = ver[v]
            if cur and key >= bestkey[v]:
                return
            seq += 1
            ver[v] = seq
            bestkey[v] = key
            heappush(h, (key, seq, v))

        def whiten(root):
            nonlocal gen
            gen += 1
            stack = [root]
            while stack:
                v = stack.pop()
                if seen_gen[v] == gen:
                    continue
                seen_gen[v] = gen
                red[v] = 0
                stack.extend(ch[v])

        insert(y, d[y])
        steps = 0
        limit = 3 * N
        while h:
            key, mark, u = heappop(h)
            if mark != ver[u]:
                continue
            ver[u] = 0
            steps += 1
            if steps > limit:
                self.recompute()
                return
            old = d[u]
            d[u] = inf
            _relink(P, ch, u, -1)

            best_v = -1
            best_d = inf
            best_is_red = 1
            for p, w in adj[u]:
                dp = d[p]
                if dp == inf or w == inf:
                    continue
                cand = dp + w
                p_red = red[p]
                if cand < best_d or (cand == best_d and best_is_red and not p_red):
                    best_d = cand
                    best_v = p
                    best_is_red = p_red

            if best_v != -1:
                if best_is_red:
                    d[u] = inf
                    insert(u, inf)
                else:
                    red[u] = 0
                    d[u] = best_d
                    _relink(P, ch, u, best_v)

            if d[u] != old:
                du = d[u]
                for v, w in adj[u]:
                    if P[v] == u:
                        insert(v, d[v])
                    elif w != inf and du + w < d[v]:
                        insert(v, du + w)
            else:
                whiten(u)

    def recompute(self):
        self.init(self.N, self.graph, self.W, self.start)

    def _set_parent(self, u, new_parent):
        _relink(self.P, self.childrens, u, new_parent)


def full_dijkstra(N, graph, W, start):
    d = [inf] * (N + 1)
    d[start] = 0
    heappush = heapq.heappush
    heappop = heapq.heappop
    h = [(0, start)]
    while h:
        dist_v, v = heappop(h)
        if d[v] != dist_v:
            continue
        nbrs = graph[v]
        for u in nbrs:
            nd = dist_v + W[(v, u)]
            if d[u] > nd:
                d[u] = nd
                heappush(h, (nd, u))
    return d
