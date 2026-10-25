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
        self._set_weight(x, y, newW)
        P = self.P
        ch = self.childrens
        d = self.d
        adj = self.adj
        N = self.N
        if P[y] != x:
            if P[x] == y:
                x, y = y, x
            else:
                return
        mark = [0] * (N + 1)
        affected = []
        heappush = heapq.heappush
        heappop = heapq.heappop
        limit = self.fallback * N
        h = [(d[y], y)]
        while h:
            _, u = heappop(h)
            du = d[u]
            alt = -1
            for p, w in adj[u]:
                if not mark[p] and d[p] + w == du:
                    alt = p
                    break
            if alt != -1:
                _relink(P, ch, u, alt)
                continue
            mark[u] = 1
            affected.append(u)
            if len(affected) > limit:
                self.recompute()
                return
            for c in ch[u]:
                heappush(h, (d[c], c))
        for v in affected:
            d[v] = inf
            _relink(P, ch, v, -1)
        h = []
        for v in affected:
            best = inf
            best_u = -1
            for u, w in adj[v]:
                if not mark[u]:
                    c = d[u] + w
                    if c < best:
                        best = c
                        best_u = u
            if best < inf:
                d[v] = best
                _relink(P, ch, v, best_u)
                heappush(h, (best, v))
        while h:
            dv, v = heappop(h)
            if dv != d[v]:
                continue
            for u, w in adj[v]:
                if mark[u]:
                    nd = dv + w
                    if nd < d[u]:
                        d[u] = nd
                        _relink(P, ch, u, v)
                        heappush(h, (nd, u))

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
