#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <tuple>
#include <cctype>

using namespace std;

static constexpr double INF = numeric_limits<double>::infinity();

struct Edge {
    int to;
    double w;
};

struct Event {
    double time;
    string op;
    int u;
    int v;
    double w;
    int line;
};

// ============================================================
// DynamicGraph
// ============================================================

class DynamicGraph {
private:
    int N = 0;
    int start = -1;

    vector<vector<Edge>> adj;

    // Parent in the shortest-path tree.
    vector<int> parent;

    // Children in the shortest-path tree.
    vector<vector<int>> children;

    // Shortest distances from start.
    vector<double> d;

    // Fraction of affected vertices above which full recomputation
    // is used. With fallback = 1.0 this condition is effectively off.
    double fallback = 1.0;

    // --------------------------------------------------------
    // Relink vertex u to new_parent in the shortest-path tree.
    // --------------------------------------------------------

    void relink(int u, int new_parent) {
        int old_parent = parent[u];

        if (old_parent == new_parent) {
            return;
        }

        if (old_parent != -1) {
            auto& kids = children[old_parent];

            for (size_t i = 0; i < kids.size(); ++i) {
                if (kids[i] == u) {
                    kids[i] = kids.back();
                    kids.pop_back();
                    break;
                }
            }
        }

        parent[u] = new_parent;

        if (new_parent != -1) {
            children[new_parent].push_back(u);
        }
    }

    // --------------------------------------------------------
    // Add or update an undirected edge in adjacency lists.
    // --------------------------------------------------------

    void assignEdge(int u, int v, double w) {
        for (auto& edge : adj[u]) {
            if (edge.to == v) {
                edge.w = w;
                return;
            }
        }

        adj[u].push_back({v, w});
    }

    void setWeight(int u, int v, double w) {
        assignEdge(u, v, w);
        assignEdge(v, u, w);
    }

    // --------------------------------------------------------
    // Initial full Dijkstra.
    // --------------------------------------------------------

    void dijkstra() {
        using State = pair<double, int>;

        priority_queue<
            State,
            vector<State>,
            greater<State>
        > pq;

        vector<char> seen(N, false);

        d[start] = 0.0;
        pq.push({0.0, start});

        while (!pq.empty()) {
            auto [dist_v, v] = pq.top();
            pq.pop();

            if (seen[v]) {
                continue;
            }

            seen[v] = true;

            for (const auto& edge : adj[v]) {
                int u = edge.to;
                double w = edge.w;

                if (!isfinite(w)) {
                    continue;
                }

                double nd = dist_v + w;

                if (d[u] > nd) {
                    d[u] = nd;
                    parent[u] = v;
                    pq.push({nd, u});
                }
            }
        }
    }

    // --------------------------------------------------------
    // Build children lists from parent array.
    // --------------------------------------------------------

    void buildChildren() {
        for (int v = 0; v < N; ++v) {
            int p = parent[v];

            if (p != -1) {
                children[p].push_back(v);
            }
        }
    }

    // --------------------------------------------------------
    // Incremental update: edge weight decreased.
    // --------------------------------------------------------

    void incrementalDijkstra(int x, int y, double newW) {
        setWeight(x, y, newW);

        if (d[x] > d[y]) {
            swap(x, y);
        }

        // The new edge does not improve the current distance.
        if (d[y] <= d[x] + newW) {
            return;
        }

        d[y] = d[x] + newW;
        relink(y, x);

        using State = pair<double, int>;

        priority_queue<
            State,
            vector<State>,
            greater<State>
        > pq;

        pq.push({d[y], y});

        while (!pq.empty()) {
            auto [dist_v, v] = pq.top();
            pq.pop();

            if (d[v] != dist_v) {
                continue;
            }

            for (const auto& edge : adj[v]) {
                int u = edge.to;
                double w = edge.w;

                if (!isfinite(w)) {
                    continue;
                }

                double nd = dist_v + w;

                if (d[u] > nd) {
                    d[u] = nd;
                    relink(u, v);
                    pq.push({nd, u});
                }
            }
        }
    }

    // --------------------------------------------------------
    // Decremental update: edge weight increased or deleted.
    //
    // Based on Algorithm 2:
    // Alshammari, Rezgui, AKCE 2020.
    // --------------------------------------------------------

    void decrementalDijkstra(int x, int y, double newW) {
        setWeight(x, y, newW);

        // The changed edge must be an edge of the shortest-path tree.
        if (parent[y] != x) {
            if (parent[x] == y) {
                swap(x, y);
            } else {
                // The changed edge is not in the shortest-path tree.
                return;
            }
        }

        // Mark the subtree rooted at y as affected (red).
        vector<char> red(N, false);

        int affected = 0;
        vector<int> stack;
        stack.push_back(y);

        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();

            if (red[v]) {
                continue;
            }

            red[v] = true;
            ++affected;

            for (int child : children[v]) {
                stack.push_back(child);
            }
        }

        // Fallback to full recomputation if too many vertices are affected.
        if (affected > fallback * N) {
            recompute();
            return;
        }

        // Version counters allow lazy deletion from the priority queue.
        vector<uint64_t> ver(N, 0);
        vector<double> bestkey(N, 0.0);

        uint64_t seq = 0;

        using State = tuple<double, uint64_t, int>;

        priority_queue<
            State,
            vector<State>,
            greater<State>
        > pq;

        auto insert = [&](int v, double key) {
            if (v == start) {
                return;
            }

            uint64_t cur = ver[v];

            if (cur != 0 && key >= bestkey[v]) {
                return;
            }

            ++seq;
            ver[v] = seq;
            bestkey[v] = key;

            pq.push({key, seq, v});
        };

        // Mark a subtree as unaffected (white).
        vector<uint64_t> seen_gen(N, 0);
        uint64_t gen = 0;

        auto whiten = [&](int root) {
            ++gen;

            vector<int> st;
            st.push_back(root);

            while (!st.empty()) {
                int v = st.back();
                st.pop_back();

                if (seen_gen[v] == gen) {
                    continue;
                }

                seen_gen[v] = gen;
                red[v] = false;

                for (int child : children[v]) {
                    st.push_back(child);
                }
            }
        };

        insert(y, d[y]);

        int steps = 0;
        const int limit = 3 * N;

        while (!pq.empty()) {
            auto [key, mark, u] = pq.top();
            pq.pop();

            if (mark != ver[u]) {
                continue;
            }

            ver[u] = 0;

            ++steps;

            if (steps > limit) {
                recompute();
                return;
            }

            double old = d[u];

            // Temporarily remove u from the shortest-path tree.
            d[u] = INF;
            relink(u, -1);

            int best_v = -1;
            double best_d = INF;

            // true means the current best predecessor is red.
            bool best_is_red = true;

            // Find the best predecessor of u.
            for (const auto& edge : adj[u]) {
                int p = edge.to;
                double w = edge.w;

                if (!isfinite(d[p]) || !isfinite(w)) {
                    continue;
                }

                double candidate = d[p] + w;
                bool p_red = red[p];

                if (
                    candidate < best_d ||
                    (
                        candidate == best_d &&
                        best_is_red &&
                        !p_red
                    )
                ) {
                    best_d = candidate;
                    best_v = p;
                    best_is_red = p_red;
                }
            }

            if (best_v != -1) {
                if (best_is_red) {
                    // The best predecessor is still affected.
                    // Wait until it gets a correct distance.
                    d[u] = INF;
                    insert(u, INF);
                } else {
                    // The best predecessor is already unaffected.
                    red[u] = false;
                    d[u] = best_d;
                    relink(u, best_v);
                }
            }

            if (d[u] != old) {
                double du = d[u];

                // Update affected children and possible alternative paths.
                for (const auto& edge : adj[u]) {
                    int v = edge.to;
                    double w = edge.w;

                    if (parent[v] == u) {
                        insert(v, d[v]);
                    } else if (
                        isfinite(w) &&
                        du + w < d[v]
                    ) {
                        insert(v, du + w);
                    }
                }
            } else {
                // The distance did not change: this subtree can be whitened.
                whiten(u);
            }
        }
    }

public:
    DynamicGraph() = default;

    // --------------------------------------------------------
    // Initialization.
    //
    // graph[v] contains all neighbors of v with their weights.
    // Vertices are numbered from 0 to N-1.
    // --------------------------------------------------------

    void init(
        int n,
        const vector<vector<Edge>>& graph,
        int source
    ) {
        if (n <= 0) {
            throw invalid_argument("N must be positive");
        }

        if (source < 0 || source >= n) {
            throw invalid_argument("Invalid source vertex");
        }

        N = n;
        start = source;

        adj = graph;

        parent.assign(N, -1);
        children.assign(N, {});
        d.assign(N, INF);

        dijkstra();
        buildChildren();
    }

    // --------------------------------------------------------
    // Change edge weight.
    //
    // newW = INF means edge deletion.
    // --------------------------------------------------------

    void changeEdgeWeight(int x, int y, double newW) {
        if (N == 0) {
            throw runtime_error("Call init() first");
        }

        if (
            x < 0 || x >= N ||
            y < 0 || y >= N
        ) {
            throw out_of_range("Invalid vertex");
        }

        if (x == y) {
            throw invalid_argument("Self-loops are not supported");
        }

        if (isnan(newW) || newW < 0.0) {
            throw invalid_argument("Weight must be non-negative");
        }

        double oldW = INF;

        for (const auto& edge : adj[x]) {
            if (edge.to == y) {
                oldW = edge.w;
                break;
            }
        }

        if (newW == oldW) {
            return;
        }

        if (newW < oldW) {
            incrementalDijkstra(x, y, newW);
        } else {
            decrementalDijkstra(x, y, newW);
        }
    }

    // --------------------------------------------------------
    // Full recomputation using the current adjacency lists.
    // --------------------------------------------------------

    void recompute() {
        vector<vector<Edge>> current = adj;
        init(N, current, start);
    }

    void setFallback(double value) {
        if (value < 0.0 || value > 1.0) {
            throw invalid_argument("Fallback must be in [0, 1]");
        }

        fallback = value;
    }

    const vector<double>& distances() const {
        return d;
    }

    const vector<int>& parents() const {
        return parent;
    }

    const vector<vector<int>>& getChildren() const {
        return children;
    }

    int size() const {
        return N;
    }
};


// ============================================================
// Full Dijkstra: reference implementation for benchmarking.
// ============================================================

vector<double> fullDijkstra(
    int N,
    const vector<vector<Edge>>& graph,
    int start
) {
    vector<double> d(N, INF);
    vector<char> seen(N, false);

    using State = pair<double, int>;

    priority_queue<
        State,
        vector<State>,
        greater<State>
    > pq;

    d[start] = 0.0;
    pq.push({0.0, start});

    while (!pq.empty()) {
        auto [dist_v, v] = pq.top();
        pq.pop();

        if (seen[v]) {
            continue;
        }

        seen[v] = true;

        for (const auto& edge : graph[v]) {
            int u = edge.to;
            double w = edge.w;

            if (!isfinite(w)) {
                continue;
            }

            double nd = dist_v + w;

            if (d[u] > nd) {
                d[u] = nd;
                pq.push({nd, u});
            }
        }
    }

    return d;
}


// ============================================================
// Report parser
// ============================================================

struct ParsedReport {
    vector<int> hosts;
    vector<Event> events;
};

ParsedReport parseReport(const string& filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        throw runtime_error("Cannot open report: " + filename);
    }

    ParsedReport result;

    int nrofHosts = -1;
    vector<int> declaredHosts;

    string line;
    int lineNo = 0;

    while (getline(file, line)) {
        ++lineNo;

        if (line.empty()) {
            continue;
        }

        istringstream iss(line);

        if (line[0] == '#') {
            string hash;
            string key;

            iss >> hash >> key;

            if (key == "nrofHosts") {
                iss >> nrofHosts;
            } else if (key == "hosts") {
                int id;

                while (iss >> id) {
                    declaredHosts.push_back(id);
                }
            }

            continue;
        }

        double t;
        string op;
        int u, v;

        if (!(iss >> t >> op >> u >> v)) {
            throw runtime_error(
                "Invalid event at line " + to_string(lineNo)
            );
        }

        for (char& c : op) {
            c = static_cast<char>(toupper(c));
        }

        double w = INF;

        if (op == "W") {
            if (!(iss >> w)) {
                throw runtime_error(
                    "Missing weight at line " + to_string(lineNo)
                );
            }

            if (!isfinite(w) || w < 0.0) {
                throw runtime_error(
                    "Invalid weight at line " + to_string(lineNo)
                );
            }
        } else if (op == "DEL") {
            w = INF;
        } else {
            throw runtime_error(
                "Unknown operation at line " + to_string(lineNo) +
                ": " + op
            );
        }

        result.events.push_back({
            t, op, u, v, w, lineNo
        });
    }

    if (result.events.empty()) {
        throw runtime_error("No events found in report");
    }

    if (!declaredHosts.empty()) {
        result.hosts = declaredHosts;
    } else if (nrofHosts > 0) {
        // The ONE host IDs are assumed to be 0..nrofHosts-1.
        result.hosts.resize(nrofHosts);

        for (int i = 0; i < nrofHosts; ++i) {
            result.hosts[i] = i;
        }
    } else {
        // Infer host IDs from events.
        for (const auto& e : result.events) {
            result.hosts.push_back(e.u);
            result.hosts.push_back(e.v);
        }

        sort(result.hosts.begin(), result.hosts.end());

        result.hosts.erase(
            unique(result.hosts.begin(), result.hosts.end()),
            result.hosts.end()
        );
    }

    sort(result.hosts.begin(), result.hosts.end());

    result.hosts.erase(
        unique(result.hosts.begin(), result.hosts.end()),
        result.hosts.end()
    );

    for (const auto& e : result.events) {
        if (
            !binary_search(result.hosts.begin(), result.hosts.end(), e.u) ||
            !binary_search(result.hosts.begin(), result.hosts.end(), e.v)
        ) {
            throw runtime_error(
                "Event contains unknown host ID at line " +
                to_string(e.line)
            );
        }
    }

    return result;
}


// ============================================================
// Graph utilities
// ============================================================

void addOrUpdateEdge(
    vector<vector<Edge>>& graph,
    int u,
    int v,
    double w
) {
    auto update = [](vector<Edge>& edges, int to, double weight) {
        for (auto& edge : edges) {
            if (edge.to == to) {
                edge.w = weight;
                return;
            }
        }

        edges.push_back({to, weight});
    };

    update(graph[u], v, w);
    update(graph[v], u, w);
}


// ============================================================
// Distance comparison
// ============================================================

pair<bool, double> compareDistances(
    const vector<double>& a,
    const vector<double>& b,
    double atol = 1e-8,
    double rtol = 1e-8
) {
    double maxError = 0.0;

    if (a.size() != b.size()) {
        return {false, INF};
    }

    for (size_t i = 0; i < a.size(); ++i) {
        double x = a[i];
        double y = b[i];

        if (isinf(x) && isinf(y)) {
            continue;
        }

        if (!isfinite(x) || !isfinite(y)) {
            return {false, INF};
        }

        double error = abs(x - y);
        maxError = max(maxError, error);

        if (error > atol + rtol * max(abs(x), abs(y))) {
            return {false, maxError};
        }
    }

    return {true, maxError};
}


// ============================================================
// Benchmark
// ============================================================

int runBenchmark(
    const string& reportFilename,
    const string& csvFilename
) {
    ParsedReport report = parseReport(reportFilename);

    const int N = static_cast<int>(report.hosts.size());

    // External The ONE ID -> internal ID [0, N).
    unordered_map<int, int> idMap;

    for (int i = 0; i < N; ++i) {
        idMap[report.hosts[i]] = i;
    }

    // --------------------------------------------------------
    // The first timestamp is treated as the initial snapshot.
    // All events at that timestamp initialize the graph.
    // --------------------------------------------------------

    double firstTime = report.events.front().time;

    for (const auto& e : report.events) {
        firstTime = min(firstTime, e.time);
    }

    vector<Event> initialEvents;
    vector<Event> updates;

    for (const auto& e : report.events) {
        if (e.time == firstTime) {
            initialEvents.push_back(e);
        } else {
            updates.push_back(e);
        }
    }

    vector<vector<Edge>> graph(N);

    // Apply initial events in their original report order.
    for (const auto& e : initialEvents) {
        int u = idMap.at(e.u);
        int v = idMap.at(e.v);

        addOrUpdateEdge(graph, u, v, e.w);
    }

    // --------------------------------------------------------
    // One DynamicGraph instance per source.
    // --------------------------------------------------------

    vector<DynamicGraph> dynamic;
    dynamic.reserve(N);

    for (int source = 0; source < N; ++source) {
        DynamicGraph instance;
        instance.init(N, graph, source);

        dynamic.push_back(move(instance));
    }

    // Graph for the full Dijkstra baseline.
    vector<vector<Edge>> baselineGraph = graph;

    // --------------------------------------------------------
    // CSV.
    // --------------------------------------------------------

    ofstream csv(csvFilename);

    if (!csv.is_open()) {
        throw runtime_error("Cannot create CSV: " + csvFilename);
    }

    csv << "event,time,operation,u,v,weight,"
        << "dynamic_ms,dijkstra_ms,speedup,correct,max_error\n";

    csv << setprecision(15);

    bool failed = false;

    double dynamicTotal = 0.0;
    double dijkstraTotal = 0.0;

    // --------------------------------------------------------
    // Process updates.
    // --------------------------------------------------------

    for (size_t index = 0; index < updates.size(); ++index) {
        const Event& e = updates[index];

        int u = idMap.at(e.u);
        int v = idMap.at(e.v);

        // Dynamic algorithm.
        // Timing includes applying the update to every instance.
        auto dynamicStart = chrono::steady_clock::now();

        for (auto& instance : dynamic) {
            instance.changeEdgeWeight(u, v, e.w);
        }

        auto dynamicEnd = chrono::steady_clock::now();

        double dynamicMs =
            chrono::duration<double, milli>(
                dynamicEnd - dynamicStart
            ).count();

        // Full Dijkstra baseline.
        // Timing includes updating the baseline graph.
        auto dijkstraStart = chrono::steady_clock::now();

        addOrUpdateEdge(baselineGraph, u, v, e.w);

        vector<vector<double>> reference;
        reference.reserve(N);

        for (int source = 0; source < N; ++source) {
            reference.push_back(
                fullDijkstra(N, baselineGraph, source)
            );
        }

        auto dijkstraEnd = chrono::steady_clock::now();

        double dijkstraMs =
            chrono::duration<double, milli>(
                dijkstraEnd - dijkstraStart
            ).count();

        // ----------------------------------------------------
        // Verify all sources.
        // ----------------------------------------------------

        bool correct = true;
        double maxError = 0.0;

        for (int source = 0; source < N; ++source) {
            auto [ok, error] = compareDistances(
                dynamic[source].distances(),
                reference[source]
            );

            correct = correct && ok;
            maxError = max(maxError, error);
        }

        if (!correct) {
            failed = true;
        }

        double speedup =
            dynamicMs > 0.0
                ? dijkstraMs / dynamicMs
                : INF;

        dynamicTotal += dynamicMs;
        dijkstraTotal += dijkstraMs;

        // ----------------------------------------------------
        // Write CSV row.
        // ----------------------------------------------------

        csv
            << index + 1 << ','
            << e.time << ','
            << e.op << ','
            << e.u << ','
            << e.v << ',';

        if (isfinite(e.w)) {
            csv << e.w;
        }

        csv << ','
            << dynamicMs << ','
            << dijkstraMs << ','
            << speedup << ','
            << (correct ? 1 : 0) << ','
            << maxError << '\n';

        // ----------------------------------------------------
        // Console output.
        // ----------------------------------------------------

        cout
            << '[' << index + 1 << '/' << updates.size() << "] "
            << "t=" << fixed << setprecision(3) << e.time
            << " " << e.op << " "
            << e.u << "-" << e.v
            << " | Dynamic: "
            << setprecision(3) << dynamicMs << " ms"
            << " | Dijkstra: "
            << dijkstraMs << " ms"
            << " | Speedup: "
            << setprecision(2) << speedup << "x"
            << " | "
            << (correct ? "OK" : "ERROR")
            << '\n';

        if (!correct) {
            cout << "    Maximum distance error: "
                 << maxError << '\n';
        }
    }

    // --------------------------------------------------------
    // Summary.
    // --------------------------------------------------------

    cout << "\n========== RESULTS ==========\n";

    cout << "Hosts: " << N << '\n';
    cout << "Initial events: " << initialEvents.size() << '\n';
    cout << "Update events: " << updates.size() << '\n';

    cout << "CSV: " << csvFilename << '\n';

    cout << fixed << setprecision(3);

    cout << "Dynamic total: "
         << dynamicTotal << " ms\n";

    cout << "Dijkstra total: "
         << dijkstraTotal << " ms\n";

    if (dynamicTotal > 0.0) {
        cout << "Overall speedup: "
             << setprecision(2)
             << dijkstraTotal / dynamicTotal
             << "x\n";
    }

    if (failed) {
        cout << "\nERROR: distance verification failed\n";
        return 1;
    }

    cout << "\nAll distance checks passed\n";

    return 0;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[]) {
    string reportFilename =
        argc >= 2 ? argv[1] : ".\\fanet_LinkBudgetReport.txt";

    string csvFilename =
        argc >= 3 ? argv[2] : ".\\results.csv";

    try {
        return runBenchmark(reportFilename, csvFilename);
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}