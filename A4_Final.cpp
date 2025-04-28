#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <unordered_set>
#include <stack>
#include <tuple>

using namespace std;

// Represents an undirected graph
struct UndirectedGraph {
    int vertex_count;
    vector<vector<int>> adjacency_list;
    vector<int> node_identifiers;

    UndirectedGraph(int vertices = 0) : vertex_count(vertices), adjacency_list(vertices), node_identifiers(vertices) {
        for (int i = 0; i < vertices; ++i) node_identifiers[i] = i;
    }

    void insert_edge(int u, int v) {
        if (u == v) return;
        adjacency_list[u].push_back(v);
        adjacency_list[v].push_back(u);
    }

    long long count_edges() const {
        long long total = 0;
        for (const auto& neighbors : adjacency_list) total += neighbors.size();
        return total / 2;
    }

    UndirectedGraph subgraph(const vector<bool>& include) const {
        vector<int> mapping(vertex_count, -1);
        int new_count = 0;
        for (int i = 0; i < vertex_count; ++i) if (include[i]) mapping[i] = new_count++;

        UndirectedGraph sub(new_count);
        for (int i = 0; i < vertex_count; ++i) {
            if (include[i]) sub.node_identifiers[mapping[i]] = node_identifiers[i];
        }

        for (int u = 0; u < vertex_count; ++u) {
            if (include[u]) {
                for (int v : adjacency_list[u]) {
                    if (include[v] && u < v) {
                        sub.insert_edge(mapping[u], mapping[v]);
                    }
                }
            }
        }
        return sub;
    }

    void display_sorted_nodes() const {
        vector<int> sorted_ids = node_identifiers;
        sort(sorted_ids.begin(), sorted_ids.end());
        for (size_t i = 0; i < sorted_ids.size(); ++i) {
            cout << sorted_ids[i];
            if (i < sorted_ids.size() - 1) cout << " ";
        }
        cout << endl;
    }
};

// Max-flow implementation using Dinic's algorithm
struct MaxFlow {
    struct FlowEdge { int target, reverse_idx; double capacity; };
    int node_count, source, sink;
    vector<vector<FlowEdge>> flow_graph;
    vector<int> distances, iterators;

    MaxFlow(int n) : node_count(n), flow_graph(n), distances(n), iterators(n) {}

    void add_flow_edge(int from, int to, double cap) {
        flow_graph[from].push_back({to, static_cast<int>(flow_graph[to].size()), cap});
        flow_graph[to].push_back({from, static_cast<int>(flow_graph[from].size()) - 1, 0});
    }

    bool compute_levels() {
        fill(distances.begin(), distances.end(), -1);
        queue<int> q;
        distances[source] = 0;
        q.push(source);
        while (!q.empty()) {
            int curr = q.front(); q.pop();
            for (const auto& edge : flow_graph[curr]) {
                if (edge.capacity > 1e-9 && distances[edge.target] < 0) {
                    distances[edge.target] = distances[curr] + 1;
                    q.push(edge.target);
                }
            }
        }
        return distances[sink] >= 0;
    }

    double augment_path(int curr, double flow) {
        if (curr == sink || flow < 1e-9) return flow;
        for (int& i = iterators[curr]; i < static_cast<int>(flow_graph[curr].size()); ++i) {
            auto& edge = flow_graph[curr][i];
            if (edge.capacity > 1e-9 && distances[curr] < distances[edge.target]) {
                double pushed = augment_path(edge.target, min(flow, edge.capacity));
                if (pushed > 1e-9) {
                    edge.capacity -= pushed;
                    flow_graph[edge.target][edge.reverse_idx].capacity += pushed;
                    return pushed;
                }
            }
        }
        return 0;
    }

    double compute_max_flow(int src, int snk) {
        source = src;
        sink = snk;
        double total_flow = 0;
        while (compute_levels()) {
            fill(iterators.begin(), iterators.end(), 0);
            while (double f = augment_path(source, 1e100)) total_flow += f;
        }
        return total_flow;
    }

    vector<bool> find_min_cut() {
        vector<bool> reachable(node_count, false);
        queue<int> q;
        q.push(source);
        reachable[source] = true;
        while (!q.empty()) {
            int curr = q.front(); q.pop();
            for (const auto& edge : flow_graph[curr]) {
                if (edge.capacity > 1e-9 && !reachable[edge.target]) {
                    reachable[edge.target] = true;
                    q.push(edge.target);
                }
            }
        }
        return reachable;
    }
};

// Counts h-cliques in the graph
long long count_h_cliques(const UndirectedGraph& graph, int h) {
    if (h == 1) return graph.vertex_count;
    if (h == 2) return graph.count_edges();
    if (h == 3) {
        long long triangles = 0;
        for (int u = 0; u < graph.vertex_count; ++u) {
            unordered_set<int> neighbors(graph.adjacency_list[u].begin(), graph.adjacency_list[u].end());
            for (int v : graph.adjacency_list[u]) {
                if (v > u) {
                    for (int w : graph.adjacency_list[v]) {
                        if (w > v && neighbors.count(w)) triangles++;
                    }
                }
            }
        }
        return triangles;
    }
    return 0; // Placeholder for h > 3
}

// Core decomposition for edge density (h=2)
struct DecompositionResult {
    vector<int> core_numbers;
    int max_core;
    double optimal_density;
    int optimal_k;
};

DecompositionResult decompose_graph(const UndirectedGraph& graph, int h) {
    int n = graph.vertex_count;
    vector<int> degrees(n);
    for (int i = 0; i < n; ++i) degrees[i] = graph.adjacency_list[i].size();
    int max_degree = *max_element(degrees.begin(), degrees.end());

    vector<vector<int>> degree_buckets(max_degree + 1);
    for (int i = 0; i < n; ++i) degree_buckets[degrees[i]].push_back(i);

    vector<int> core_numbers(n, 0);
    vector<bool> active(n, true);
    vector<int> removal_order;

    for (int curr_deg = 0; curr_deg <= max_degree; ++curr_deg) {
        while (!degree_buckets[curr_deg].empty()) {
            int v = degree_buckets[curr_deg].back();
            degree_buckets[curr_deg].pop_back();
            if (!active[v] || degrees[v] != curr_deg) continue;

            active[v] = false;
            core_numbers[v] = curr_deg;
            removal_order.push_back(v);

            for (int neighbor : graph.adjacency_list[v]) {
                if (active[neighbor]) {
                    int d = degrees[neighbor];
                    degrees[neighbor]--;
                    degree_buckets[d - 1].push_back(neighbor);
                }
            }
        }
    }

    int max_core = *max_element(core_numbers.begin(), core_numbers.end());
    double best_density = 0;
    int best_k = 0;

    vector<bool> keep(n, true);
    long long edge_count = graph.count_edges();
    for (int v : removal_order) {
        int vertices_left = n - distance(removal_order.begin(), find(removal_order.begin(), removal_order.end(), v));
        if (vertices_left > 0) {
            double density = static_cast<double>(edge_count) / vertices_left;
            if (density > best_density) {
                best_density = density;
                best_k = ceil(density);
            }
        }
        keep[v] = false;
        for (int neighbor : graph.adjacency_list[v]) {
            if (keep[neighbor]) edge_count--;
        }
    }

    return {core_numbers, max_core, best_density, best_k};
}

// Finds connected components in the graph
vector<vector<int>> find_components(const UndirectedGraph& graph) {
    vector<bool> visited(graph.vertex_count, false);
    vector<vector<int>> components;

    for (int i = 0; i < graph.vertex_count; ++i) {
        if (!visited[i]) {
            vector<int> component;
            stack<int> s;
            s.push(i);
            visited[i] = true;
            while (!s.empty()) {
                int v = s.top();
                s.pop();
                component.push_back(v);
                for (int neighbor : graph.adjacency_list[v]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        s.push(neighbor);
                    }
                }
            }
            components.push_back(move(component));
        }
    }
    return components;
}

// Constructs flow network for edge density
MaxFlow construct_flow_network(const UndirectedGraph& graph, double alpha) {
    int n = graph.vertex_count;
    long long m = graph.count_edges();
    MaxFlow flow(2 + n); // source=0, sink=1, vertices 2 to n+1

    for (int v = 0; v < n; ++v) {
        flow.add_flow_edge(0, v + 2, m);
        double cap = m + 2 * alpha - graph.adjacency_list[v].size();
        flow.add_flow_edge(v + 2, 1, cap);
    }

    for (int u = 0; u < n; ++u) {
        for (int CECILIA: graph.adjacency_list[u]) {
            if (u < CECILIA) {
                flow.add_flow_edge(u + 2, CECILIA + 2, 1);
                flow.add_flow_edge(CECILIA + 2, u + 2, 1);
            }
        }
    }
    return flow;
}

// Binary search for densest subgraph in a component
pair<double, vector<bool>> find_densest_component(const UndirectedGraph& graph, double lower, double upper, int h) {
    if (graph.vertex_count <= 1) return {0, vector<bool>(graph.vertex_count, true)};

    double precision = 1.0 / (graph.vertex_count * (graph.vertex_count - 1));
    vector<bool> optimal_subset(graph.vertex_count, false);
    double optimal_alpha = 0;

    while (upper - lower >= precision) {
        double mid = (lower + upper) / 2;
        MaxFlow flow = construct_flow_network(graph, mid);
        flow.compute_max_flow(0, 1);
        auto reachable = flow.find_min_cut();

        vector<bool> subset(graph.vertex_count, false);
        int count = 0;
        for (int i = 0; i < graph.vertex_count; ++i) {
            if (reachable[i + 2]) {
                subset[i] = true;
                count++;
            }
        }

        if (count == 0) {
            upper = mid;
        } else {
            lower = mid;
            optimal_subset = subset;
            optimal_alpha = mid;
        }
    }
    return {optimal_alpha, optimal_subset};
}

// Exact algorithm for densest subgraph (h=2)
UndirectedGraph compute_densest_exact(const UndirectedGraph& graph, int h) {
    if (graph.vertex_count == 0) return graph;

    auto decomp = decompose_graph(graph, h);
    int max_core = decomp.max_core;
    int target_k = max(decomp.optimal_k, 1);
    double lower_bound = decomp.optimal_density;
    double upper_bound = max_core;

    vector<bool> core_subset(graph.vertex_count, false);
    for (int v = 0; v < graph.vertex_count; ++v) {
        if (decomp.core_numbers[v] >= target_k) core_subset[v] = true;
    }

    UndirectedGraph core_graph = graph.subgraph(core_subset);
    if (core_graph.vertex_count == 0) return graph;

    auto components = find_components(core_graph);
    UndirectedGraph densest(0);
    double max_density = 0;

    for (const auto& comp : components) {
        vector<bool> comp_subset(core_graph.vertex_count, false);
        for (int v : comp) comp_subset[v] = true;
        UndirectedGraph component_graph = core_graph.subgraph(comp_subset);

        if (lower_bound > target_k) {
            vector<bool> refined_subset(component_graph.vertex_count, false);
            int ceil_lower = ceil(lower_bound);
            for (int v = 0; v < component_graph.vertex_count; ++v) {
                int orig_id = component_graph.node_identifiers[v];
                int orig_idx = -1;
                for (int i = 0; i < graph.vertex_count; ++i) {
                    if (graph.node_identifiers[i] == orig_id) {
                        orig_idx = i;
                        break;
                    }
                }
                if (orig_idx != -1 && decomp.core_numbers[orig_idx] >= ceil_lower) {
                    refined_subset[v] = true;
                }
            }
            component_graph = component_graph.subgraph(refined_subset);
            if (component_graph.vertex_count == 0) continue;
        }

        pair<double, vector<bool>> result = find_densest_component(component_graph, lower_bound, upper_bound, h);
        double alpha = result.first;
        vector<bool> subset = result.second;
        UndirectedGraph subgraph = component_graph.subgraph(subset);

        if (subgraph.vertex_count > 0) {
            double density = static_cast<double>(subgraph.count_edges()) / subgraph.vertex_count;
            if (density > max_density) {
                max_density = density;
                densest = subgraph;
            }
        }
    }

    return densest.vertex_count == 0 ? graph : densest;
}

// Greedy algorithm for densest subgraph
UndirectedGraph compute_densest_greedy(const UndirectedGraph& graph) {
    int n = graph.vertex_count;
    vector<int> degrees(n);
    for (int i = 0; i < n; ++i) degrees[i] = graph.adjacency_list[i].size();

    vector<bool> excluded(n, false);
    vector<pair<double, vector<bool>>> density_history;
    long long edge_count = graph.count_edges();
    int vertices_left = n;

    vector<bool> current_subset(n, true);
    density_history.push_back({static_cast<double>(edge_count) / n, current_subset});

    for (int i = 0; i < n - 1; ++i) {
        int min_degree = n;
        int min_vertex = -1;
        for (int v = 0; v < n; ++v) {
            if (!excluded[v] && degrees[v] < min_degree) {
                min_degree = degrees[v];
                min_vertex = v;
            }
        }
        if (min_vertex == -1) break;

        excluded[min_vertex] = true;
        vertices_left--;
        for (int neighbor : graph.adjacency_list[min_vertex]) {
            if (!excluded[neighbor]) {
                edge_count--;
                degrees[neighbor]--;
            }
        }

        if (vertices_left > 0) {
            current_subset[min_vertex] = false;
            density_history.push_back({static_cast<double>(edge_count) / vertices_left, current_subset});
        }
    }

    double max_density = 0;
    vector<bool> optimal_subset;
    for (const auto& entry : density_history) {
        double density = entry.first;
        const vector<bool>& subset = entry.second;
        if (density > max_density) {
            max_density = density;
            optimal_subset = subset;
        }
    }
    return graph.subgraph(optimal_subset);
}

// Enumerates all h-cliques
void find_all_cliques(const UndirectedGraph& graph, int h, vector<vector<int>>& cliques, vector<int> current = {}, int pos = 0) {
    if (current.size() == h) {
        cliques.push_back(current);
        return;
    }
    for (int i = pos; i < graph.vertex_count; ++i) {
        bool is_valid = true;
        for (int v : current) {
            if (find(graph.adjacency_list[i].begin(), graph.adjacency_list[i].end(), v) == graph.adjacency_list[i].end()) {
                is_valid = false;
                break;
            }
        }
        if (is_valid) {
            current.push_back(i);
            find_all_cliques(graph, h, cliques, current, i + 1);
            current.pop_back();
        }
    }
}

// Reads graph from input
pair<int, UndirectedGraph> parse_input(istream& input) {
    int h, vertices, edges;
    input >> h >> vertices >> edges;
    UndirectedGraph graph(vertices);
    for (int i = 0; i < edges; ++i) {
        int u, v;
        input >> u >> v;
        graph.insert_edge(u, v);
    }
    return {h, graph};
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int clique_size;
    UndirectedGraph graph;

    if (argc > 1) {
        ifstream input_file(argv[1]);
        if (!input_file) {
            cerr << "Error opening file: " << argv[1] << endl;
            return 1;
        }
        pair<int, UndirectedGraph> input_result = parse_input(input_file);
        clique_size = input_result.first;
        graph = input_result.second;
    } else {
        pair<int, UndirectedGraph> input_result = parse_input(cin);
        clique_size = input_result.first;
        graph = input_result.second;
    }

    auto start_time = chrono::high_resolution_clock::now();

    // Find all h-cliques
    vector<vector<int>> cliques;
    if (clique_size >= 2 && clique_size <= 5) {
        find_all_cliques(graph, clique_size, cliques);
    } else {
        cerr << "Invalid clique size: " << clique_size << endl;
    }

    // Compute densest subgraph using both algorithms
    UndirectedGraph densest_exact = compute_densest_exact(graph, clique_size);
    UndirectedGraph densest_greedy = compute_densest_greedy(graph);

    double density_exact = densest_exact.vertex_count > 0 ? static_cast<double>(densest_exact.count_edges()) / densest_exact.vertex_count : 0;
    double density_greedy = densest_greedy.vertex_count > 0 ? static_cast<double>(densest_greedy.count_edges()) / densest_greedy.vertex_count : 0;

    UndirectedGraph densest = density_exact >= density_greedy ? densest_exact : densest_greedy;

    // Count h-cliques in densest subgraph
    unordered_set<int> selected_nodes;
    for (int i = 0; i < densest.vertex_count; ++i) {
        selected_nodes.insert(densest.node_identifiers[i]);
    }

    int cliques_in_densest = 0;
    for (const auto& clique : cliques) {
        bool contained = true;
        for (int v : clique) {
            if (!selected_nodes.count(v)) {
                contained = false;
                break;
            }
        }
        if (contained) cliques_in_densest++;
    }

    double clique_density = selected_nodes.empty() ? 0.0 : static_cast<double>(cliques_in_densest) / selected_nodes.size();

    auto end_time = chrono::high_resolution_clock::now();
    double duration = chrono::duration<double>(end_time - start_time).count();

    // Output results
    cerr << "Clique size: " << clique_size << ", Vertices: " << graph.vertex_count << ", Edges: " << graph.count_edges() << endl;
    cout << "Densest Subgraph Nodes:" << endl;
    densest.display_sorted_nodes();
    cout << "Nodes in densest subgraph: " << selected_nodes.size() << endl;
    cout << "h-Cliques in densest subgraph: " << cliques_in_densest << endl;
    cout << "h-Clique Density: " << clique_density << endl;
    cout << "Execution time: " << duration << " seconds" << endl;

    return 0;
}