#include <iostream>
#include <vector>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <cmath>
#include <fstream>
#include <functional>
#include <chrono>
#include <set>  

using namespace std;
using namespace chrono;

struct AlgorithmOutput {
    vector<int> densestGraph;
    double graphDensity;
    int vertexCount;
    int edgeCount;
    int selfLoops;
    int iterationCount;
    double timeTaken;
};

struct FlowEdge {
    int start, end;
    double capacity, currentFlow;
    size_t reverseIndex;
};

vector<vector<FlowEdge>> flowNetwork;
vector<int> nodeLevels;
vector<size_t> edgePointers;

struct SetEquality {
    bool operator()(const unordered_set<int>& firstSet, const unordered_set<int>& secondSet) const {
        if (firstSet.size() != secondSet.size()) {
            return false;
        }
        return all_of(firstSet.begin(), firstSet.end(),
                    [&secondSet](int val) { return secondSet.count(val) > 0; });
    }
};

void initializeNetwork(int vertices) {
    flowNetwork.clear();
    flowNetwork.resize(vertices);
    nodeLevels.resize(vertices);
    edgePointers.resize(vertices);
}

void addEdgeToNetwork(int start, int end, double capacity) {
    size_t fromIndex = flowNetwork[start].size();
    size_t toIndex = flowNetwork[end].size();
    flowNetwork[start].emplace_back(FlowEdge{start, end, capacity, 0, toIndex});
    flowNetwork[end].emplace_back(FlowEdge{end, start, 0, 0, fromIndex});
}

bool buildLevelGraph(int source, int sink) {
    fill(nodeLevels.begin(), nodeLevels.end(), -1);
    nodeLevels[source] = 0;
    
    deque<int> queue = {source};
    
    while (!queue.empty() && nodeLevels[sink] == -1) {
        int currentNode = queue.front();
        queue.pop_front();
        
        for (const auto& edge : flowNetwork[currentNode]) {
            if (nodeLevels[edge.end] == -1 && edge.capacity > edge.currentFlow) {
                nodeLevels[edge.end] = nodeLevels[currentNode] + 1;
                queue.push_back(edge.end);
            }
        }
    }
    
    return nodeLevels[sink] != -1;
}

double pushFlow(int currentNode, int sink, double flowCap) {
    if (currentNode == sink) return flowCap;
    
    for (; edgePointers[currentNode] < flowNetwork[currentNode].size(); ++edgePointers[currentNode]) {
        auto& edge = flowNetwork[currentNode][edgePointers[currentNode]];
        
        if (nodeLevels[edge.end] == nodeLevels[currentNode] + 1 && edge.capacity > edge.currentFlow) {
            double pushed = pushFlow(edge.end, sink, min(flowCap, edge.capacity - edge.currentFlow));
                
            if (pushed > 0) {
                edge.currentFlow += pushed;
                flowNetwork[edge.end][edge.reverseIndex].currentFlow -= pushed;
                return pushed;
            }
        }
    }
    
    return 0;
}

struct CustomHasher {
    template <typename T>
    void combineHash(size_t& seed, const T& value) const {
        seed ^= hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    size_t operator()(const unordered_set<int>& set) const {
        size_t seed = 0;
        for (const auto& element : set) {
            combineHash(seed, element);
        }
        return seed;
    }
};

double calculateMaxFlow(int source, int sink) {
    double totalFlow = 0;
    
    while (buildLevelGraph(source, sink)) {
        fill(edgePointers.begin(), edgePointers.end(), 0);
        
        while (double pushed = pushFlow(source, sink, numeric_limits<double>::infinity())) {
            totalFlow += pushed;
        }
    }
    
    return totalFlow;
}

vector<bool> getMinimumCutSet(int source) {
    vector<bool> reachable(flowNetwork.size(), false);
    vector<bool> visited(flowNetwork.size(), false);
    deque<int> queue = {source};
    
    reachable[source] = visited[source] = true;
    
    while (!queue.empty()) {
        int currentNode = queue.front();
        queue.pop_front();
        
        for (const auto& edge : flowNetwork[currentNode]) {
            if (!visited[edge.end] && edge.capacity > edge.currentFlow) {
                reachable[edge.end] = visited[edge.end] = true;
                queue.push_back(edge.end);
            }
        }
    }
    
    return reachable;
}

vector<unordered_set<int>> findCliques(const vector<vector<int>>& adjList, int cliqueSize) {
    vector<unordered_set<int>> allCliques;
    vector<bool> inClique(adjList.size(), false);
    
    function<void(int, vector<int>&)> buildClique = [&](int start, vector<int>& currentClique) {
        if (currentClique.size() == cliqueSize - 1) {
            allCliques.emplace_back(currentClique.begin(), currentClique.end());
            return;
        }
        
        for (int v = start; v < adjList.size(); v++) {
            if (inClique[v]) continue;
            
            bool canAdd = true;
            for (int u : currentClique) {
                if (find(adjList[v].begin(), adjList[v].end(), u) == adjList[v].end()) {
                    canAdd = false;
                    break;
                }
            }
            
            if (canAdd) {
                currentClique.push_back(v);
                inClique[v] = true;
                buildClique(v + 1, currentClique);
                inClique[v] = false;
                currentClique.pop_back();
            }
        }
    };
    
    vector<int> currentClique;
    for (int i = 0; i < adjList.size(); i++) {
        currentClique.push_back(i);
        inClique[i] = true;
        buildClique(i + 1, currentClique);
        inClique[i] = false;
        currentClique.pop_back();
    }
    
    return allCliques;
}

struct Graph {
    int numVertices;  
    vector<vector<int>> adjacencyList; 
    
    Graph(int vertices) : numVertices(vertices), adjacencyList(vertices) {}
};

Graph readGraph(const string& filename) {
    ifstream file(filename);
    set<pair<int, int>> edges;
    int maxVertex = -1;
    int u, v;
    while (file >> u >> v) {
        if (u == v) continue;
        if (u > v) swap(u, v);
        maxVertex = max(maxVertex, max(u, v));
        edges.insert({u, v});
    }
    
    Graph graph(maxVertex + 1);
    for (const auto& edge : edges) {
        graph.adjacencyList[edge.first].push_back(edge.second);
        graph.adjacencyList[edge.second].push_back(edge.first);
    }
    
    return graph;
}

void displayGraphInfo(const Graph& graph) {
    cout << "Graph loaded successfully:" << endl;
    cout << "Number of vertices: " << graph.numVertices << endl;
    
    int edgeCount = 0;
    for (const auto& edges : graph.adjacencyList) {
        edgeCount += edges.size();
    }
    cout << "Number of edges: " << edgeCount / 2 << endl << endl;
}

int getCliqueSize(int maxVertices) {
    int h;
    cout << "Enter the value of h (clique size): ";
    while (!(cin >> h) || h <= 0 || h > maxVertices) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    return h;
}

vector<unordered_set<int>> processCliques(const vector<vector<int>>& graph, int h, vector<int>& degree) {
    vector<unordered_set<int>> cliques = findCliques(graph, h);
    
    for (const auto& sigma : cliques) {
        if (sigma.empty()) continue;
        
        unordered_set<int> commonNeighbors;
        bool first = true;
        
        for (int v : sigma) {
            if (first) {
                for (int neighbor : graph[v]) {
                    if (sigma.count(neighbor) == 0) {
                        commonNeighbors.insert(neighbor);
                    }
                }
                first = false;
            } else {
                unordered_set<int> newCommon;
                for (int neighbor : graph[v]) {
                    if (sigma.count(neighbor) == 0 && 
                        commonNeighbors.count(neighbor) != 0) {
                        newCommon.insert(neighbor);
                    }
                }
                commonNeighbors = move(newCommon);
            }
        }
        
        for (int v : commonNeighbors) {
            degree[v]++;
        }
    }
    return cliques;
}

AlgorithmOutput findDensestSubgraph(const vector<vector<int>>& graph, const vector<unordered_set<int>>& cliques, 
                                   const vector<int>& degree, int h, int n) {
    auto start = high_resolution_clock::now();
    
    int maxDegree = *max_element(degree.begin(), degree.end());
    double lowerBound = 0.0, upperBound = static_cast<double>(maxDegree);
    double epsilon = (n > 1) ? 1.0 / (n * (n - 1)) : 1e-9;
    
    vector<int> densestSubgraph;
    int iterations = 0;
    
    while (upperBound - lowerBound > epsilon) {
        iterations++;
        double lambda = (lowerBound + upperBound) / 2.0;
        int numNodes = 3 + n + cliques.size() - 1; 
        
        int source = 0;
        int sink = 1;  
        int vertexOffset = 2;
        int cliqueOffset = vertexOffset + n;
        initializeNetwork(numNodes);
        
        for (int v = 0; v < n; v++) {
            if (degree[v] > 0) {
                addEdgeToNetwork(source, vertexOffset + v, degree[v]);
            }
        }
        
        for (int v = 0; v < n; ++v) {
            addEdgeToNetwork(vertexOffset + v, sink, h * lambda);
        }
        
        for (size_t i = 0; i < cliques.size(); ++i) {
            const auto& sigma = cliques[i];
            int cliqueNode = cliqueOffset + i;
            
            for (int v : sigma) {
                addEdgeToNetwork(cliqueNode, vertexOffset + v, numeric_limits<double>::infinity());
            }
            
            unordered_set<int> neighbors;
            bool first = true;
            
            for (int v : sigma) {
                if (first) {
                    for (int neighbor : graph[v]) {
                        if (sigma.count(neighbor) == 0) {
                            neighbors.insert(neighbor);
                        }
                    }
                    first = false;
                } else {
                    unordered_set<int> newCommon;
                    for (int neighbor : graph[v]) {
                        if (sigma.count(neighbor) == 0 && 
                            neighbors.count(neighbor) != 0) {
                            newCommon.insert(neighbor);
                        }
                    }
                    neighbors = move(newCommon);
                }
            }
            
            for (int v : neighbors) {
                addEdgeToNetwork(vertexOffset + v, cliqueNode, 1.0);
            }
        }
        
        double maxFlow = calculateMaxFlow(source, sink);
        vector<bool> minCutSet = getMinimumCutSet(source);
        
        bool onlySourceInS = true;
        for (int i = 1; i < numNodes; ++i) {
            if (minCutSet[i]) {
                onlySourceInS = false;
                break;
            }
        }
        
        if (onlySourceInS) {
            upperBound = lambda;
        } else {
            lowerBound = lambda;
            densestSubgraph.clear();
            for (int v = 0; v < n; ++v) {
                if (minCutSet[vertexOffset + v]) {
                    densestSubgraph.push_back(v);
                }
            }
        }
        
        if (abs(upperBound - lowerBound) < 1e-12) {
            break;
        }
    }
    
    int numEdges = 0;
    int selfLoops = 0;
    int numVertices = densestSubgraph.size();
    
    for (int i = 0; i < densestSubgraph.size(); ++i) {
        for (int j = i; j < densestSubgraph.size(); ++j) {
            bool isEdge = false;
            for (int neighbor : graph[densestSubgraph[i]]) {
                if (neighbor == densestSubgraph[j]) {
                    isEdge = true;
                    break;
                }
            }
            
            if (isEdge) {
                numEdges++;
                if (densestSubgraph[i] == densestSubgraph[j]) {
                    selfLoops++;
                }
            }
        }
    }
    
    auto end = high_resolution_clock::now();
    double elapsedTime = duration<double>(end - start).count();
    
    return {densestSubgraph, lowerBound, numVertices, numEdges, selfLoops, iterations, elapsedTime};
}

void displayResults(const AlgorithmOutput& results, int h) {
    cout << "\nDensest Subgraph: ";
    vector<int> orderedSubgraph = results.densestGraph;
    sort(orderedSubgraph.begin(), orderedSubgraph.end());
    
    for (int i = 0; i < orderedSubgraph.size(); ++i) {
        cout << orderedSubgraph[i];
        if (i < orderedSubgraph.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "Number of vertices: " << results.vertexCount << endl;
    cout << "Number of edges: " << results.edgeCount << endl;
    cout << (h - 1) << "-Clique Density: " << round(results.graphDensity * 10000) / 10000 << endl;
    cout << "Execution Time: " << results.timeTaken << " s" << endl;
}

int main() {
    Graph graphData = readGraph("CA-HepTh.txt");
    
    vector<int> vertexDegree(graphData.numVertices, 0);
    displayGraphInfo(graphData);
    
    int h = getCliqueSize(graphData.numVertices);
    
    vector<unordered_set<int>> cliques = processCliques(graphData.adjacencyList, h, vertexDegree);
    
    AlgorithmOutput output = findDensestSubgraph(graphData.adjacencyList, cliques, vertexDegree, h, graphData.numVertices);
    displayResults(output, h);
    
    return 0;
}
