#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <chrono>

class Graph {
public:
    Graph(int numVertices) : numVertices(numVertices) {
        adjacencyMatrix.resize(numVertices, std::vector<int>(numVertices, 0));
    }

    virtual ~Graph() {}

    void addEdge(int v1, int v2, int weight = 1) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
            adjacencyMatrix[v2][v1] = weight;
        }
    }

    virtual void addVertex() {
        numVertices++;
        for (int i = 0; i < numVertices; ++i) {
            adjacencyMatrix[i].push_back(0);
        }
        adjacencyMatrix.push_back(std::vector<int>(numVertices, 0));
    }

    void removeEdge(int v1, int v2) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = 0;
            adjacencyMatrix[v2][v1] = 0;
        }
    }

    virtual void removeVertex(int v) {
        if (v >= 0 && v < numVertices) {
            numVertices--;
            adjacencyMatrix.erase(adjacencyMatrix.begin() + v);

            for (int i = 0; i < numVertices; ++i) {
                adjacencyMatrix[i].erase(adjacencyMatrix[i].begin() + v);
            }
        }
    }

    void print() const {
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                std::cout << adjacencyMatrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void generateRandomGraph() {
        srand(static_cast<unsigned>(time(0)));

        for (int i = 0; i < numVertices; ++i) {
            for (int j = i + 1; j < numVertices; ++j) {
                if (rand() % 2 == 1) {
                    int weight = 1 + rand() % 10;
                    addEdge(i, j, weight);
                }
            }
        }
    }

    void convertToListGraph(std::vector<std::vector<int>>& adjacencyList) {
        adjacencyList.clear();
        adjacencyList.resize(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j] != 0) {
                    adjacencyList[i].push_back(j);
                }
            }
        }
    }

protected:
    int numVertices;
    std::vector<std::vector<int>> adjacencyMatrix;
};

class DirectedGraph : public Graph {
public:
    DirectedGraph(int numVertices) : Graph(numVertices) {}

    void addEdge(int v1, int v2, int weight = 1) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
        }
    }

    void removeEdge(int v1, int v2) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = 0;
        }
    }
};

class UndirectedGraph : public Graph {
public:
    UndirectedGraph(int numVertices) : Graph(numVertices) {}

    void removeVertex(int v) {
        if (v >= 0 && v < numVertices) {
            throw std::logic_error("Cannot remove vertex in an undirected graph");
        }
    }

    void addVertex() {
        throw std::logic_error("Cannot add vertex to an undirected graph");
    }
};

class WeightedGraph : public Graph {
public:
    WeightedGraph(int numVertices) : Graph(numVertices) {}

    void addEdge(int v1, int v2, int weight) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
        }
    }

    void generateRandomGraph() {
        Graph::generateRandomGraph();
    }

    class DijkstraAlgorithm {
    public:
        DijkstraAlgorithm(WeightedGraph& graph) : numVertices(graph.getNumVertices()), graph(graph) {}

        void shortestPath(int source, int destination) {
            std::vector<int> distance(numVertices, INT_MAX);
            distance[source] = 0;
            std::vector<int> visited(numVertices, 0);

            for (int i = 0; i < numVertices - 1; ++i) {
                int minDistance = INT_MAX;
                int minIndex = -1;
                for (int v = 0; v < numVertices; ++v) {
                    if (!visited[v] && distance[v] < minDistance) {
                        minDistance = distance[v];
                        minIndex = v;
                    }
                }

                visited[minIndex] = true;

                for (int v = 0; v < numVertices; ++v) {
                    if (!visited[v] && graph.getEdgeWeight(minIndex, v) != 0 && distance[minIndex] != INT_MAX &&
                        distance[minIndex] + graph.getEdgeWeight(minIndex, v) < distance[v]) {
                        distance[v] = distance[minIndex] + graph.getEdgeWeight(minIndex, v);
                    }
                }
            }

            // if (distance[destination] == INT_MAX) {
            //     std::cout << "No path from " << source << " to " << destination << std::endl;
            // }
            // else {
            //     std::cout << "Shortest path from " << source << " to " << destination << ": " << distance[destination] << std::endl;
            // }
        }

    private:
        int numVertices;
        WeightedGraph& graph;
    };

    class FloydWarshallAlgorithm {
    public:
        FloydWarshallAlgorithm(WeightedGraph& graph) : numVertices(graph.getNumVertices()), graph(graph) {}

        void shortestPaths() {
            std::vector<std::vector<int>>& adjacencyMatrix = graph.getAdjacencyMatrix();

            for (int k = 0; k < numVertices; ++k) {
                for (int i = 0; i < numVertices; ++i) {
                    for (int j = 0; j < numVertices; ++j) {
                        if (adjacencyMatrix[i][k] != 0 && adjacencyMatrix[k][j] != 0 &&
                            adjacencyMatrix[i][k] + adjacencyMatrix[k][j] < adjacencyMatrix[i][j]) {
                            adjacencyMatrix[i][j] = adjacencyMatrix[i][k] + adjacencyMatrix[k][j];
                        }
                    }
                }
            }

            // for (int i = 0; i < numVertices; ++i) {
              //   for (int j = 0; j < numVertices; ++j) {
                    // if (adjacencyMatrix[i][j] == INT_MAX) {
                    //     std::cout << "No path from " << i << " to " << j << std::endl;
                   //  }
                    // else {
                      //   std::cout << "Shortest path from " << i << " to " << j << ": " << adjacencyMatrix[i][j] << std::endl;
                    // }
               //  }
            // }
        }

    private:
        int numVertices;
        WeightedGraph& graph;
    };

    std::vector<std::vector<int>>& getAdjacencyMatrix() {
        return adjacencyMatrix;
    }

    int getNumVertices() const {
        return numVertices;
    }

    int getEdgeWeight(int v1, int v2) const {
        return adjacencyMatrix[v1][v2];
    }
};

int main() {
    std::cout << "DirectedGraph" << "\n";
    DirectedGraph directedGraph(3);
    directedGraph.generateRandomGraph();
    directedGraph.print();

    std::cout << "UndirectedGraph" << "\n";
    UndirectedGraph undirectedGraph(4);
    undirectedGraph.generateRandomGraph();
    undirectedGraph.print();

    std::cout << "WeightedGraph" << "\n";
    WeightedGraph weightedGraph(3);
    weightedGraph.generateRandomGraph();
    weightedGraph.print();

    std::cout << "DijkstraAlgorithm" << "\n";
    WeightedGraph::DijkstraAlgorithm dijkstraGraph(weightedGraph);
    dijkstraGraph.shortestPath(0, 2);

    std::cout << "FloydWarshallAlgorithm" << "\n";
    WeightedGraph::FloydWarshallAlgorithm floydWarshallGraph(weightedGraph);
    floydWarshallGraph.shortestPaths();

    std::vector<int> graphSizes = { 10, 20, 50, 100, 200, 300 };

    for (int size : graphSizes) {
        WeightedGraph weightedGraph(size);
        weightedGraph.generateRandomGraph();


        auto startDijkstra = std::chrono::high_resolution_clock::now();
        WeightedGraph::DijkstraAlgorithm dijkstraGraph(weightedGraph);
        dijkstraGraph.shortestPath(0, size - 1);
        auto endDijkstra = std::chrono::high_resolution_clock::now();
        auto durationDijkstra = std::chrono::duration_cast<std::chrono::microseconds>(endDijkstra - startDijkstra);

        auto startFloydWarshall = std::chrono::high_resolution_clock::now();
        WeightedGraph::FloydWarshallAlgorithm floydWarshallGraph(weightedGraph);
        floydWarshallGraph.shortestPaths();
        auto endFloydWarshall = std::chrono::high_resolution_clock::now();
        auto durationFloydWarshall = std::chrono::duration_cast<std::chrono::microseconds>(endFloydWarshall - startFloydWarshall);

        std::cout << "Graph Size: " << size << "\n";
        std::cout << "Dijkstra Algorithm Execution Time: " << durationDijkstra.count() << " microseconds\n";
        std::cout << "Floyd-Warshall Algorithm Execution Time: " << durationFloydWarshall.count() << " microseconds\n";
    }

    return 0;
}

