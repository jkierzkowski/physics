#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <ctime>
#include <thread>
#include <cstdlib>
#include <cmath>
#include <set>
#include <utility>
#include <numeric>
#include <algorithm>
#include <map>
#include <list>
#include <type_traits>

using namespace std;
using Graph = map<pair<int, int>, double>;
using Matrix = map<int, vector<double>>;
using Neighbors = map<int, set<int>>;
using EdgeNeighbors = map<int, map<int, double>>;
using Timeline = map<int, Graph>;

using namespace std;

template <typename T>
void insert_value(vector<T> &container, const T &value)
{
    container.push_back(value);
}

template <typename T>
void insert_value(set<T> &container, const T &value)
{
    container.insert(value);
}

template <typename T, typename Container>
vector<Container> import_from_file(const string &filename)
{
    vector<Container> out;
    ifstream OpenFile(filename);
    string input;

    while (getline(OpenFile, input))
    {
        if (input.empty())
            continue;

        istringstream iss(input);

        if constexpr (is_same<Container, pair<T, T>>::value)
        {
            T a, b;
            iss >> a >> b;
            out.emplace_back(a, b);
        }
        else
        {
            Container row;
            T value;

            while (iss >> value)
            {
                insert_value(row, value);
            }

            out.push_back(row);
        }
    }

    OpenFile.close();
    return out;
}

auto import_d(const string &path)
{
    vector<double> out;
    ifstream OpenFile(path);
    string input;
    while (getline(OpenFile, input))
    {
        out.push_back(stod(input));
    }
    return out;
}

double length(const vector<double> &pos1, const vector<double> &pos2)
{
    return sqrt(pow(pos2[0] - pos1[0], 2) + pow(pos2[1] - pos1[1], 2));
}

pair<int, int> normalizePair(const int &a, const int &b)
{
    return (a <= b) ? make_pair(a, b) : make_pair(b, a);
}

EdgeNeighbors edge_neighbors_init(const vector<pair<int, int>> &edges, const vector<double> &diameters, const Neighbors &neighbors, const vector<vector<double>> &pos)
{
    EdgeNeighbors edge_neighbors;
    for (auto edge : edges)
    {
        set<int> commonNeighbors;
        int n0 = edge.first;
        int n1 = edge.second;
        auto itA = neighbors.find(n0);
        auto itB = neighbors.find(n1);
        if (itA == neighbors.end() || itB == neighbors.end())
        {
            cerr << "Edge nodes not in neighbor map\n";
            return {};
        }
        const auto &nbrA = itA->second;
        const auto &nbrB = itB->second;
        // finding nodes neigboring the edge
        set_intersection(
            nbrA.begin(), nbrA.end(),
            nbrB.begin(), nbrB.end(),
            inserter(commonNeighbors, commonNeighbors.begin()));
        for (auto n : commonNeighbors)
        {

            pair<int, int> e1 = normalizePair(n0, n);
            pair<int, int> e2 = normalizePair(n1, n);
            auto num1 = find(edges.begin(), edges.end(), e1);
            auto num2 = find(edges.begin(), edges.end(), e2);
            auto num = find(edges.begin(), edges.end(), edge);
            int e_indx;
            if (num != edges.end())
            {
                e_indx = distance(edges.begin(), num);
            }
            else
            {
                cerr << "Warning cannot find edge (" << edge.first << "," << edge.second
                     << ") in edges";
                break;
            }
            if (num1 != edges.end() && num2 != edges.end())
            {
                int index = distance(edges.begin(), num1);
                edge_neighbors[e_indx][index] = length(pos[n0], pos[n]);
                index = distance(edges.begin(), num2);
                edge_neighbors[e_indx][index] = length(pos[n1], pos[n]);
            }
            else
            {
                cerr << "ERROR problem with neighbors, cannot create edge_neighbors";
            }
        }
    }
    return edge_neighbors;
}

Neighbors convertToMap(const vector<set<int>> &neighbors_vec)
{
    Neighbors neighbors_map;
    for (size_t i = 0; i < neighbors_vec.size(); ++i)
    {
        neighbors_map[i] = neighbors_vec[i];
    }
    return neighbors_map;
}

void import(EdgeNeighbors &edge_neighbors, vector<double> &diameters, vector<pair<int, int>> &edges)
{
    Neighbors neighbors;
    vector<vector<double>> pos = import_from_file<double, vector<double>>("pos.txt");
    neighbors = convertToMap(import_from_file<int, set<int>>("neighbors.txt"));
    edges = import_from_file<int, pair<int, int>>("edges.txt");
    diameters = import_d("D.txt");

    if (edges.size() == 0 || pos.size() == 0 || diameters.size() == 0 || neighbors.size() == 0)
    {
        cerr << "Import error";
        return;
    }
    edge_neighbors = edge_neighbors_init(edges, diameters, neighbors, pos);
}

int addRandom(vector<double> &diameters, const vector<int> &deleted_edges, const vector<int> &merged_edges, const double &constant)
{
    set<int> excluded(deleted_edges.begin(), deleted_edges.end());
    excluded.insert(merged_edges.begin(), merged_edges.end());

    // Collect valid indices
    vector<int> valid_indices;
    for (int i = 0; i < static_cast<int>(diameters.size()); ++i)
    {
        if (excluded.count(i) == 0)
        {
            valid_indices.push_back(i);
        }
    }

    if (!valid_indices.empty())
    {
        int rand_index = valid_indices[rand() % valid_indices.size()];
        diameters[rand_index] *= 1+constant;
        return rand_index;
    }
    cout << "No valid indices to modify.\n";
    return -1;
}

int common_neighbor(map<int, double> N1, map<int, double> N2)
{
    if(N2.size()==0)return -1;
    for (const auto &[key, _] : N1)
    {
        if (N2.count(key))
        {
            return key;
        }
    }
    cerr << "ERROR cannot find common edge neighbor\n";
    return -1;
}
void edge_neighbors_update(EdgeNeighbors &edge_neighbors, int id, int id_del, int id_merge)
{

    double d = edge_neighbors[id][id_merge];
    // adding new neighbors to bolded edge
    for (auto &N : edge_neighbors[id_merge])
    {
        int n_neighhbor = N.first;
        double l = N.second;
        if (!edge_neighbors[id].count(n_neighhbor) && n_neighhbor != id && n_neighhbor != id_del)
        {
            edge_neighbors[id][n_neighhbor] = l + d;
        }
    }
    // deleting the edge that has been merged and the egde that has been deleted
    edge_neighbors.erase(id_merge);
    edge_neighbors.erase(id_del);
}
bool check(const int id, EdgeNeighbors &edge_neighbors, vector<double> &diameters, vector<int> &deleted_edges,
           vector<int> &merged_edges)
{
    bool out = false;
    map<int, double> N = edge_neighbors[id];
    double d = diameters[id];
    for (auto &e1 : N)
    {
        if (find(deleted_edges.begin(), deleted_edges.end(), e1.first) != deleted_edges.end())
            break;
        int id1 = e1.first;
        double l1 = e1.second;
        double d1 = diameters[id1];
        if (d1 + d > l1)
        {
            out = true;
            diameters[id] += d1;
            merged_edges.push_back(id1);
            int id2 = common_neighbor(N, edge_neighbors[id1]);
            if (id2 != -1 && !(find(merged_edges.begin(), merged_edges.end(), id2) != merged_edges.end()))
                deleted_edges.push_back(id2);
            else
                return false;
            edge_neighbors_update(edge_neighbors, id, id1, id2);
        }
    }
    return out;
}
Graph generate_graph(const vector<double> &diameters, const vector<pair<int, int>> &edges)
{
    Graph graph;
    int it = 0;
    for (auto e : edges)
    {
        graph[e] = diameters[it];
        ++it;
    }
    return graph;
}
void update_graph(Graph &graph, vector<int> &deleted_edges,
                  vector<int> &merged_edges, const vector<pair<int, int>> &edges, const int &id, const vector<double> &diameters)
{
    set<int> excluded(deleted_edges.begin(), deleted_edges.end());
    excluded.insert(merged_edges.begin(), merged_edges.end());
    for (auto ex : excluded)
    {
        graph.erase(edges[ex]);
    }
    graph[edges[id]] = diameters[id];
}
void loop(EdgeNeighbors &edge_neighbors, vector<double> &diameters,
          int n, double constant, vector<Graph> &timeline, vector<int> &deleted_edges,
          vector<int> &merged_edges, const vector<pair<int, int>> &edges)
{
    Graph graph = generate_graph(diameters, edges);
    timeline.push_back(graph);
    int id;
    for (int it = 0; it < n; it++)
    {
        id = addRandom(diameters, deleted_edges, merged_edges, constant);
        if (id != -1)
        {
            if (check(id, edge_neighbors, diameters, deleted_edges,
                      merged_edges))
            {
                
                update_graph(graph, deleted_edges, merged_edges, edges, id, diameters);
                timeline.push_back(graph);
            }
            
        }
    }
}
void export_file(const vector<Graph> &timeline,const vector<int> &deleted_edges,const vector<int> &merged_edges,const vector<pair<int, int>> &edges)
{
    ofstream Output("timeline.txt");
    for (const auto &graph : timeline)
    {
        for (const auto &[key, diameter] : graph)
        {
            Output << "(" << key.first << ", " << key.second << "), ";
            Output << diameter << "\n";
        }
        Output<<"------------------------------------------\n";
    }
    Output.close();
    ofstream Inf("infinite_edges.txt");
    for (const auto &edge : deleted_edges)
    {
        Inf <<  edges[edge].first << ", " << edges[edge].second<<"\n";
    }
    ofstream Mer("merged_edges.txt");
    for (const auto &edge : merged_edges)
    {
        Mer <<  edges[edge].first << ", " << edges[edge].second<<"\n";
    }
}
int main(int argc, char *argv[])
{

    srand(time(nullptr));
    EdgeNeighbors edge_neighbors;
    vector<Graph> timeline;
    int num = stoi(argv[1]);
    double constant = stod(argv[2]);
    vector<double> diameters;
    vector<pair<int, int>> edges;
    vector<int> deleted_edges;
    vector<int> merged_edges;
    import(edge_neighbors, diameters, edges);
    loop(edge_neighbors,diameters,num,constant,timeline,deleted_edges,merged_edges,edges);
    export_file(timeline,deleted_edges,merged_edges,edges);

    return (0);
}