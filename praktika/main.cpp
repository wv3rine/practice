#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <queue>

class Graph {
public:
    int n, k, s, t;
    // Здесь возможно третий int сделать double, потому что веса double
    std::unordered_map<int, std::unordered_map<int, int>> g;
    // ДОЛЖЕН БЫТЬ private И const
    std::vector<std::pair<int, int>> blocked_edges;

    bool isBlocked(int u, int v) const {
        for (const auto& [first, second] : blocked_edges) {
            if ((first == u && second == v) || (first == v && second == u)) {
                return true;
            }
        }
        return false;
    }
};

struct Path {
    int64_t c;
    std::vector<int> path;

    Path(const std::vector<int>& path, int c) {
        this->path = path;
        this->c = c;
    }

    bool operator<(const Path& other) const {
        return c < other.c;
    }
};

// Готовый алгоритм Дейкстры, которые я сам сделал и раньше использовал для учебы
static int64_t dijkstra(Graph& graph) {
    int64_t kInf = 1000000000000, n = graph.n, s = graph.s, t = graph.t;
    std::vector<int64_t> dist = std::vector<int64_t>(n, kInf);
    dist[s] = 0;
    std::priority_queue<std::pair<int64_t, int>> q;
    q.emplace(0, s);
    while (!q.empty()) {
        int v = q.top().second;
        int64_t cur_d = -q.top().first;
        q.pop();
        if (cur_d > dist[v]) {
            continue;
        }

        for (const auto& [to, w] : graph.g[v]) {
            if (!graph.isBlocked(v, to)) {
                if (dist[v] + w < dist[to]) {
                    dist[to] = dist[v] + w;
                    q.emplace(-dist[to], to);
                }
            }
        }
    }
    return dist[t];
}

// Костыль - тут std::vector<PATH>
bool isSimilarCostsPropertyFulfilled(const std::vector<Path>& vec, int64_t x) {
    // Скопировано для удобства
    std::vector<int64_t> c(vec.size() + 1);
    for (int i = 0; i < vec.size(); ++i) {
        c[i] = vec[i].c;
    }
    c[vec.size()] = x;
    int n = vec.size() + 1;

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += 2/static_cast<double>(n) * c[j];
        }
        if (static_cast<double>(c[i]) > sum) {
            return false;
        }
    }
    return true;
}

// bool - дошел ли, int - длина пути
// Придется копировать, чтобы удалять
std::pair<bool, int64_t> mainAlgorithm(std::vector<Path> path, int k, Graph& graph) {
    std::random_device rd;
    std::mt19937 gen(rd());

    int64_t d = 0;
    while (!path.empty()) {
        // Рассчитываем вероятности (НЕНОРМИРОВАННЫЕ)
        std::vector<double> p(path.size(), 0);
        for (int i = 0; i < path.size(); ++i) {
            for (int j = 0; j < path.size(); ++j) {
                if (j == i) {
                    p[i] += (1 - k) * path[i].c;
                } else {
                    p[i] += 2 * path[j].c;
                }
            }
            p[i] /= static_cast<double>(k + 1) * (k + 1) * path[i].c;
        }

        // Распределение
        std::discrete_distribution<> distribution(p.begin(), p.end());

        // Просто идем по пути, и когда видим заблокированное ребро,
        // добавляем пройденное расстояние (то есть возвращаемся в s)
        int path_id = distribution(gen), current_d = 0;
        size_t i, path_length = path[path_id].path.size();
        std::cout << path[path_id].path[0] << " ";
        for (i = 1; i < path_length; ++i) {
            int u = path[path_id].path[i - 1];
            int v = path[path_id].path[i];
            if (!graph.isBlocked(u, v)) {
                d += graph.g[u][v];
                current_d += graph.g[u][v];
                std::cout << v << " ";
            } else {
                d += current_d;
                std::cout << " ...returning to " << graph.s << "\n";
                break;
            }
        }
        // дошли до t
        if (i == path[path_id].path.size()) {
            std::cout << "\n";
            break;
        }
        path.erase(path.begin() + path_id);
    }
    // Если path.empty(), то путешественник не дошел до t
    return {path.empty(), d};
}

void printGroups(const std::vector<std::vector<Path>>& group) {
    std::cout << "--------GROUPS--------\n";
    for (auto i = 0; i < group.size(); ++i) {
        std::cout << "Group " << i << "\n";
        for (auto x : group[i]) {
            for (auto p : x.path) {
                std::cout << p << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "----------------------\n";
}

int64_t algo(Graph& graph) {
    int s = graph.s, t = graph.t;
    std::vector<Path> path;

    // Заполнение path
    for (const auto& [start, weight] : graph.g[s]) {
        path.emplace_back(std::vector<int>({s, start}), weight);
        int u = start;
        int prev = s;
        while (u != t) {
            // Берем соседа u, который не вершина, из которой мы пришли
            for (const auto& [v, w] : graph.g[u]) {
                if (v != prev) {
                    prev = u;
                    u = v;
                    path[path.size() - 1].path.push_back(u);
                    path[path.size() - 1].c += w;
                    break;
                }
            }
            // здесь может быть какая-то отработка вершины u
        }
    }

    // Сортировка по возрастанию весов
    std::sort(path.begin(), path.end());

    // Группы
    std::vector<std::vector<Path>> group;
    // Номер текущей группы
    int j = 0, i = 0;
    while (i < path.size()) {
        group.emplace_back();
        group[j].push_back(path[i]);
        ++i;
        while (i < path.size() && isSimilarCostsPropertyFulfilled(group[j], path[i].c)) {
            group[j].push_back(path[i]);
            ++i;
        }
        ++j;
    }

    // основной алгоритм
    int l = 0;
    int64_t answer = 0;
    auto pair = mainAlgorithm(group[l], graph.k, graph);
    // не будет выхода за границы у l, потому что гарантированно есть k + 1 путь
    while (pair.first) {
        answer += pair.second;
        ++l;
        pair = mainAlgorithm(group[l], graph.k, graph);
    }
    answer += pair.second;
    return answer;
}

// из 0 в 1 вершину
Graph generateSpecialGraph(int min_c, int max_c, int k, int path_count) {
    // должны быть разные проверки

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> gen_c(min_c, max_c);
    std::uniform_int_distribution<int> gen_w(1, 1000000);

    Graph graph = Graph();
    graph.k = k;
    graph.s = 0;
    graph.t = 1;
    int cur = 2;

    for (auto i = 0; i < path_count; ++i) {
        // заполнение i-го путя
        int prev = graph.s;
        for (auto j = 0; j < gen_c(mt); ++j) {
            int w = gen_w(mt);
            graph.g[prev][cur] = w;
            graph.g[cur][prev] = w;
            prev = cur;
            ++cur;
        }
        int w = gen_w(mt);
        graph.g[cur - 1][graph.t] = w;
        graph.g[graph.t][cur - 1] = w;
    }
    graph.n = cur;

    // Распределение blocked ребер
    // Могут быть совпадения, но blocked ребер все равно < k
    for (auto i = 0; i < k; ++i) {
        std::uniform_int_distribution<int> gen(0, cur - 1);
        int u = gen(mt);
        for (const auto& [v, w] : graph.g[u]) {
            if (!graph.isBlocked(u, v)) {
                graph.blocked_edges.emplace_back(u, v);
                graph.blocked_edges.emplace_back(v, u);
                break;
            }
        }
    }
    return graph;
}

void printGraph(Graph& g) {
    std::cout << "--------GRAPH-------\n";
    for (const auto& [u, other] : g.g) {
        for (const auto& [v, w] : other) {
            std::cout << u << " " << v << " - " << w << "\n";
        }
    }
    std::cout << "\n\n";
    for (const auto& [u, v] : g.blocked_edges) {
        std::cout << u << " " << v << "\n";
    }
    std::cout << "\n\n";
    std::cout << "--------------";
}

int main() {
    // mode == 0 - ввод с клавиатуры, mode == 1 - случайная генерция
    int mode = 0;


    Graph g;
    int64_t alg_sum = 0, opt_sum;
    if (mode == 0) {
        g = Graph();
        int m;
        std::cout << "Enter n, m, k, s, t\n";
        std::cin >> g.n >> m >> g.k >> g.s >> g.t;
        std::cout << "\nEnter <u v w> m times\n";
        for (auto i = 0; i < m; ++i) {
            int u, v, w;
            std::cin >> u >> v >> w;
            g.g[u][v] = w;
            g.g[v][u] = w;
        }
        std::cout << "\nEnter <u v> k times\n";
        for (auto i = 0; i < g.k; ++i) {
            int u, v;
            std::cin >> u >> v;
            g.blocked_edges.emplace_back(u, v);
            g.blocked_edges.emplace_back(v, u);
        }
        std::cout << "\n";
    } else if (mode == 1) {
        g = generateSpecialGraph(100, 200, 500, 501);
    }
    opt_sum = dijkstra(g);
    std::cout << opt_sum << "\n";
    int N = 1;
    for (auto i = 0; i < N; ++i) {
        alg_sum += algo(g);
        //std::cout << algo(g) << "\n";
        //std::cout << alg_sum << " " << opt_sum << "\n";
    }
    std::cout << "OPT - " << opt_sum << "\nALG - " << alg_sum << "\n";
    std::cout << "competitive-ratio - (k + 1) - " << static_cast<double>(alg_sum) / N / opt_sum - g.k - 1;
    return 0;
}
