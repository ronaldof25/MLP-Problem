#include "data.h"
#include "tsp.h"
#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char** argv) {
    srand(time(NULL));

    // leitura da inst�ncia usando o Data do kit
    Data data(argc, argv[1]);
    data.read();

    // preencho as vari�veis globais que o tsp.cpp usa
    n = data.getDimension();
    dist.assign(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dist[i][j] = data.getDistance(i + 1, j + 1); // Data usa �ndice 1..n

    // par�metros conforme o kit orienta
    int maxIterIls = min(100, n);

    Solution s = ILS(10, maxIterIls);

    cout << "Instancia: " << data.getInstanceName() << endl;
    cout << "Custo final: " << s.cost << endl;

    return 0;
}
