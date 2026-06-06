#include "tsp.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <iostream>

int n;
vector<vector<double>> dist;

inline Subsequence Subsequence::Concatenate(Subsequence &sigma_1, Subsequence &sigma_2) {
    Subsequence sigma;
    double temp = dist[sigma_1x.last][sigma_2.first];
    
    sigma.W = sigma_1.W + sigma_2.W;
    sigma.T = sigma_1.T + temp + sigma_2.T;
    sigma.C = sigma_1.C + sigma_2.W * (sigma_1.T + temp) + sigma_2.C;
    sigma.first = sigma_1.first;
    sigma.last = sigma_2.last;

    return sigma;
}

double calculateCost(const vector<int>& tour) {
    double cost = 0.0;
    double tempo_acumulado = 0.0;
    for(int i = 0; i < tour.size()-1; i++) {
        tempo_acumulado += dist[tour[i]][tour[i+1]];
        cost += tempo_acumulado;
    }
    return cost;
}

vector <int> nosRestantes(const Solution& s) {
    vector<int> CL;
    vector<bool> visitado(n, false);
    
    for (int i = 0; i < s.tour.size(); i++) { 
        visitado[s.tour[i]] = true;
    }
    
    for (int i = 0; i < n; i++) {
        if (!visitado[i]) CL.push_back(i);
    }
    return CL;
}

Solution Construcao() {
    Solution s;
    s.tour.push_back(0);
    vector <int> CL = nosRestantes(s); //linha 4 5
    double alpha = (rand() % 26) / 100.0; //linha 2
    while (!CL.empty()) { //linha 7
        int r = s.tour.back(); 
        sort(CL.begin(), CL.end(), [&](int a, int b) { //linha 8
            return dist[r][a] < dist[r][b];
        });
        int limit = max(1, ((int) ceil(alpha*CL.size()))); //linha 9
        int selecionado = rand() % limit; //linha 9
        s.tour.push_back(CL[selecionado]); //linha 10
        CL.erase(CL.begin() + selecionado); //linha 12
    }
    s.tour.push_back(0);
    s.cost = calculateCost(s.tour);
    
    return s;
}

bool bestImprovementSwap(Solution *s, vector<vector<Subsequence>> &subseq_matrix) {
    double best_cost = s->cost;
    int best_i = -1, best_j = -1;
    int sz = s->tour.size();
    bool improved = false;
    for(int i = 1; i < sz - 1; i++) {
        for(int j = i + 1; j < sz - 1; j++) {
            Subsequence sigma_total;
            if (j == i+1) {
                Subsequence s1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j]);
                Subsequence s2 = Subsequence::Concatenate(s1, subseq_matrix[i][i]);
                sigma_total = Subsequence::Concatenate(s2, subseq_matrix[j+1][sz-1]);
            }
            else {
                Subsequence s1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j]);
                Subsequence s2 = Subsequence::Concatenate(s1, subseq_matrix[i+1][j-1]);
                Subsequence s3 = Subsequence::Concatenate(s2, subseq_matrix[i][i]);
                sigma_total = Subsequence::Concatenate(s3, subseq_matrix[j+1][sz-1]);
            }
            if (sigma_total.C < best_cost) {
            best_cost = sigma_total.C;
            best_i = i;
            best_j = j;
            improved = true;
    }
    }
}  
if (improved) {
    swap(s->tour[best_i], s->tour[best_j]);
    s->cost = best_cost;
    UpdateAllSubseq(s, subseq_matrix);
}
return improved;
}

bool bestImprovement20pt(Solution *s, vector<vector<Subsequence>> &subseq_matrix) {
    double best_cost = s->cost;
    int best_i = -1, best_j = -1;
    int sz = s->tour.size();
    bool improved = false;

    for (int i = 1; i<sz-2;i++) {
        for (int j = i+1; j<sz-1;j++) {
        Subsequence s1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][i]);
        Subsequence sigma_total = Subsequence::Concatenate(s1, subseq_matrix[j+1][sz-1]);
        if (sigma_total.C < best_cost) {
        best_cost = sigma_total.C;
        best_i = i;
        best_j = j;
        improved = true;
    }
    }
}
    if (improved) {
        reverse(s->tour.begin() + best_i, s->tour.begin() + best_j + 1);
        s->cost = best_cost;
        UpdateAllSubseq(s, subseq_matrix);
    }
    return improved; // nenhuma melhoria
}

bool bestImprovement0r0pt(Solution *s, int blocksize, vector<vector<Subsequence>> &subseq_matrix) {
    double best_cost = s->cost;
    int best_i = -1, best_j = -1;
    int sz = s->tour.size();
    bool improved = false;

    for (int i = 1; i<=sz-blocksize-2;i++) {
        for (int j = 0;j<sz-1;j++) {
            if (j >= i-1 && j <= i+blocksize-1) continue;
            Subsequence sigma_total;
            if (j <i-1) {
                Subsequence s1 = Subsequence::Concatenate(subseq_matrix[0][j], subseq_matrix[i][i+blocksize-1]);
                Subsequence s2 = Subsequence::Concatenate(s1, subseq_matrix[j+1][i-1]);
                sigma_total = Subsequence::Concatenate(s2,subseq_matrix[i+blocksize][sz-1]);
            }
            else {
                Subsequence s1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+blocksize][j]);
                Subsequence s2 = Subsequence::Concatenate(s1, subseq_matrix[i][i+blocksize-1]);
                sigma_total = Subsequence::Concatenate(s2,subseq_matrix[j+1][sz-1]);
            }
           if (sigma_total.C < best_cost) {
            best_cost = sigma_total.C;
            best_i = i;
            best_j = j;
            improved = true;
        }
        }
    }
    if (improved) {
        vector<int> block(s->tour.begin() + best_i, s->tour.begin() + best_i + blocksize);

        s->tour.erase(s->tour.begin() + best_i, s->tour.begin() + best_i + blocksize);

        if (best_j > best_i) best_j -= blocksize;

        s->tour.insert(s->tour.begin() + best_j + 1,block.begin(),block.end());
        s->cost = best_cost;
        UpdateAllSubseq(s, subseq_matrix);
    }
    return improved;
}

void BuscaLocal(Solution *s, vector<vector<Subsequence>> &subseq_matrix) {
    vector <int> NL = {1, 2, 3, 4, 5};
    bool improved = false;
    while (NL.empty() == false) {
        int n = rand() % NL.size();
        switch(NL[n]) {
             case 1:
                improved = bestImprovementSwap(s, subseq_matrix);
                break;
            case 2:
                improved = bestImprovement20pt(s, subseq_matrix);
                break;
            case 3:
                improved = bestImprovement0r0pt(s, 1, subseq_matrix);
                break;
            case 4:
                improved = bestImprovement0r0pt(s, 2, subseq_matrix);
                break;
            case 5:
                improved = bestImprovement0r0pt(s, 3, subseq_matrix);
                break;
        }
        if (improved) NL = {1, 2, 3, 4, 5};
        else NL.erase(NL.begin() + n);
    }
}

Solution Perturbacao(Solution s) {
    int sz = s.tour.size() - 1;
    int maxSize = max(2, sz/10);

    int size1 = 2 + rand() % (maxSize-1);
    int size2 = 2 + rand() % (maxSize-1);

    int i, j;

    while (true) {
        i = 1 + rand() % (sz - size1);
        j = 1 + rand() % (sz - size2);
        if (i + size1 < j || j+size2<i) break;
        }
        if (j<i) {
            swap(i,j);
            swap(size1,size2);
        }
        vector <int> block1(s.tour.begin() + i, s.tour.begin() + i + size1);
        vector <int> block2(s.tour.begin() + j, s.tour.begin() + j + size2);

        s.tour.erase(s.tour.begin() + j, s.tour.begin() + j + size2);
        s.tour.erase(s.tour.begin() + i, s.tour.begin() + i + size1);

        s.tour.insert(s.tour.begin() + i, block2.begin(), block2.end());
        s.tour.insert(s.tour.begin() + j - size1 + size2, block1.begin(), block1.end());

        return s;
    }

void UpdateAllSubseq(Solution *s, vector <vector<Subsequence>> &subseq_matrix) {
    int n = s->tour.size();

    for (int i=0;i<n;i++) {
        subseq_matrix[i][i].W = (i>0);
        subseq_matrix[i][i].C = 0;
        subseq_matrix[i][i].T = 0;
        subseq_matrix[i][i].first = s->tour[i];
        subseq_matrix[i][i].last = s->tour[i];
    }
    for (int i=0;i<n;i++) 
        for (int j=i+1;j<n;j++) 
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);
    
    for (int i=n-1;i>=0;i--) 
        for (int j=i-1;j>=0;j--) 
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]); 
} 

Solution ILS(int maxIter, int maxIterIls) {
    Solution bestOfAll;
    bestOfAll.cost = INFINITY;
    vector <vector <Subsequence>> subseq_matrix(n+2, vector <Subsequence>(n+2));
    for (int i=0;i<maxIter;i++) {
        Solution s = Construcao();
        UpdateAllSubseq(&s, subseq_matrix);
        s.cost = subseq_matrix[0][n].C;
        Solution best = s;
        int iterIls = 0;
        while (iterIls <= maxIterIls) {
            BuscaLocal(&s, subseq_matrix);
            if (s.cost < best.cost) {
                best = s;
                iterIls = 0;
            }
            s = Perturbacao(best);
            UpdateAllSubseq(&s, subseq_matrix);
            s.cost = subseq_matrix[0][n].C;
            iterIls++;
        }
        if (best.cost < bestOfAll.cost) bestOfAll = best;
        cout << "iteracao " << i+1 << "melhor custo: " << bestOfAll.cost << endl;
    }

    return bestOfAll;
    }
