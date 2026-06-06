#ifndef TSP_H
#define TSP_H
#include <vector>

using namespace std;

struct Solution {
    vector<int> tour;
    double cost;
};

struct Subsequence {
    double T, C;
    int W;
    int first, last;
       
    static Subsequence Concatenate (Subsequence &sigma_1, Subsequence &sigma_2);
};

extern int n;
extern vector<vector<double>> dist;

double calculateCost(const vector<int>& tour);
Solution Construcao();
void UpdateAllSubseq(Solution *s, vector<vector<Subsequence>> &subseq_matrix);
void BuscaLocal(Solution *s, vector <vector<Subsequence>> &subseq_matrix);
Solution ILS(int maxIter, int maxIterIls);

//att?
#endif
