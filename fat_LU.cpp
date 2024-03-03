// LU Factorization using numeric methods

#include <iostream>
#include <vector> 
#include <iomanip>
#include <cmath>
using namespace std;

pair<vector<vector<double>>, vector<vector<double>>> LU(vector<vector<double>> matrix) {
    int n = matrix.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U = matrix;

    for (int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(U[i][i]) < 1e-8) {
                std::cout << "Error: Dvision by zero found.\n";
                return {};
            }
            double factor = U[k][i] / U[i][i];
            L[k][i] = factor;
            for (int j = i; j < n; ++j) {
                U[k][j] -= factor * U[i][j];
            }
        }
    }

    return make_pair(L, U);
}

vector<double> solveLyb(const vector<vector<double>>& L, const vector<double>& b) {
    int n = L.size();
    vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    return y;
}

vector<double> solveUxy(const vector<vector<double>>& U, const vector<double>& y) {
    int n = U.size();
    vector<double> x(n);

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

int main() {
    int n_var, n_eq;
    cout << "Number of variables: \n";
    cin >> n_var;

    cout << "Number of equations: \n";
    cin >> n_eq;

    vector<vector<double>> matrix(n_eq, vector<double>(n_var + 1));
    for (int i = 0; i < n_eq; ++i) {
        cout << "Equation " << i + 1 << "\n";
        for (int j = 0; j < (n_var + 1); ++j) {
            if (j == n_var) {
                cout << "Independent term: \n";
            } else {
                cout << "Variable " << j + 1 << "\n";
            }
            cin >> matrix[i][j];
        }
    }

    auto result = LU(matrix);
    vector<vector<double>> L = result.first;
    vector<vector<double>> U = result.second;

    std::cout << "\nMatrix L:\n";
    for (const auto& row : L) {
        for (double value : row) {
            std::cout << std::setw(10) << value << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nMatrix U:\n";
    for (const auto& row : U) {
        for (double value : row) {
            std::cout << std::setw(10) << value << " ";
        }
        std::cout << "\n";
    }

    // Resolving Ly = b
    vector<double> b(n_eq);
    for (int i = 0; i < n_eq; ++i) {
        b[i] = matrix[i][n_var];
    }

    vector<double> y = solveLyb(L, b);

    // Resolving Ux = y
    vector<double> x = solveUxy(U, y);

    // vetor resposta
    std::cout << "\nRoot Vector (x, y, z...):\n";
    for (double value : x) {
        std::cout << std::setw(10) << value << " ";
    }
    std::cout << "\n";

    return 0;
}
