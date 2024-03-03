// Implementando o método da Eliminação de Gauss
// Implementing the Gaussian Elimination method

#include <iostream>
#include <vector> 
#include <iomanip> // library for output handling
#include <cmath>
using namespace std;
 
// Gauss Elimination
vector<double> gaussElimination(vector<vector<double>> matrix) {
    int row = matrix.size();
    int col = matrix[0].size() - 1; // The last element of each row is the independent term

    for (int i = 0; i < row - 1; ++i) {
        if (matrix[i][i] == 0) {
            // swap with a non-null row below, if necessary
            for (int j = i + 1; j < col; ++j) {
                if (matrix[j][i] != 0) {
                    swap(matrix[i], matrix[j]);
                    break;
                }
            }
        }

        if (matrix[i][i] != 0) {
            for (int k = i + 1; k < row; ++k) {
                double factor = matrix[k][i] / matrix[i][i];
                for (int j = i; j < col + 1; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }
    }

    // Solving the upper triangular system
    vector<double> solution(row);
    for (int i = row - 1; i >= 0; --i) {
        solution[i] = matrix[i][col] / matrix[i][i];
        for (int k = i - 1; k >= 0; --k) {
            matrix[k][row] -= matrix[k][i] * solution[i];
        }
    }

    return solution;
}

int main() {
    cout << std::setprecision(8);
    int n_var, n_eq;
    cout << "Number of variables: \n";
    cin >> n_var;

    cout << "Number of equations: \n";
    cin >> n_eq;

    vector<vector<double>> matrix(n_eq, vector<double>(n_var+1));
    for (int i = 0; i < n_eq; ++i) {
        cout << "Equation " << i+1 << "\n";
        for (int j = 0; j < (n_var+1); ++j) {
            if (j == n_var) {
                cout << "Independent term: \n";
            }
            else {
                cout << "Variable " << j+1 << "\n";
            }
            cin >> matrix[i][j];
        }
    }

    vector<double> result = gaussElimination(matrix);
    
    cout << "Results: ";
    const double almost_zero = 1.0e-4;

    for (size_t i = 0; i < result.size(); ++i) {
        if (abs(result[i]) < almost_zero) { // Checks if the value is "almost zero"  
            cout << 0 << " "; // Shows 0 if "almost zero"
        } else {
            cout << result[i] << " ";
        }
    }
    
    return 0;
}
