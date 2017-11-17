#include <iostream>
#include <algorithm>
#include <cmath>
#include <valarray>

using namespace std;

double f(valarray<double>& x){
    // Funcao objetivo
    return x[0]**3 + x[1]**3 + x[2]**3;
}

valarray<double> grad(valarray<double>& x){
    vector<double> resultado;
    resultado.push_back(3*(x[0]**2));
    resultado.push_back(3*(x[1]**2));
    resultado.push_back(3*(x[2]**2));
    return resultado;
}

double armijo(vector<double>& x, vector<double>& d, double g, double n ){
    double t = 1;
    while (f(soma(x, multiplica_escalar(t, d))) > f(x) + n*t*(multiplica_vetor(grad(x), d){
        t = g*t;
    }
    return t;
}

int main(){

cout << "Hello World" << endl;

return 0;
}
