#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;


vector<double> soma(vector<double>& a, vector<double>& b){
    if (a.size() != b.size()) return -1;
    vector<double> resultado;
    for (int i = 0; i < a.size(); i++){
        resultado.push_back(a[i]+b[i]);
    }
    return resultado;
}

vector<double> multiplica_escalar(int t, vector<double>& x){
    vector<double> resultado;
    for (int i = 0; i < x.size(); i++){
        resultado.push_back(t*x[i]);
    }
    return resultado;
}

double multiplica_vetor(vector<double>& a, vector<double>& b){
    if (a.size() != b.size()) return -1;
    double resultado;
    for (int i = 0; i < a.size(); i++){
        resultado += a[i]*b[i];
    }
    return resultado;
}

double f(vector<double>& x){
    // Funcao objetivo
    return x[0]**3 + x[1]**3 + x[2]**3;
}

vector<double> grad(vector<double>& x){
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
