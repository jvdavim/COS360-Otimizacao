#include <iostream>
#include <valarray>


double f(std::valarray<double>& x){
    // Funcao objetivo
    return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
}

std::valarray<double> grad(std::valarray<double>& x){
    // Funcao gradiente de f() no ponto x
    std::valarray<double> resultado (3);
    resultado[0] = 3*pow(x[0], 2);
    resultado[1] = 3*pow(x[1], 2);
    resultado[2] = 3*pow(x[2], 2);
    return resultado;
}

double armijo(std::valarray<double>& x, std::valarray<double>& d, double g, double n ){
    // Funcao de armijo. Retorna um tamanho de passo
    double t = 1;
    std::valarray<double> arg = x + t*d;
    while (f(arg) > f(x) + n*t*((grad(x)*d).sum())){
        t = g*t;
        arg = x + t*d;
    }
    return t;
}

int main(){
double initx[] = {0, -1, 1};
double initd[] = {1, 0, 0};
std::valarray<double> x (initx, 3);
std::valarray<double> d (initd, 3);
std::cout << armijo(x, d, 0.8, 0.25) << std::endl;
return 0;
}
