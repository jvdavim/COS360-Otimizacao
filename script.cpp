#include <iostream>
#include <valarray>


double norma(std::valarray<double> x){
    // Funcao norma. Calcula a norma do vetor x
    return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

double f(std::valarray<double> x){
    // Funcao objetivo
    return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
}

std::valarray<double> grad(std::valarray<double> x){
    // Funcao gradiente de f() no ponto x
    std::valarray<double> resultado (3);
    resultado[0] = 3*pow(x[0], 2);
    resultado[1] = 3*pow(x[1], 2);
    resultado[2] = 3*pow(x[2], 2);
    return resultado;
}

double armijo(std::valarray<double> x, std::valarray<double> d, double g, double n ){
    // Funcao de armijo. Retorna um tamanho de passo
    double t = 1;
    while (f(x + t*d) > f(x) + n*t*((grad(x)*d).sum())){
        t = g*t;
    }
    return t;
}

std::valarray<double> gradiente(std::valarray<double> x){
    // Metodo do gradiente. Retorna um ponto estacionario usando o passo calculado por armijo
    double epslon = 0.000000927;
    double t = 1;
    std::valarray<double> d;
    while (t > epslon){
        d = grad(x)*(-1.0);
        t = armijo(x, d, 0.8, 0.25);
        x = x + t*d;
        if (grad(x).sum() < 1.8*pow(10,307)){
            std::cout << "ERRO: Funcao nao converge para ponto estacionario" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    return x;
}

int main(){
    // Testes
    double initx[] = {-2, 1, 1};
    std::valarray<double> x (initx, 3);
    std::cout << "( " << gradiente(x)[0] << ", " << gradiente(x)[1] << ", " << gradiente(x)[2] << " )" << std::endl;
    return 0;
}
