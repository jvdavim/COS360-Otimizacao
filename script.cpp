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
    std::cout << "Chamei o armijo" << std::endl;
    while (f(arg) > f(x) + n*t*((grad(x)*d).sum())){
        std::cout << "Armijo t = " << t << std::endl;
        t = g*t;
        arg = x + t*d;
    }
    return t;
}

std::valarray<double> gradiente(std::valarray<double>& x){
    int k = 0;
    double epslon = 0.00927;
    std::valarray<double> px = x;
    while (grad(px).sum() > epslon){
        std::cout << "x = ( " << grad(x)[0] << ", " << grad(x)[1] << ", " << grad(x)[2] << " )" << std::endl;
        std::cout << "grad(px).sum() = " << grad(px).sum() << std::endl;
        std::valarray<double> d = grad(px)*(-1.0);
        std::cout << "-grad(px).sum() = " << d.sum() << std::endl;
        double t = armijo(px, d, 0.8, 0.25);
        std::cout << "t = " << t << std::endl;
        x = px + t*d;
        k++;
        px = x;
        getchar();
    }
    return x;
}

int main(){
double initx[] = {0, -1, 1};
// double initd[] = {1, 0, 0};
std::valarray<double> x (initx, 3);
// std::valarray<double> d (initd, 3);
// std::cout << armijo(x, d, 0.8, 0.25) << std::endl;
std::cout << "( " << gradiente(x)[0] << ", " << gradiente(x)[1] << ", " << gradiente(x)[2] << " )" << std::endl;
return 0;
}
