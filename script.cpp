#include <iostream>
#include "op.hpp"

/* Função objetivo */
double f(std::valarray<double> x)
{
    
    return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
}


/* Função de penalidade */
double P(std::valarray<double> x)
{
	double fator = pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2);

	return fator-2*sqrt(fator) + 1;
}


/* Função objetivo com penalidade exterior */
double fi(std::valarray<double> x, double rho=1.0)
{
	return f(x) + rho*P(x);
}


/* Operador norma */
double norma(std::valarray<double> x)
{
    return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}


/* Função gradiente de f() no ponto x */
std::valarray<double> grad(std::valarray<double> x, double rho=1.0, bool penExt=false)
{
    std::valarray<double> resultado (3);

    if (penExt==true){
    	double fator = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
    	resultado[0] = 3*pow(x[0],2) + 2*rho*x[0]*(1 - 1/fator);
    	resultado[1] = 3*pow(x[1],2) + 2*rho*x[1]*(1 - 1/fator);
    	resultado[2] = 3*pow(x[2],2) + 2*rho*x[2]*(1 - 1/fator);
    }

    else {
		resultado[0] = 3*pow(x[0], 2);
	    resultado[1] = 3*pow(x[1], 2);
	    resultado[2] = 3*pow(x[2], 2);
    }

    return resultado;
}


/* Método de armijo */
double armijo(std::valarray<double>& x, std::valarray<double>& d, double g, double n, double rho )
{
    double t = 1.0;

    while (fi(x + t*d, rho) > fi(x, rho) + n*t*((grad(x, rho, true)*d).sum())){
        t = g*t;
    }

    return t;
}


/* Método de BFGS -- Atualiza a matriz hessiana no método de Quase-Newton */
std::valarray< std::valarray<double> > bfgs(std::valarray< std::valarray<double> > h, std::valarray<double> p, std::valarray<double> q)
{
	//Utiliza o metodo BFGS para calcular a aproximacao da inversa da hessiana
	//multVetor(pk,qk) retorna uma matriz que é multiplicada por outra matriz
	//Depois efetua-se uma soma de matrizes e por fim, dividi-se a matriz por um escalar
	double denominador = (p*q).sum();
	double numerador1 = (multVetMatriz(q,h)*q).sum();
	std::valarray< std::valarray<double> > numerador2 = multEscalar(multVetor(p,p), (1+numerador1/denominador));
	std::valarray< std::valarray<double> > m1 = multMatriz(multVetor(p,q),h);
	std::valarray< std::valarray<double> > m2 = multMatriz(h,multVetor(p,q));
	std::valarray< std::valarray<double> > numerador3 = somaMatriz(multMatriz(multVetor(p,q),h), multMatriz(h, multVetor(q,p)));

	return subMatriz(somaMatriz(h, divMatriz(numerador2, denominador)), divMatriz(numerador3, denominador));

}


/* Método de Quase-Newton */
std::valarray<double> qNewton(std::valarray<double> x, double rho = 1.0, double epslon = 10e-6, int maxIter = 10000)
{
    std::valarray< std::valarray<double> > hess (3);
    hess = getIdentidade(3); //hessiana inicial assume valor matriz identidade
    std::valarray<double> gf (3); //vetor gradiente
    std::valarray<double> d (3); //direção de descida
    std::valarray<double> x0 (3); //ponto x anterior
    std::valarray<double> p (3);
    std::valarray<double> q (3);
    int nIter = 0;

    while (norma(grad(x, rho, true)) > epslon || nIter < maxIter){
        gf = grad(x, rho, true);
		d = (-1.0)*multMatVet(hess, gf);
        x0 = x; //guarda o x anterior antes de atualizar
        x+= armijo(x, d, 0.8, 0.25, rho)*d;

        // Calcular 'p' e 'q'
        p = x-x0; //diferença entre x atual e x anterior (x0)
        q = grad(x, rho, true)-grad(x0, rho, true); //diferença entre gradiente calculado no ponto x atual e no ponto x anterior

        // Atualizar hessiana
        hess = bfgs(hess, p, q);

        // Atualizar rho
        rho= rho*3; // Beta = 3

        nIter++; //incrementa número de iterações
    }

    return x;
}


/* Método do gradiente -- Retorna um ponto estacionário usando o passo calculado por armijo */
std::valarray<double> gradiente(std::valarray<double> x, double rho = 1.0, double epslon = 10e-6, int maxIter = 10000)
{
    double t = 1;
    std::valarray<double> d;

    while (norma(grad(x)) > epslon){
        d = grad(x)*(-1.0);
        t = armijo(x, d, 0.8, 0.25, rho);
        x = x + t*d;
        
        if (norma(grad(x)) > 1.8*pow(10,307)){
            std::cout << "Função não converge para ponto estacionário" << std::endl;
            break;
        }
    }

    return x;
}


int main(){
    /* Testes Quase-Newton */
    double ninitx[] = {1, -1, 1};
    std::valarray<double> nx (ninitx, 3);

    printVetor(qNewton(nx, 1.0));


    /* Testes Gradiente */
    double ginitx[] = {1, -1, 1};
    std::valarray<double> gx (ginitx, 3);

    return 0;
}
