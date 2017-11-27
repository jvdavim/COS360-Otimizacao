#include <iostream>
#include <valarray>


double f(std::valarray<double> x){
    // Funcao objetivo
    return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
}

double P(std::valarray<double> x){
	double fator = pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2);
	return fator-2*sqrt(fator) + 1;
}

double fi(std::valarray<double> x, double rho=1.0){
	// Funcao com penalidade exterior
	return f(x) + rho*P(x);
}

double norma(std::valarray<double> x){
    return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

std::valarray<double> grad(std::valarray<double> x, double rho=1.0, bool penExt=false){
    // Funcao gradiente de f() no ponto x
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

double armijo(std::valarray<double>& x, std::valarray<double>& d, double g, double n ){
    // Funcao de armijo. Retorna um tamanho de passo
    double t = 1.0;
    // std::valarray<double> arg = x + t*d;
    // std::cout << "Chamei o armijo" << std::endl;
    while (f(x + t*d) > f(x) + n*t*((grad(x)*d).sum())){
        // std::cout << "Armijo t = " << t << std::endl;
        t = g*t;
    }
    return t;
}

// std::valarray< std::valarray<double> > hessiana(std::valarray<double> x){
//     // Calcula a Hessiana
//     double l1[3] = {6*x[0], 0, 0}; //linha 1 da matriz
//     double l2[3] = {0, 6*x[1], 0}; //linha 2 da matriz
//     double l3[3] = {0, 0, 6*x[2]}; //linha 3 da matriz
//     std::valarray< std::valarray<double> > M (3);
//     std::valarray<double> a (l1, 3);
//     std::valarray<double> b (l2, 3);
//     std::valarray<double> c (l3, 3);
//     //===
//     M[0] = a;
//     M[1] = b;
//     M[2] = c;

//     return M;
// }

// std::valarray< std::valarray<double> > getHessiana(std::valarray<double> x, double rho=1.0){
//     std::valarray< std::valarray<double> > M (3);
//     double l1[3];
//     double l2[3];
//     double l3[3];
//     double div = pow((pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)), 1.5);
//     l1[0] = 6*x[0] + 2*rho*(1 - (pow(x[1],2) + pow(x[2],2))/div);
//     l1[1] = 2*rho*x[0]*x[1]/div;
//     l1[2] = 2*rho*x[0]*x[2]/div;
//     l2[0] = 2*rho*x[1]*x[0]/div;
//     l2[1] = 6*x[1] + 2*rho*(1 - (pow(x[0],2) + pow(x[2],2))/div);
//     l2[2] = 2*rho*x[1]*x[2]/div;
//     l3[0] = 2*rho*x[2]*x[0]/div;
//     l3[1] = 2*rho*x[2]*x[1]/div;
//     l3[2] = 6*x[2] + 2*rho*(1 - (pow(x[0],2) + pow(x[1],2))/div);

//     std::valarray<double> a (l1, 3);
//     std::valarray<double> b (l2, 3);
//     std::valarray<double> c (l3, 3);
//     M[0] = a;
//     M[1] = b;
//     M[2] = c;

//     return M;
// }

std::valarray< std::valarray<double> > multVetor(
    std::valarray<double> V1, std::valarray<double> V2){
    //Multipica V1 por V2 transposto, retorna uma matriz
    double l1[3] = {V1[0]*V2[0], V1[0]*V2[1], V1[0]*V2[2]}; //linha 1 da matriz
    double l2[3] = {V1[1]*V2[0], V1[1]*V2[1], V1[1]*V2[2]}; //linha 2 da matriz
    double l3[3] = {V1[2]*V2[0], V1[2]*V2[1], V1[2]*V2[2]}; //linha 3 da matriz
    std::valarray< std::valarray<double> > M (3);
    std::valarray<double> a (l1, 3);
    std::valarray<double> b (l2, 3);
    std::valarray<double> c (l3, 3);
    //===
    M[0] = a;
    M[1] = b;
    M[2] = c;

    return M;
}

std::valarray<double> multVetMatriz(
     std::valarray<double> V1,
     std::valarray< std::valarray<double> > M1,
     int n =3){
     //Multiplica um vetor linha por uma matriz
     //Resulta em um vetor linha
     std::valarray<double> V (n);
     for (int i=0; i<n; i++){
        V[i] = 0;
        for (int j=0; j<n; j++){
            V[i] += V1[j]*M1[j][i];
        }
     }
     return V;
}


std::valarray<double> multMatVet(
	std::valarray< std::valarray<double> > M1,
    std::valarray<double> V1,
    int n =3){
	std::valarray<double> V (n);
	for (int i=0; i<n; i++){
		V[i] = (M1[i]*V1).sum();
	}
	return V;
}


std::valarray< std::valarray<double> > multMatriz(
    std::valarray< std::valarray<double> > M1,
    std::valarray< std::valarray<double> > M2,
    int n =3){
    //Multipica duas matrizes M1 e M2
    std::valarray< std::valarray<double> > M (n);
    std::valarray<double> temp (n);
    double aux[n];
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            aux[j] = 0;
            for (int k=0; k<n; k++){
                aux[j] += M1[i][k]*M2[k][j];
            }
        temp = {aux[0], aux[1], aux[2]};
        M[j] = temp;
        }
    }
    return M;
}

std::valarray< std::valarray<double> > divMatriz(
    std::valarray< std::valarray<double> > M1,
    double divisor=1, int n=3){
    for (int i=0; i<n; i++){
        M1[i] = M1[i]/divisor;
    }
    return M1;
}

std::valarray< std::valarray<double> > somaMatriz(
    std::valarray< std::valarray<double> > M1,
    std::valarray< std::valarray<double> > M2,
    int n=3){
    std::valarray< std::valarray<double> > M (3);
    for (int i=0; i<n; i++){
        M[i] = M1[i]+M2[i];
    }
    return M;
}

std::valarray< std::valarray<double> > subMatriz(
    std::valarray< std::valarray<double> > M1,
    std::valarray< std::valarray<double> > M2,
    int n=3){
    std::valarray< std::valarray<double> > M (3);
    for (int i=0; i<n; i++){
        M[i] = M1[i]-M2[i];
    }
    return M;
}

std::valarray< std::valarray<double> > multEscalar(
    std::valarray< std::valarray<double> > M1,
    double a,
    int n=3){
    for (int i=0; i<n; i++){
        M1[i] = M1[i]*a;
    }
    return M1;
}

std::valarray< std::valarray<double> > getIdentidade(int n = 3){
	std::valarray< std::valarray<double> > M (n);
	std::valarray<double> l (0.0,n);
	for (int i=0; i<n; i++){
		M[i] = l;
		M[i][i] = 1.0;
	}
	return M;
}

void printMatriz(std::valarray< std::valarray<double> > M, int n=3){
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){ 
        std::cout << M[i][j] <<" ";
        }
        std::cout <<""<<std::endl;
    }
}

void printVetor(std::valarray <double> V, int n=3){
    std::cout<<"[";
    for (int i=0; i<n; i++){
        std::cout << V[i] << " ";
    }
    std::cout << "]"<<std::endl;
}


std::valarray< std::valarray<double> > bfgs(
	std::valarray< std::valarray<double> > h,
	std::valarray<double> pk,
	std::valarray<double> qk){
	//Utiliza o metodo BFGS para calcular a aproximacao da inversa da hessiana
	//multVetor(pk,qk) retorna uma matriz que é multiplicada por outra matriz
	//Depois efetua-se uma soma de matrizes e por fim, dividi-se a matriz por um escalar 
	double denominador = (pk*qk).sum();
	double numerador1 = (multVetMatriz(qk,h)*qk).sum();
	std::valarray< std::valarray<double> > numerador2 = multEscalar(multVetor(pk,pk), (1+numerador1/denominador));
	std::valarray< std::valarray<double> > numerador3 = somaMatriz(multMatriz(multVetor(pk,qk),h), multMatriz(h, multVetor(qk,pk)));
	return subMatriz(somaMatriz(h, divMatriz(numerador2, denominador)), divMatriz(numerador3, denominador));

}

std::valarray<double> qNewton(std::valarray<double> x, double rho = 1.0, double epslon = 10e-6, int maxIter = 10000){
    std::valarray< std::valarray<double> > hess (3);

    // hess = hessiana(x);
    // hess = getHessiana(x);
    hess = getIdentidade(3); //hessiana inicial assume valor matriz identidade
    double t;
    int nIter = 0;
    std::valarray<double> gf (3); //vetor gradiente
    std::valarray<double> d (1.0, 3); //direção de descida
    std::valarray<double> x0 (3); //ponto x anterior
    std::valarray<double> p (3);
    std::valarray<double> q (3);
    while (norma(grad(x)) > epslon || nIter < maxIter){
        gf = grad(x);
        
		d = (-1.0)*multMatVet(hess, gf);
        x0 = x; //guarda o x anterior antes de atualizar
        x+= armijo(x, d, 0.8, 0.25)*d;
        // printMatriz(hess);
        // printVetor(d);
        // printVetor(x);
        // getchar();

        // //Calcular 'p' e 'q'
        p = x-x0; //diferença entre x atual e x anterior (x0)
        q = grad(x)-grad(x0); //diferença entre gradiente calculado no ponto x atual e no ponto x anterior
        hess = bfgs(hess, p, q);
        nIter++;
    }
    printVetor(x);
    return x;

}


int main(){
    double initx[] = {0, -1, 1};
    std::valarray<double> x (initx, 3);

    qNewton(x);

return 0;
}
