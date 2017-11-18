#include <iostream>
#include <valarray>


double f(std::valarray<double> x){
    // Funcao objetivo
    return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
}

double norma(std::valarray<double> x){
    return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

std::valarray<double> grad(std::valarray<double> x){
    // Funcao gradiente de f() no ponto x
    std::valarray<double> resultado (3);
    resultado[0] = 3*pow(x[0], 2);
    resultado[1] = 3*pow(x[1], 2);
    resultado[2] = 3*pow(x[2], 2);
    return resultado;
}

double armijo(std::valarray<double>& x, std::valarray<double>& d, double g, double n ){
    // Funcao de armijo. Retorna um tamanho de passo
    double t = 1.0;
    // std::valarray<double> arg = x + t*d;
    std::cout << "Chamei o armijo" << std::endl;
    while (f(x + t*d) > f(x) + n*t*((grad(x)*d).sum())){
        std::cout << "Armijo t = " << t << std::endl;
        t = g*t;
    }
    return t;
}

std::valarray< std::valarray<double> > hessiana(std::valarray<double> x){
    // Calcula a Hessiana
    double l1[3] = {6*x[0], 0, 0}; //linha 1 da matriz
    double l2[3] = {0, 6*x[1], 0}; //linha 2 da matriz
    double l3[3] = {0, 0, 6*x[2]}; //linha 3 da matriz
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

std::valarray< std::valarray<double> > multVetor(std::valarray<double> V){
    //Multipica V por V transposto, retorna uma matriz
    double l1[3] = {V[0]*V[0], V[0]*V[1], V[0]*V[2]}; //linha 1 da matriz
    double l2[3] = {V[1]*V[0], V[1]*V[1], V[1]*V[2]}; //linha 2 da matriz
    double l3[3] = {V[2]*V[0], V[2]*V[1], V[2]*V[2]}; //linha 3 da matriz
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

// std::valarray< std::valarray<double> > multMatriz(
//     std::valarray< std::valarray<double> > M1,
//     std::valarray< std::valarray<double> > M2){
//     //Multipica duas matrizes M1 e M2
//     std::valarray< std::valarray<double> > M (3);
//     std::valarray<double> temp (3);
//     double aux[3];
//     for (int i=0; i<3; i++){
//         for (int j=0; j<3; j++){
//             aux[j] = M1[i]*M2[j];
//         }

//         temp = {aux[0], aux[1], aux[2]};
//         M[0] = temp;
//     }
//     return M;
// }

void printMatrix3x3(std::valarray< std::valarray<double> > M){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){ 
        std::cout << M[i][j] <<" ";
        }
        std::cout <<""<<std::endl;
    }
}

// std::valarray<double> hessiana(valarray<double> x){
//     // Atualiza a Hessiana
//     valarray<double> p = 
// }

std::valarray<double> qNewton(std::valarray<double> x){
    std::valarray< std::valarray<double> > hess (3);
    std::valarray< std::valarray<double> > m1 (3);
    hess = hessiana(x);

    std::valarray<double> gf (3);
    std::valarray<double> d (3);
    std::valarray<double> x0 (3);
    double t;
    double epslon = 10e-6;
    std::valarray<double> p (3);
    std::valarray<double> q (3);
    while (norma(grad(x)) > epslon){
        gf = grad(x);
         d = {hess[0][0]*gf[0]+hess[0][1]*gf[1]+hess[0][2]*gf[2],
              hess[1][0]*gf[0]+hess[1][1]*gf[1]+hess[1][2]*gf[2],
              hess[2][0]*gf[0]+hess[2][1]*gf[1]+hess[2][2]*gf[2]};
        t = armijo(x, d, 0.8, 0.25);
        x0 = x;
        x = x + t*d;
        //Calcular 'p' e 'q'
        p = x-x0;
        q = grad(x)-grad(x0);

        //Atualizar a Hessiana

        m1 = multVetor(p); 





    }
}


int main(){
double initx[] = {0, -1, 1};
    std::valarray<double> x (initx, 3);
    // std::valarray< std::valarray<double> > hess (3);
    // hess = hessiana(x);
    // printMatrix3x3(hess);
    // std::valarray<double> g (3);
    // g = grad(x);
    // std::cout<< "HELLO"<<std::endl;
    // std::cout<< hess[1][1]*g[1]<<std::endl;

return 0;
}
