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
    // std::cout << "Chamei o armijo" << std::endl;
    while (f(x + t*d) > f(x) + n*t*((grad(x)*d).sum())){
        // std::cout << "Armijo t = " << t << std::endl;
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
        std::cout << V[i] << ", ";
    }
    std::cout << "]"<<std::endl;
}

std::valarray<double> qNewton(std::valarray<double> x){
    std::valarray< std::valarray<double> > hess (3);
    std::valarray< std::valarray<double> > m1 (3);
    std::valarray< std::valarray<double> > m2 (3);
    std::valarray< std::valarray<double> > m3 (3);
    std::valarray< std::valarray<double> > m4 (3);
    std::valarray< std::valarray<double> > m5 (3);
    std::valarray<double> v1 (3);
    hess = hessiana(x);

    std::valarray<double> gf (3);
    std::valarray<double> d (3);
    std::valarray<double> x0 (3);
    double t;
    double div;
    double aux;
    double epslon = 10e-6;
    std::valarray<double> p (3);
    std::valarray<double> q (3);
    while (norma(grad(x)) > epslon){
        // std::cout<<"entrei no loop"<<std::endl;
        gf = grad(x);
         d = {hess[0][0]*gf[0]+hess[0][1]*gf[1]+hess[0][2]*gf[2],
              hess[1][0]*gf[0]+hess[1][1]*gf[1]+hess[1][2]*gf[2],
              hess[2][0]*gf[0]+hess[2][1]*gf[1]+hess[2][2]*gf[2]};
        t = armijo(x, d, 0.8, 0.25);
        // std::cout <<"Saí do Armijo"<<std::endl;
        x0 = x;
        x = x + t*d;
        // std::cout <<"Calculei x0 e x"<<std::endl;

        //Calcular 'p' e 'q'
        p = x-x0;
        q = grad(x)-grad(x0);
        // std::cout <<"Calculei p e q"<<std::endl;

        //Atualizar a Hessiana; agrupar as operações
        div = (p*q).sum();
        v1 = multVetMatriz(q,hess);
        m1 = multVetor(p,p);
        m2 = multVetor(p,q);
        m3 = multVetor(q,p);
        m4 = multMatriz(m2,hess);
        m5 = multMatriz(hess,m3);
        // printMatriz(m2);
        // std::cout <<"Calculei matrizes"<<std::endl;
        //Calcular a equação
        aux = (1 + ((v1*q).sum())/div);
        // std::cout <<"Calculei aux"<<std::endl;
        hess = subMatriz(somaMatriz(hess, multEscalar(divMatriz(m1,div), aux)), divMatriz(somaMatriz(m4,m5),div));
        // std::cout <<"Calculei hessiana final"<<std::endl;
    }
    printVetor(x);
    return x;

}


int main(){
    double initx[] = {0, -1, 1};
    double aa[3] = {2,2,2};
    std::valarray<double> x (initx, 3);
    
    // std::valarray< std::valarray<double> > m1 (3);
    // std::valarray< std::valarray<double> > m2 (3);
    // std::valarray<double> aaa (aa, 3);
    // m1[0]=aaa;
    // m1[1]=aaa;
    // m1[2]=aaa;
    // m2[0]=aaa;
    // m2[1]=aaa;
    // m2[2]=aaa;
    // std::cout<<"===M1==="<<std::endl;
    // printMatriz(m1);
    // std::cout<<"===M2==="<<std::endl;
    // printMatriz(m2);
    // std::cout<<"========"<<std::endl;
    // printMatriz(multMatriz(m1,m2));
    // printMatriz(somaMatriz(m1,m2));
    // printMatriz(subMatriz(m1,m2));
    // printMatriz(multEscalar(m1,100));
    // printMatriz(divMatriz(m1,100));
    // printVetor(multVetMatriz(aaa,m1));


    // std::valarray< std::valarray<double> > m3 (3);
    // m3 = multMatriz(m1,m2);
    // printMatriz(m3);

    // std::valarray< std::valarray<double> > hess (3);
    // hess = hessiana(x);
    // printMatriz(hess);
    // std::valarray<double> g (3);
    // g = grad(x);
    // std::cout<< "HELLO"<<std::endl;
    // std::cout<< hess[1][1]*g[1]<<std::endl;

    qNewton(x);

return 0;
}
