#include <valarray>

std::valarray< std::valarray<double> > multVetor(std::valarray<double> V1, std::valarray<double> V2)
{
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

std::valarray<double> multVetMatriz(std::valarray<double> V1, std::valarray< std::valarray<double> > M1, int n =3)
{
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


std::valarray<double> multMatVet(std::valarray< std::valarray<double> > M1, std::valarray<double> V1, int n =3)
{
	std::valarray<double> V (n);

	for (int i=0; i<n; i++){
		V[i] = (M1[i]*V1).sum();
	}

	return V;
}


std::valarray< std::valarray<double> > multMatriz(std::valarray< std::valarray<double> > M1, std::valarray< std::valarray<double> > M2, int n =3)
{
    //Multipica duas matrizes M1 e M2
    std::valarray< std::valarray<double> > M (n);
    std::valarray<double> zero (0.0,n);

    for (int i=0; i<n; i++)
    {
        M[i] = zero;
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                M[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }

    return M;
}

std::valarray< std::valarray<double> > divMatriz(std::valarray< std::valarray<double> > M1, double divisor=1, int n=3)
{
    for (int i=0; i<n; i++){
        M1[i] = M1[i]/divisor;
    }

    return M1;
}

std::valarray< std::valarray<double> > somaMatriz(std::valarray< std::valarray<double> > M1, std::valarray< std::valarray<double> > M2, int n=3)
{
    std::valarray< std::valarray<double> > M (3);

    for (int i=0; i<n; i++){
        M[i] = M1[i]+M2[i];
    }

    return M;
}

std::valarray< std::valarray<double> > subMatriz(std::valarray< std::valarray<double> > M1, std::valarray< std::valarray<double> > M2, int n=3)
{
    std::valarray< std::valarray<double> > M (3);

    for (int i=0; i<n; i++){
        M[i] = M1[i]-M2[i];
    }

    return M;
}

std::valarray< std::valarray<double> > multEscalar(std::valarray< std::valarray<double> > M1, double a, int n=3)
{
    for (int i=0; i<n; i++){
        M1[i] = M1[i]*a;
    }

    return M1;
}

std::valarray< std::valarray<double> > getIdentidade(int n = 3)
{
	std::valarray< std::valarray<double> > M (n);
	std::valarray<double> l (0.0,n);

	for (int i=0; i<n; i++){
		M[i] = l;
		M[i][i] = 1.0;
	}

	return M;
}

void printMatriz(std::valarray< std::valarray<double> > M, int n=3)
{
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){ 
        std::cout << M[i][j] <<" ";
        }
        std::cout <<""<<std::endl;
    }
}

void printVetor(std::valarray <double> V, int n=3)
{
    std::cout<<"[ ";

    for (int i=0; i<n; i++){
        std::cout << V[i] << " ";
    }

    std::cout << "]"<<std::endl;
}
