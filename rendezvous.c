#include <stdio.h>
#include <math.h>

//-------------------Functions-------------------
double vX(double A, double w, double t, double B, double E, double F, int N, int gama);
double vY(double H, double w, double t, double I, double J, int N, int gama);
double vZ(double H, double w, double t, double I, double J, int N, int gama);
double A, B, C, D, E, F, G, H, I, J;

//#define int N = 20;

void main() {
	printf("Resultado: %f", vY(10, 5, 2, 1, 1, 20, 1));
}

/**
* Calcular coeficiente A do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de A
* @param X Chi - Variável física Chi a ser calculado o valor de A
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de A
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de A
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de A
* @returns O coeficiênte A dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_A (int N, double y0, double xl0, double gama, double X, double w, double vex, double vey) {
    double result = 0;
    int n;
    double aux;
    double sum = 0;

    result += (2*xl0)/w - 3*y0 +((2*vex)/w)*log((X+1)/X);

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
//		aux = (1/(n*pow(X, n)))*(1/(1+pow(((n*Y)/w),2)))*(((2*vex)/w)+((n*Y*vey)/(w*w)));
		aux = (1/(n*pow(X, n)))*(1/(1+((n*gama)/w)*((n*gama)/w)))*(((2*vex)/w)+((n*gama*vey)/(w*w)));
        if (n%2 == 0) {//iteração Par
            aux = -aux;
        }
        sum += aux;
    }

    result+= sum;

    return result;
}

/**
* Calcular coeficiente B do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de B
* @param X Chi - Variável física Chi a ser calculado o valor de B
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de B
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de B
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de B
* @returns O coeficiênte B dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_B (int N, double yl0, double gama, double X, double w, double vex, double vey){
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result += yl0/w + (vey/w)*log((X+1)/X);

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
//      aux = (1/(n*pow(X,n)))*(1/(1+pow(((n*Y)/w),2)))*(vey/w + (2*n*Y*vex)/(w*w));
        aux = (1/(n*pow(X,n)))*(1/(1+pow(((n*gama)/w),2)))*(vey/w + (n*gama*vex)/(w*w));
        if (n%2 == 0) {//iteração Par
            aux = -aux;
        }
        sum += aux;
    }

    result+= sum;

    return result;
}

/**
* Calcular o somatório dos coeficientes Cn do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de C
* @param X Chi - Variável física Chi a ser calculado o valor de C
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de C
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de C
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de C
* @returns O somatório dos coeficiêntes Cn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_C (int n, double gama, double X, double w, double vex, double vey){
    double result = 0;
    double aux;

    //Calculo do somatorio Cn
    aux = 1/(n*pow(X, n)) * (vex + (n * gama * vey)/(w*w)) * (1/(1 + (n*gama/w)*(n*gama/w) ));
    if (n%2 == 0) {
        aux = -aux;
    }

    result +=aux;

    return result;
}

/**
* Calcular coeficiente D do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de D
* @param X Chi - Variável física Chi a ser calculado o valor de D
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de D
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de D
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de D
* @returns O coeficiênte D dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_D (int N, double y0, double xl0, double Y, double X, double w, double vex) {
    double result = 0;

    result -= (vex* log((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;

    return result;
}


/**
* Calcular coeficiente E do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de E
* @param X Chi - Variável física Chi a ser calculado o valor de E
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de E
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de E
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de E
* @returns O coeficiênte E dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_E (int N, double y0, double xl0, double X, double w, double vex) {
    double result = 0;

    result -=  3*vex*log((X+1)/X);
    result +=  6*w*y0 - 3*xl0;

    return result;
}

/**
* Calcular o somatório dos coeficientes Fn do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de Fn
* @param X Chi - Variável física Chi a ser calculado o valor de Fn
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de Fn
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de Fn
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Fn
* @returns O somatório coeficiênte Fn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;


    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*pow(X,n)))*((2*vey)/w + (4*vex)/(n*Y))/((1+pow((n*Y)/w,2)));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux -= vex/(n*Y);
        aux *= pow(n,a);
        sum += aux;
    }

    result = sum;

    return result;
}

/**
* Calcular coeficiente G do Rendezvous
* @author Iago, Filipe e João
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de G
* @param X Chi - Variável física Chi a ser calculado o valor de G
* @param w
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de G
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de G
* @returns o coeficiênte G dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_G (int N, double x0, double yl0, double X, double w, double vex, double vey) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;

		for (n = 1; n <= N; n++) {
        aux = 3*vex/(pow(n,2)*pow(X,n)*w);

				if (n%2 == 0) {
            aux = -aux;
        }

        sum +=aux;
    }

    result-=sum;

    return result;
}

/**
* Authors: Gledson e Jhone
* Calcular coeficiente H do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de H
* @param X Chi - Variável física Chi a ser calculado o valor de H
* @param w
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de H
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de H
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de H
* @returns O coeficiênte H dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_H (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result = z0;
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vex*Y)/(pow(Y,n)*(w*w)))/(1+((n*Y)/w)*((n*Y)/w));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }
    result += sum;

    return result;
}

/**
* Calcular coeficiente I do Rendezvous
* @author Iago, Filipe e João
* @param N Número de iterações no somatório interno
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de I
* @param X Chi - Variável física Chi a ser calculado o valor de I
* @param w - velocidade angular
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de I
* @returns o coeficiênte I dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_I (int N, double zl0, double Y, double X, double w, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result = zl0/w - (vez/w)*(log((X+1)/X));

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez)/(pow(n,2)*pow(X,n)*w))/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;

    return result;
}

/**
* Calcular o somatório dos coeficientes Jn do Rendezvous
* @author Iago, Filipe e João
* @param Y Gama - Variável física Gama a ser calculado o valor de Jn
* @param X Chi - Variável física Chi a ser calculado o valor de Jn
* @param w - velocidade angular
* @param n - indíce do somatório
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Jn
* @returns o somatório coeficiênte Jn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_J(double Y, double X, double w, double vez, int n){
    double result = 0;

    result = vez/(n*pow(X,n)*w)/(1+pow((n*Y)/w,2));

    if (n%2 == 0) {
        result = -result;
    }

    return result;
}

// vetor X da velocidade
double vX(double A, double w, double t, double B, double E, double F, int N, int gama) {
	double result1 = 2 * ( (A * w * cos(w * t)) + (B * w * sin(w * t)) ) + E;
	double result2 = 0;
	for (int n = 1; n <= N; n++) {
		result2 += F * ((-n) * gama * pow(M_E, -(n * gama * t)));
	}
	return result1 + result2;
}

// vetor Y da velocidade
double vY(double A, double w, double t, double B, double C, int N, int gama) {
	double result1 = (-A) * w * sin(w * t);
	double result2 = B * w * cos(w * t);
	double result3 = 0;
	for (int n = 1; n <= N; n++) {
		result3 += C * ((-n) * gama * pow(M_E, -(n * gama * t)));
	}
	return result1 + result2 + result3;
}

// vetor Z da velocidade
double vZ(double H, double w, double t, double I, double J, int N, int gama) {
	double result1 = (-H) * w * sin(w * t);
	double result2 = I * w * cos(w * t);
	double result3 = 0;

	for (int n = 1; n <= N; n++) {
		result3 += J * ((-n) * gama * pow(M_E, -(n * gama * t)));
	}

	return result1 + result2 + result3;
}

double dY (double A, double B, double C, double D, double Y, int N, double w, double t) {
	double result1 = A*cos(w*t)+B*sin(w*t);
	double result2 = 0;
	for (int n = 1; n < N; ++n){
		result2 = C*pow(M_E, -(n*w*t)) + D;
	}
	return result1 + result2;
}
