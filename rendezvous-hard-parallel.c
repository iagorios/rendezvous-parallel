#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

clock_t startTime, finalTime;

//-------------------Functions-------------------
double vX(int t, double X, double gama, double vex, double vey, double A, double B, double E);
double vY(int t, double X, double gama, double vex, double vey, double A, double B);
double vZ(int t, double X, double gama, double vez, double H, double I);

double rT(double X, double gama, double Z);
double vT(double dx, double dy, double dz);

double dX(int t, double vey, double vex, double gama, double X, double A, double B, double E, double G);
double dY(int t, double vex, double vey, double gama, double X, double A, double B, double D);
double dZ(int t, double X, double gama, double vez, double G, double H, double I);

double brute_A (double y0, double xl0, double gama, double X, double vex, double vey);
double brute_B (double yl0, double gama, double X, double vex, double vey);
double brute_D (double y0, double xl0, double X, double vex);
double brute_E (double y0, double xl0, double X, double vex);
double brute_G (double x0, double yl0, double X, double vex, double vey);
double brute_H (double z0, double gama, double vex);
double brute_I(double zl0, double gama, double X, double vez);

double x=0, y=0, z=0, xl0=0, yl0=0, zl0=0;
int Alt= 220;
double w = 398600.4418/sqrt((6378.0 + 220)*(6378.0 + 220)*(6378.0 + 220));
//otimizacao ---------------------
double ww;
//--------------------------------

static int const N = 20;
double nave = 0;
/*
 * main
 */
void main(int argc, char *argv[]) {
	//otimizacao ----------------------
	ww = w*w;

	//---------------------------------
	startTime = clock();
	int Tmax = 86400;
	int NPI = atoi(argv[1]); // numero de posicoes iniciais
	FILE *arq, *out;
	char url[] = "in.dat";
	arq = fopen(url, "r");
	out = fopen("parallel-out.txt", "w");
	double var1;

	//printf("Numero de posicoes iniciais: %d\n", NPI);

	for(int np = 1; np <= NPI; np++) {
		//printf("Problema %d\n", np);
		if(arq == NULL) {
			//printf("Erro, nao foi possivel abrir o arquivo\n");
			exit(EXIT_FAILURE);
		} else {
			fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &var1, &var1, &var1, &x, &y, &z, &var1, &xl0, &yl0, &zl0, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
			//printf("%lf %lf %lf %lf %lf %lf\n", x, y, z, xl0, yl0, zl0);
		}
		#pragma omp parallel for
		for(int VeAux = 1; VeAux<=10; VeAux++) {
			//printf("Ve %d\n", VeAux);
			double Ve =VeAux;
            Ve = Ve/2;
            double vex, vey, vez;
            vex = vey = vez =Ve*Ve/3;
            #pragma omp parallel for
            for(int aux = -14; aux<=2; aux++){
            	//printf("Gama %d\n", aux);
				double gama = pow(10, aux);
				//int tid = omp_get_thread_num();
        		//printf("Hello world from omp thread %d\n", tid);
                #pragma omp parallel for firstprivate(z, x, y, zl0, xl0, yl0)
				for(int Xaux=1; Xaux<=100; Xaux++) {
					//printf("X %d\n", Xaux);
					double X = Xaux;

					//printf("P A: y%lf xl0%lf gama%lf X%lf w%lf vex%lf", y, xl0, gama, X, w, vex);

					double A = brute_A (y, xl0, gama, X, vex, vey);
					double B = brute_B (yl0, gama, X, vex, vey);
					double D = brute_D (y, xl0, X, vex);
					double E = brute_E (y, xl0, X, vex);
					double G = brute_G (x, yl0, X, vex, vey);
					double H = brute_H (z, gama, vex);
					double I = brute_I (zl0, gama, X, vez);

					//printf("\nA:%lf \nB:%lf \nD:%lf \nE:%lf \nG:%lf \nH:%lf \nI:%lf \n", A, B, D, E, G, H, I);
					nave = nave + 1;
					int ID = omp_get_thread_num();
					//printf("Thread %d: Simulando nave %.1f\n", ID, nave);
					#pragma omp parallel for
					for(int t = 0; t <= Tmax; t++) {
						//printf("t %d\n", t);
						double fx = dX(t, vey, vex, gama, X, A, B, E, G);
                        double fy = dY(t, vex, vey, gama, X, A, B, D);
                        double fz = dZ(t, X, gama, vez, G, H, I);

						double r = rT(fx, fy, fz);

						if(r <= 0.003f && r >= -0.003f) {
							double fdx =  vX(t, X, gama, vex, vey, A, B, E);
                            double fdy =  vY(t, X, gama, vex, vey, A, B);
                            double fdz =  vZ(t, X, gama, vez, H, I);

							double v = vT( fdx, fdy, fdz);
							if(v <= 0.003f && v >= -0.003f) {
								int ID = omp_get_thread_num();
								//printf("Thread %d: RANDEZVOUS %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", ID, x, y, z, xl0, yl0, zl0, X, gama, Ve, t);
								fprintf(out, "RANDEZVOUS %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", x, y, z, xl0, yl0, zl0, X, gama, Ve, t);
								//printf("R: %lf V:%lf\n", r,v);
								t = Tmax + 1;
							}
						}
						//printf(" ======== Saídas Finais ========\n");
						//printf("R: %lf V:%lf\n", r,v);
					}
				}
			}
		}
	}
	finalTime = clock();
   	double excecutionTime = (finalTime-startTime)/CLOCKS_PER_SEC;
	fprintf(out, "Tempo em segundos: %lf", excecutionTime);
    fclose(out);
}

/**
* Calcular coeficiente A do Rendezvous
* @author Weverson, Jhone, Gledson
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param y0 valor no eixo Y da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param z0 valor no eixo z da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de A
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de A
* @param w
* @param a - 0 - PotÃªncia utilizada dentro do somÃ¡torio para casos em que o indice do somÃ¡torio Ã© utilizado elevado a potÃªncias diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de A
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de A
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de A
* @returns O coeficiÃªnte A dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_A (double y0, double xl0, double gama, double X, double vex, double vey) {
    double result = 0;
    double aux;
    double sum = 0;

    result = (2*xl0)/w - 3*y0 +((2*vex)/w)*log((X+1)/X);

    //otimizacao
    double gamavey_ww = (gama*vey)/ww;
    double gama_w = gama/w;
    double vex2_w = (2*vex)/w;
    //Calculo do somatorio
    //#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
    	//aux = (1/(n*pow(X, n)))*(1/(1+((n*gama)/w)*((n*gama)/w)))*(((2*vex)/w)+((n*gama*vey)/(w*w)));
		aux = (1/(n*pow(X, n)))*(1/(1+(n*gama_w)*(n*gama_w)))*((vex2_w)+(n*gamavey_ww));
        if (n%2 == 0) {//iteraÃ§Ã£o Par
            aux = -aux;
        }
        sum += aux;
    }

    result-= sum;

    return result;
}

/**
* Calcular coeficiente B do Rendezvous
* @author Weverson, Jhone, Gledson
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param y0 valor no eixo Y da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param z0 valor no eixo z da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de B
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de B
* @param w
* @param a - 0 - PotÃªncia utilizada dentro do somÃ¡torio para casos em que o indice do somÃ¡torio Ã© utilizado elevado a potÃªncias diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de B
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de B
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de B
* @returns O coeficiÃªnte B dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_B (double yl0,  double gama, double X, double vex, double vey) {
    double result = 0;
    double sum = 0;
    double aux;

    result = yl0/w + (vey/w)*log((X+1)/X);

    //otimizacao
    double gamavex_ww = (gama*vex)/ww;
    double gama_wpow = (gama/w)*(gama/w);
    double vey_w = vey/w;
    //Calculo do somatorio
    //#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
 //     aux = (1/(n*pow(X,n)))*(1/(1+pow(((n*gama)/w),2)))*(vey/w + (n*gama*vex)/(ww));
        aux = (1/(n*pow(X,n)))*(1/(1+(n*n*gama_wpow)))*(vey_w + (n*gamavex_ww));
        if (n%2 == 0) {//iteraÃ§Ã£o Par
            aux = -aux;
        }
        sum += aux;
    }

    result+= sum;

    return result;
}


/**
* Calcular coeficiente D do Rendezvous
* @author Weverson, Jhone, Gledson
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param y0 valor no eixo Y da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param z0 valor no eixo z da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de D
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de D
* @param w
* @param a - 0 - PotÃªncia utilizada dentro do somÃ¡torio para casos em que o indice do somÃ¡torio Ã© utilizado elevado a potÃªncias diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de D
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de D
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de D
* @returns O coeficiÃªnte D dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_D (double y0, double xl0, double X, double vex) {
    double result = 0;

    result -= (vex* log((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;

    return result;
}


/**
* Calcular coeficiente E do Rendezvous
* @author Weverson, Jhone, Gledson
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param y0 valor no eixo Y da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param z0 valor no eixo z da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de E
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de E
* @param w
* @param a - 0 - PotÃªncia utilizada dentro do somÃ¡torio para casos em que o indice do somÃ¡torio Ã© utilizado elevado a potÃªncias diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de E
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de E
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de E
* @returns O coeficiÃªnte E dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_E (double y0, double xl0, double X, double vex) {
    double result = 0;

    result -=  3*vex*log((X+1)/X);
    result +=  6*w*y0 - 3*xl0;

    return result;
}

/**
* Calcular coeficiente G do Rendezvous
* @author Iago, Filipe e JoÃ£o
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de G
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de G
* @param w
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de G
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de G
* @returns o coeficiÃªnte G dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_G (double x0, double yl0, double X, double vex, double vey) {
    double result = 0;
    double sum = 0;
    double aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    //otimizacao
    double vex3 = vex*3;
    //#pragma omp parallel for reduction(+:sum) private(aux)
	for (int n = 1; n <= N; n++) {
    	aux = vex3/(n*n*pow(X,n)*w);

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
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param x0 valor no eixo X da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param y0 valor no eixo Y da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param z0 valor no eixo z da posiÃ§Ã£o relativa inicial entre o satÃ©lite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Y Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de H
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de H
* @param w
* @param a - 0 - PotÃªncia utilizada dentro do somÃ¡torio para casos em que o indice do somÃ¡torio Ã© utilizado elevado a potÃªncias diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo X a ser calculado o valor de H
* @param vey VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Y a ser calculado o valor de H
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de H
* @returns O coeficiÃªnte H dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_H (double z0, double gama, double vex) {
    double result = 0;
    double sum = 0;
    double aux;

    result = z0;
    //otimizacao
    double vexgama = vex*gama;
    double gama_wpow = (gama/w)*(gama/w);
    //Calculo do somatorio
    //#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
    	//aux = ((vex*gama)/(pow(gama,n)*(ww)))/(1+((n*gama)/w)*((n*gama)/w));
        aux = ((vexgama)/(pow(gama,n)*(ww)))/(1+(n*n*gama_wpow));
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
* @author Iago, Filipe e JoÃ£o
* @param N NÃºmero de iteraÃ§Ãµes no somatÃ³rio interno
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satÃ©lite e o detrito
* @param Gama - VariÃ¡vel fÃ­sica Gama a ser calculado o valor de I
* @param X Chi - VariÃ¡vel fÃ­sica Chi a ser calculado o valor de I
* @param w - velocidade angular
* @param vez VariÃ¡vel fÃ­sica da Velocidade de exaustÃ£o no eixo Z a ser calculado o valor de I
* @returns o coeficiÃªnte I dado os valores iniciais e as variÃ¡veis fÃ­sicas a serem testadas
*/
double brute_I(double zl0, double gama, double X, double vez) {
    double result = 0;
    double sum = 0;
    double aux;

    result = zl0/w - (vez/w)*(log((X+1)/X));

    //otimizacao
    double gama_wpow = (gama/w)*(gama/w);

    //Calculo do somatorio
    //#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
        aux = ((vez)/(n*n*pow(X,n)*w))/(1+(n*n*gama_wpow));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;

    return result;
}

/* @author Weverson, Iago
 * vetor X da distancia
 */
double dX(int t, double vey, double vex, double gama, double X, double A, double B, double E, double G) {
	//otimizacao
	double wt = w*t;

	double resultFn = 0;
	double result1 = 2 * (A*sin(wt)-B*cos(wt))+E*t;
	double result2 = G;
	//otimizacao
	double vey2_w = (2*vey)/w;
	double vex4 = vex*4;
	double gama_wpow = (gama/w)*(gama/w);
	double gamat = gama*t;

	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n <= N; n++) {
		// brute_F
		resultFn = (1/(n*pow(X,n)))*(vey2_w + vex4/(n*gama))/(1+(n*n*gama_wpow));

		if (n%2 == 0) {
            resultFn = - resultFn;
        }
		resultFn -= vex/(n*gama);
		//brute_F

		result2 += resultFn * pow(M_E, -(n * gamat));
	}
	return result1 + result2;
}

double dY (int t, double vex, double vey, double gama, double X, double A, double B, double D) {
	//otimizacao
	double wt = w*t;

	double resultCn = 0;
	double result1 = A*cos(wt)+B*sin(wt);
	double result2 = D;

	//otimizacao
	double gama_wpow = (gama/w)*(gama/w);
    double gamavey_ww = (gama*vey)/ww;
	double gamat = gama*t;
	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n < N; ++n){
		//brute_C
		resultCn = 1/(n*pow(X,n))*(vex+(n*gamavey_ww))*(1/(1+(n*n*gama_wpow)));

		if (n%2 == 0) {
	        resultCn = -resultCn;
	    }
		//brute_C

		result2 += resultCn*pow(M_E, -(n*gamat));
	}
	return result1 + result2;
}

/* @author Weverson, Iago
 * vetor Z da distancia
 */
double dZ(int t, double X, double gama, double vez, double G, double H, double I) {
	//otimizacao
	double wt = w*t;

	double resultJn = 0;
	double result1 = H * cos(wt) + I * sin(wt);
	double result2 = 0;

	//otimizacao
	double gama_wpow = (gama/w)*(gama/w);
	double gamat = gama*t;
	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n <= N; n++) {
		//brute_J
		resultJn = vez/(n*pow(X,n)*w)/(1+(n*n*gama_wpow));
		if (n%2 == 0) {
	        resultJn = -resultJn;
	    }
		//brute_J

		result2 += resultJn * pow(M_E, -(n * gamat));
	}
	return result1 - result2;
}

/* @author Weverson, Jhone, Gledson
 * vetor X da velocidade
 */
double vX(int t, double X, double gama, double vex, double vey, double A, double B, double E) {
	//otimizacao
	double wt = w*t;

	double resultFn = 0;
	double result1 = 2 * ((A*w*sin(wt)) + (B*w*cos(wt)));
	double result2 = E;

	//otimizacao
	double vey2_w = (2*vey)/w;
	double vex4 = vex*4;
	double gama_wpow = (gama/w)*(gama/w);
	double gamat = gama*t;
	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n <= N; n++) {
		// brute_F
		//resultFn = (1/(n*pow(X,n)))*((2*vey)/w + (4*vex)/(n*gama))/(1+pow((n*gama)/w,2));
		resultFn = (1/(n*pow(X,n)))*(vey2_w + vex4/(n*gama))/(1+(n*n*gama_wpow));

		if (n%2 == 0) {
	        resultFn = - resultFn;
	    }
		resultFn -= vex/(n*gama);
		//brute_F

		result2 += resultFn*((-n)*gama*pow(M_E, -(n*gamat)));
	}
	return result1 + result2;
}

/* @author Weverson, Jhone, Gledson
 * vetor Y da velocidade
 */
double vY(int t, double X, double gama, double vex, double vey, double A, double B) {
	//otimizacao
	double wt = w*t;

	double resultCn = 0;
	double result1 = (-A)*w*sin(wt)+B*w*cos(wt);
	double result2 = 0;

	//otimizacao
    double gamavey_ww = (gama*vey)/ww;
	double gama_wpow = (gama/w)*(gama/w);
	double gamat = gama*t;
	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n <= N; n++) {
		//brute_C
		//resultCn = 1/(n*pow(X, n))*(vex+(n*gama*vey)/(ww))*(1/(1+(n*gama/w)*(n*gama/w)));
		resultCn = 1/(n*pow(X, n))*(vex+(n*gamavey_ww))*(1/(1+(n*n*gama_wpow)));

		if (n%2 == 0) {
	        resultCn = -resultCn;
	    }
		//brute_C

		result2 += resultCn*((-n)*gama*pow(M_E, -(n*gamat)));
	}
	return result1 + result2 ;
}

/* @author Weverson
 * vetor Z da velocidade
 */
double vZ(int t, double X, double gama,  double vez, double H, double I) {
	//otimizacao
	double wt = w*t;

	double resultJn = 0;
	double result1 = (-H)*w*sin(wt)+I*w*cos(wt);
	double result2 = 0;

	//otimizacao
	double gama_wpow = (gama/w)*(gama/w);
	double gamat = gama*t;
	//#pragma omp parallel for reduction(+:result2)
	for (int n = 1; n <= N; n++) {
		//brute_J
		resultJn = vez/(n*pow(X,n)*w)/(1+(n*n*gama_wpow));
		if (n%2 == 0) {
	        resultJn = -resultJn;
	    }
		//brute_J

		result2 += resultJn*((-n)*gama*pow(M_E, -(n*gamat)));
	}

	return result1  - result2;
}

/* @author Weverson
 * Modificada por Gledson
 * Velocidade
 */
double rT(double X, double gama, double Z) {
	return sqrt(X*X + gama*gama + Z*Z);
}

/* @author Weverson
 * Modificada por Gledson
 * Posicao
 */
double vT(double dx, double dy, double dz) {
	return sqrt(dx*dx + dy*dy + dz*dz);
}
