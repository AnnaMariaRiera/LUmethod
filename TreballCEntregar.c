//Anna Maria Riera Escandell, Física-Matemàtiques. NIU: 1359293

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TOL 0.000001


void PintaMatriu ( float** matriu, int n)
{
	int i, j;
	printf("\n");
	for (i=0 ; i<n ; i++)
	{ 
		for ( j=0 ; j<n ; j++)
		{
			printf("%f ", matriu[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

int LlegeixMatriuISolucions ( float *** m, float ** v )
/*Posam float *** m perquè és un punter a una matriu, i una matriu, al seu torn, és un punter de punters*/
{
	int i, j, n;
	float ** matriu;
	float *  vector;
	
	printf ("\n\nDigues la mida de la matriu:\n");
	scanf("%d", &n);
	
	/*Reservem memòria*/
	matriu = (float**) malloc (n*sizeof(float*));
	vector = (float*) malloc (n*sizeof(float));
	
	for (i=0 ; i<n ; i++)
	{ 
		matriu[i] = (float*) malloc (n*sizeof(float));
		
		for ( j=0 ; j<n ; j++)
		{
			printf("Dona el valor del coeficient (%d, %d):\n", i+1, j+1);
			scanf("%f", &matriu[i][j]);
		}
		
		printf("Dona l'element %d de la solucio:\n", i+1);
		scanf("%f", &vector[i]);
	}
	
	*m = matriu;
	*v = vector;
	
return n;
} 

void Pivotar ( float** matriu, float** P, float** Q, int on_som, int n)
/*El paràmetre "on_som" és perquè no voldrem que ens busqui l'element (en valor absolut) més gran de tota la matriu, sinó que voldrem que
 ens el busqui sense mirar les files anteriors. */
{
	int i, j, max_i, max_j;
	float intermig;

	max_i = on_som;
	max_j = on_som;
	
	for (i=on_som; i<n; i++) {
		for (j=on_som; j<n ; j++) {
			if ( abs( matriu[i][j] ) > abs( matriu[max_i][max_j] ))  {
				max_i=i;    
				max_j=j;
			}
		}
	}
	
	/*Volem pivotar no des de la primera fila, sinó des d'on ens interessi: des de "on_som"*/
	for (i=0; i<n; i++) {
		intermig = matriu[on_som][i];
		matriu[on_som][i] = matriu[max_i][i];
		matriu[max_i][i] = intermig;
		
		intermig = P[on_som][i];
		P[on_som][i] = P[max_i][i];
		P[max_i][i] = intermig;
	}
	/*Com al joc del solitari, necessitem una variable intermitja que ens ajudi a canviar de lloc les files (i les columnes) */
	
	for (j=0; j<n; j++) {
		intermig = matriu[j][on_som];
		matriu[j][on_som] = matriu[j][max_j];
		matriu[j][max_j] = intermig;
		
		intermig = Q[j][on_som];
		Q[j][on_som] = Q[j][max_j];
		Q[j][max_j] = intermig;
	}
}			

float** LU_SensePivotatge ( float** matriu, float*** lu, int n)
{
	int i, j, k;
	float multiplicador;
	float **LU, **L;
	
	//Reservem memòria per a una fila
	LU = (float**) malloc (n*sizeof(float*));
	L = (float**) malloc (n*sizeof(float*));
	
	for (i = 0; i < n; i++){   //Demanem memòria per a cada columna de cada element de la fila
		LU[i] = (float*) malloc (n*sizeof(float));
		L[i] = (float*) malloc (n*sizeof(float));
		for (j = 0; j < n; j++){
			LU[i][j] = matriu[i][j];
		}
	}
	

	for (i=0 ; i<n ; i++)
	{
		if ( fabs(LU[i][i])<TOL) {
			return NULL;
		}
		for (j=i+1 ; j<n ; j++)
		{
			multiplicador = LU[j][i] / LU[i][i];
			for (k=0 ; k<n ; k++)
			{
				LU[j][k]= LU[j][k] - (LU[i][k]*multiplicador);
			}
			L[j][i]=multiplicador;
		}
	}
	
	/*Imprimirem L i U per separat per a què l'usuari les pugui veure i emplenem la part inferior de LU*/
	//Emplenem LU
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<i ; j++){
			LU[i][j] = L[i][j];
		}
	}
	
	//Imprimim U
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = LU[i][j];
			if (j<i) {
				L[i][j] = 0;
			}
		}
	}
	printf("La matriu U es:");
	PintaMatriu(L,n);
	
	//Imprimim L
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = 0;
			if (j==i){
				L[i][i] = 1 ;
			}
			else if (j<i) {
				L[i][j] = LU[i][j];
			}
		}
	}
	printf("\n La matriu L es:");
	PintaMatriu(L, n);
	
	//Alliberem la memòria de L
	for (i=0; i<n; i++){
		free(L[i]);
		L[i] = NULL;
	}
	free(L);
	L = NULL;
	
	*lu = LU;
	
return (*lu);
}

float** CreaIdentitat (int n)
{
	int i, j;
	float** Id;
	
	Id = (float**) malloc (n*sizeof(float*));
	
	for (i=0; i < n; i++){
		Id[i] = (float*) malloc (n*sizeof(float));
		
		for (j = 0; j < n; j++){
			if ( i==j ) {
				Id[i][j] = 1;
			}
			else {
				Id[i][j]=0;
			}
		}
	}
	
return(Id);
}

float** LU_AmbPivotatge ( float** matriu, float*** lu, float** P, float** Q, int n)
{
	int i, j, k;
	float multiplicador;
	float **LU, **L;
	
	LU = (float**) malloc (n*sizeof(float*));
	L = (float**) malloc (n*sizeof(float*));
	
	for (i = 0; i < n; i++){
		LU[i] = (float*) malloc (n*sizeof(float));
		L[i] = (float*) malloc (n*sizeof(float));
		for (j = 0; j < n; j++){
			LU[i][j] = matriu[i][j];
		}
	}
	
	for (i=0 ; i<n ; i++)
	{
		Pivotar ( LU, P, Q, i, n );
		if ( fabs(LU[i][i])<TOL ) {
			printf ("\nEl determinant de la matriu es nul o algun pivot es menor que %f.\nEl metode no es pot executar.\n", TOL);
			return NULL;
		}
		for (j=i+1 ; j<n ; j++)
		{
			multiplicador = LU[j][i] / LU[i][i];
			for (k=0 ; k<n ; k++)
			{
				LU[j][k]= LU[j][k] - (LU[i][k]*multiplicador);
			}
			L[j][i]=multiplicador;
		}
	}
	
	//Emplenem LU
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<i ; j++){
			LU[i][j] = L[i][j];
		}
	}
	
	//Imprimim U
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = LU[i][j];
			if (j<i) {
				L[i][j] = 0;
			}
		}
	}
	printf("La matriu U es:");
	PintaMatriu(L,n);
	
	//Imprimim L
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = 0;
			if (j==i){
				L[i][i] = 1 ;
			}
			else if (j<i) {
				L[i][j] = LU[i][j];
			}
		}
	}
	printf("\n La matriu L es:");
	PintaMatriu(L, n);
	
	*lu = LU;
	
return (*lu);
}

float* TriangularInferior (float** LU, float* vector, int n)
{
	float *y;
	float sumatori; 
	int i, j;
	
	y = (float*) malloc (n*sizeof(float));
	
	y[0] = vector[0];
	for ( i=1; i<n; i++ ) {
		sumatori=0;
		for ( j=0; j<i; j++ ) {
			sumatori += LU[i][j]*y[j];
		}
		y[i] = vector[i] - sumatori;
	}
	
return y;
}

float* TriangularSuperior (float** LU, float* y, int n)
{
	float *x;
	float sumatori; 
	int i, j;
	
	x = (float*) malloc (n*sizeof(float));
	
	x[n-1] = y[n-1] / LU[n-1][n-1];  /*Es conta de 0 fins a n-1*/
	for (i=n-2; i>=0; i--) {
		sumatori=0;
		for (j=n-1; j>=i; j--) {
			sumatori += LU[i][j]*x[j];
		}
		x[i] = ( y[i] - sumatori ) / LU[i][i]; 
	}
	
return x;
}

float* MatriuPerVector (float** m, float* v, int n)
{
	int i, j;
	float* VectFinal;
	
	VectFinal = (float*) malloc (n*sizeof(float));
	
	for (i=0; i<n; i++) {
		VectFinal[i] = 0;
		for (j=0; j<n; j++) {
			VectFinal[i] += m[i][j]*v[j];
		}
	}
	
return VectFinal;
}

float* MetodeComplet_SensePivotatge (float** matriu, float* vector, int n)
{
	int i;
	float** LU; 
	float *y, *resultat;
	
	LU = LU_SensePivotatge (matriu, &LU, n);
	if (LU == NULL){
		return NULL;
	}
	
	y = TriangularInferior (LU, vector, n);
	resultat = TriangularSuperior (LU, y, n);
	
	//Alliberem memòria:
	for (i=0; i<n; i++){
		free(LU[i]);
		LU[i] = NULL;
	}
	free(LU);
	LU = NULL;
	
	free(y);
	y = NULL;
	
return resultat;
}

float* MetodeComplet_AmbPivotatge (float** matriu, float* vector, int n)
{
	int i;
	float **LU, **P, **Q; 
	float *y, *resultat2, *vector2, *resultat;
	
	P = CreaIdentitat (n);
	Q = CreaIdentitat (n);
	LU = LU_AmbPivotatge( matriu, &LU, P, Q, n);  /* Li passem l'adreça de LU perquè la passem per referència: a la funció és *** */
	
	if (LU == NULL){
		return NULL;
	}
	
	vector2 = MatriuPerVector ( P, vector, n);
	
	y = TriangularInferior (LU, vector2, n);
	resultat2 = TriangularSuperior (LU, y, n);
	
	resultat = MatriuPerVector (Q, resultat2, n);
	
	//Alliberem memòria:
	for (i=0; i<n; i++){
		free(LU[i]);
		LU[i] = NULL;
	}
	free(LU);
	LU = NULL;
	
	for (i=0; i<n; i++){
		free(P[i]);
		P[i] = NULL;
	}
	free(P);
	P = NULL;
	
	for (i=0; i<n; i++){
		free(Q[i]);
		Q[i] = NULL;
	}
	free(Q);
	Q = NULL;
	
	free(y);
	y = NULL;
	
	free(resultat2);
	resultat2 = NULL;
	
	free(vector2);
	vector2 = NULL;
	
return resultat;
}

float* PivotatgeCondicional (float** matriu, float* vector, int n)
{
	int resposta;
	float *resultat;
	
	printf ("\nEl determinant de la matriu es nul o algun pivot es menor que %f.\nProva amb pivotatge:\n", TOL);
	
	printf("\nQue vols fer?\n 1. Vull provar amb pivotatge.\n 2. Vull sortir del programa.\n");
	scanf("%d", &resposta);

	if (resposta==1) {
		resultat = MetodeComplet_AmbPivotatge (matriu, vector, n);
		if (resultat == NULL) {
			return NULL;
		}
	}
	else if (resposta==2) {
		return NULL;
	}
	else {
		printf ("Has triat una opcio que no es valida\n");
		return 0;
	}
	
return resultat;
}

int main()
{
	int resposta, n, i;
	float** matriu;
	float *vector, *resultat;
	
	n = LlegeixMatriuISolucions( &matriu, &vector);
	
	printf("\nVols trobar el resultat utilitzant el metode amb pivotatge?\n 1. Si\n 2. No\n");
	scanf("%d", &resposta);
	
	if (resposta==1) {
		resultat = MetodeComplet_AmbPivotatge (matriu, vector, n);
		if (resultat == NULL) {
			return 1;
		}
	}
	else if (resposta==2) {
		resultat = MetodeComplet_SensePivotatge (matriu, vector, n);
		if (resultat == NULL) {
			resultat = PivotatgeCondicional (matriu, vector, n);
			if (resultat == NULL) {
			return 1;
			}
		}
	}
	else {
	printf ("Has triat una opcio que no es valida");
	return 0;
	}
	
	printf ("\n\nResultat:\n");
	for (i=0; i<n; i++) {
		printf("\tx%d = %f\n", i+1, resultat[i]);
	}
	
	//Alliberem memòria:
	for (i=0; i<n; i++){
		free(matriu[i]);
		matriu[i] = NULL;
	}
	free(matriu);
	matriu = NULL;
	
	free(vector);
	vector = NULL;
	
	free(resultat);
	resultat = NULL;
	
return 0;
}