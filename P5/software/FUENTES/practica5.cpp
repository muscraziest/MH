#include "arff_parser.h"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <math.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
// Funciones para generar números pseudoaleatorios
/////////////////////////////////////////////////////////////////////////////////

unsigned long Seed;

#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9

/* Inicializa la semilla al valor x.
   Solo debe llamarse a esta funcion una vez en todo el programa */
void Set_random (unsigned long x) {

    Seed = (unsigned long) x;
}

/* Genera un numero aleatorio real en el intervalo [0,1[
   (incluyendo el 0 pero sin incluir el 1) */
float Rand(void) {

    return (( Seed = ( (Seed * PRIME) & MASK) ) * SCALE );
}

/* Genera un numero aleatorio entero en {low,...,high} */
int Randint(int low, int high){

    return (int) (low + (high-(low)+1) * Rand());
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular la distancia euclídea entre dos vectores
/////////////////////////////////////////////////////////////////////////////////
template<class Iter_T, class Iter2_T>
double distanciaVectores(Iter_T first, Iter_T last, Iter2_T first2) {
  double ret = 0.0;
  while (first != last) {
    double dist = (*first++) - (*first2++);
    ret += dist * dist;
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
}


/////////////////////////////////////////////////////////////////////////////////
// Lector ARFF
/////////////////////////////////////////////////////////////////////////////////
void lectorARFF(int idArchivo,	vector<vector <double> > &v,int &num_filas,int &num_atributos){
	string nombreArchivo;
	switch(idArchivo){
	    case 1:			
	    	nombreArchivo = "./datos/libras00.arff";		    
	    	break;

	    case 2:			
	    	nombreArchivo = "./datos/libras01.arff";			
	    	break;

	    case 3:			
	    	nombreArchivo = "./datos/libras02.arff";			
	    	break;

	    case 4:			
	    	nombreArchivo = "./datos/libras03.arff";			
	    	break;
	    
	    case 5:			
	    	nombreArchivo = "./datos/libras04.arff";			
	    	break;
	    
	    case 6:			
	    	nombreArchivo = "./datos/libras05.arff";			
	    	break;
	    
	    case 7:			
	    	nombreArchivo = "./datos/libras06.arff";			
	    	break;
	    
	    case 8:			
	    	nombreArchivo = "./datos/libras07.arff";			
	    	break;
	    
	    case 9:			
	    	nombreArchivo = "./datos/libras08.arff";			
	    	break;
	    
	    case 10:			
	    	nombreArchivo = "./datos/libras09.arff";			
	    	break;

	    case 11:			
	    	nombreArchivo = "./datos/arrhythmia00.arff";		    
	    	break;
	    
	    case 12:			
	    	nombreArchivo = "./datos/arrhythmia01.arff";			
	    	break;
	    
	    case 13:			
	    	nombreArchivo = "./datos/arrhythmia02.arff";			
	    	break;
	    
	    case 14:			
	    	nombreArchivo = "./datos/arrhythmia03.arff";			
	    	break;
	    
	    case 15:			
	   		nombreArchivo = "./datos/arrhythmia04.arff";			
	    	break;
	    
	    case 16:			
	    	nombreArchivo = "./datos/arrhythmia05.arff";			
	    	break;
	    
	    case 17:			
	    	nombreArchivo = "./datos/arrhythmia06.arff";			
	    	break;
	    
	    case 18:			
	    	nombreArchivo = "./datos/arrhythmia07.arff";			
	    	break;
	    
	    case 19:			
	    	nombreArchivo = "./datos/arrhythmia08.arff";			
	    	break;
	    
	    case 20:			
	    	nombreArchivo = "./datos/arrhythmia09.arff";			
	    	break;

	    case 21:			
	    	nombreArchivo = "./datos/wdbc00.arff";		    
	    	break;
	    
	    case 22:			
	    	nombreArchivo = "./datos/wdbc01.arff";			
	    	break;
	    
	    case 23:			
	    	nombreArchivo = "./datos/wdbc02.arff";			
	    	break;
	    
	    case 24:			
	    	nombreArchivo = "./datos/wdbc03.arff";			
	    	break;
	    
	    case 25:			
	    	nombreArchivo = "./datos/wdbc04.arff";			
	    	break;
	    
	    case 26:			
	    	nombreArchivo = "./datos/wdbc05.arff";			
	    	break;
	    
	    case 27:			
	    	nombreArchivo = "./datos/wdbc06.arff";			
	    	break;
	    
	    case 28:			
	    	nombreArchivo = "./datos/wdbc07.arff";			
	    	break;
	    
	    case 29:			
	    	nombreArchivo = "./datos/wdbc08.arff";			
	    	break;
	    
	    case 30:			
	    	nombreArchivo = "./datos/wdbc09.arff";			
	    	break;
	}
	
	ArffParser parser(nombreArchivo);
	ArffData *data = parser.parse();
	
	num_filas = data->num_instances();
	num_atributos = data->num_attributes();

	ArffValue *valor;
	ArffInstance *instancia;

	for(int j=0;j<num_filas;++j){
		instancia = data->get_instance(j);
		vector<double> aux;

		for (int i = 0; i < data->num_attributes(); ++i) {
			aux.push_back(stod(instancia->get(i)->operator string()) );
		}
		v.push_back(aux);
		aux.clear();
	}

}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular el mínimo y el máximo de cada columna 
/////////////////////////////////////////////////////////////////////////////////
void calcularMinMax(vector<vector <double> > &v,vector<double> &minimoColumnas,vector<double>  &maximoColumnas,int num_filas,int num_atributos){
	maximoColumnas = minimoColumnas = v[0];
	for(int i=0; i < num_filas; ++i){
		for(int j=0; j < num_atributos; ++j){
			if(v[i][j] > maximoColumnas[j])
				maximoColumnas[j] = v[i][j];
			if(v[i][j] < minimoColumnas[j])
				minimoColumnas[j] = v[i][j];
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////
// Función para normalizar el 
/////////////////////////////////////////////////////////////////////////////////
void normalizar(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int num_filas,int num_atributos,vector<double> &minimoColumnas,vector<double>  &maximoColumnas){
	
	vector<double> aux;
	for(int i=0; i < num_filas; ++i){
		for(int j=0; j < num_atributos-1; ++j){
			if(isnan((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j])))
				aux.push_back(0);
			else
			aux.push_back((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j]));
		}

		vNormalizado.push_back(aux);
		aux.clear();
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para aplicar la máscara a los datos 
/////////////////////////////////////////////////////////////////////////////////
void aplicarMascara(vector<vector <double> > &vNormalizado,vector<vector <double> > &vConMascara,vector<int> &mascara,int num_filas,int num_atributos){
	
	vector<double> aux;

	for(int i=0; i < num_filas; ++i){
		for(int j=0; j < num_atributos; ++j){
			aux.push_back(vNormalizado[i][j]*mascara[j]);
		}

		vConMascara.push_back(aux);
		aux.clear();
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular la distancia euclidea entre todas las filas
/////////////////////////////////////////////////////////////////////////////////
void distanciaEuclidea(vector<double> &ejemploElegido,vector<vector <double> > &vConMascara,vector<double> &vectorDistancias,int num_filas,int num_atributos){
	
	vector<double> aux;
	for(int i=0; i < num_filas; ++i){
		aux = vConMascara[i];
		vectorDistancias.push_back(distanciaVectores(aux.begin(), aux.end(), ejemploElegido.begin()));
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular los 3 valores más cercanos (de menor distancia)
/////////////////////////////////////////////////////////////////////////////////
void calcularMenores(vector<double> &vectorDistancias,int &indice1, int &indice2, int &indice3){
	
	indice1 = 0;

	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice1])
			indice1 = i;		
	}

	if(indice1 != 0)
		indice2 = 0;

	else
		indice2 = 1;


	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice2] && i != indice1)
			indice2 = i;		
	}

	if(indice1 != 0 && indice2 != 0)
		indice3 = 0;
	
	else{

		if(indice1 != vectorDistancias.size()-1 && indice2 != vectorDistancias.size()-1)
			indice3 = vectorDistancias.size()-1;
		else if(indice1 == vectorDistancias.size()-1)
			indice3 = indice2+1;
		else
			indice3 = indice1+1;

	}

	
	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice3] && i != indice1 && i != indice2)
			indice3 = i;		
	}
	
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular el índice de máximo valor
/////////////////////////////////////////////////////////////////////////////////
int indiceMax(vector<int> &resultado){

	int indice = 0;
	for(int i=0; i < resultado.size(); ++i)
		if(resultado[i] > resultado[indice])
			indice = i;
	return indice;
}

/////////////////////////////////////////////////////////////////////////////////
// Clasificador knn
/////////////////////////////////////////////////////////////////////////////////
int knn(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int ejemplo, vector<int> &mascara,int num_filas,int num_atributos){
	
	//Aplicamos la máscara
	vector<vector <double> > vConMascara;
	aplicarMascara(vNormalizado,vConMascara,mascara,num_filas,num_atributos);

	vector<double> ejemploElegido = vConMascara[ejemplo];

	//Calculamos las distancias 
	vector<double> vectorDistancias;
	distanciaEuclidea(ejemploElegido,vConMascara,vectorDistancias,num_filas,num_atributos);

	//Calculamos las 3 distancias menores
	int indice1, indice2,  indice3;
	calcularMenores(vectorDistancias,indice1, indice2, indice3);

	//Devolvemos la clase en la que más coincidan
	if(v[indice1][num_atributos-1] == v[indice2][num_atributos-1])
		return v[indice1][num_atributos-1];
	if(v[indice2][num_atributos-1] == v[indice3][num_atributos-1])
		return v[indice2][num_atributos-1];
	if(v[indice3][num_atributos-1] == v[indice1][num_atributos-1])
		return v[indice3][num_atributos-1];

	return v[indice1][num_atributos-1];
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo SFS
/////////////////////////////////////////////////////////////////////////////////
void sfs(const vector<vector <double> > &v, vector<vector <double> > &vNormalizado, vector<int> &mascara, int num_filas, int num_atributos, vector<int> &resultado){

	bool fin = false;
	int aciertosAnt = 0;
	vector<int> mascaraAuxiliar = resultado;
	mascara.assign(mascara.size(),0);

	while(!fin){

		//Inicializamos los valores del bucle a 0
		resultado.assign(vNormalizado.size(),0);

		//Calculamos los aciertos de cada clase
		for(int i=0; i < num_atributos-1; ++i){
				mascara[i] = 1;
			for(int j=0; j < num_filas; ++j){
				if(knn(v,vNormalizado,j,mascara,num_filas,num_atributos) == v[j][num_atributos-1])
					resultado[i] += 1;
			}
			if(mascaraAuxiliar[i] != 1)
				mascara[i] = 0;
		}

		//Cogemos la característica más prometedora
		if(resultado[indiceMax(resultado)] > aciertosAnt){
	        mascaraAuxiliar[indiceMax(resultado)] = 1;
	        mascara[indiceMax(resultado)] = 1;
	        aciertosAnt = resultado[indiceMax(resultado)];  
      	}

      	//Si no mejoramos, terminamos
	    else
	    	fin = true;
		
	}

	resultado = mascaraAuxiliar;
}


/////////////////////////////////////////////////////////////////////////////////
// Algoritmo LS
/////////////////////////////////////////////////////////////////////////////////
int funcionObjetivo(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mascara,int num_filas,int num_atributos){

	int coste = 0;

	for(int i=0; i < num_filas; ++i)
		if(knn(v,vNormalizado,i,mascara,num_filas,num_atributos) == v[i][num_atributos-1])
			coste++;
	
	return coste;
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar una máscara inicial aleatoria
/////////////////////////////////////////////////////////////////////////////////
void generarMascaraInicial(int num_atributos, vector<int> & mascara){

	for(int i=0; i < num_atributos-1; ++i){
		if(Rand() >= 0.5)
			mascara.push_back(1);
		else
			mascara.push_back(0);
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar vecinos
/////////////////////////////////////////////////////////////////////////////////
void generarVecinos(vector<int> mascara,vector< vector<int> >&resultado){

	for(int i=0; i < mascara.size(); ++i){
		if(mascara[i] == 0)
			mascara[i] = 1;
		else
			mascara[i] = 0;

		resultado.push_back(mascara);
		
		if(mascara[i] == 0)
			mascara[i] = 1;
		else
			mascara[i] = 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar un vecino de forma aleatoria
/////////////////////////////////////////////////////////////////////////////////
void generarVecinoAleatorio(vector<int> &mascara){

	int pos = Randint(0,mascara.size()-1);
	if(mascara[pos] == 0)
		mascara[pos] = 1;
	else
		mascara[pos] = 0;
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar entornos de tamaño tam
/////////////////////////////////////////////////////////////////////////////////
void generarVecindario(vector<int> &mascara,vector< vector<int> > &resultado,int tam){
  
	resultado.clear();

	for(int i=0; i < tam; ++i){
		auto aux = mascara;
		generarVecinoAleatorio(aux);
		resultado.push_back(aux);
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para encontrar el primer índice en el que los valores de dos máscaras dadas difieren
/////////////////////////////////////////////////////////////////////////////////
int indiceDistinto(vector<int> &mascara1,vector<int> &mascara2){

  	for(int i=0; i < mascara1.size(); ++i)
    	if(mascara1[i] != mascara2[i])
      		return i;
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo LS
/////////////////////////////////////////////////////////////////////////////////
void ls(const vector<vector <double> > &v, vector<vector <double> > &vNormalizado, vector<int> &mascara, int numFilas, int numAtributos, int &resultado, int &n_eval){
	
	int coste;
	bool fin = false;
	int aciertosAnt = resultado;
	int i;

	//Generamos los vecinos
	vector<vector<int> > vecindario;
	generarVecinos(mascara,vecindario);

	for(i=0; i < vecindario.size(); ++i){

		//Calculamos el coste de cada vecino usando KNN
		coste = funcionObjetivo(v,vNormalizado,vecindario[i],numFilas,numAtributos);
		n_eval+=1;

		//Si encontramos una solución de mejor coste, dejamos de evaluar los vecinos y actualizamos la solución
		if(aciertosAnt < coste){
			aciertosAnt = coste;
			mascara = vecindario[i];
			i = vecindario.size()+1;
		}				
	}

	resultado = aciertosAnt;
}

/////////////////////////////////////////////////////////////////////////////////
// Función del operador de mutación 
/////////////////////////////////////////////////////////////////////////////////
void mutacion(vector<int> &mascara, int num_mutaciones){

	for(int i=0; i < num_mutaciones; ++i){

		generarVecinoAleatorio(mascara);
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para encontrar los n mejores elementos
/////////////////////////////////////////////////////////////////////////////////
void nMejores(vector<int> v,vector<int> &resultado, int n){

	int mayor = 0;
	int i_mayor;

	for(int i=0; i < n; ++i){
		for(int j=0; j < v.size(); ++j){
			if(mayor < v[j]){
				mayor = v[j];
				i_mayor = j;
			}
		}
		resultado.push_back(i_mayor);
		v[i_mayor] = -1;
		mayor = 0;
	}
}


/////////////////////////////////////////////////////////////////////////////////
// Función seleccionar los mejores cromosomas
/////////////////////////////////////////////////////////////////////////////////
void torneo(vector<vector<int>>& poblacion,vector<int>costes,vector<vector<int>> & poblacion_seleccionada,vector<int>& costes_seleccionada, int num_torneos){

	for(int i=0; i < num_torneos; ++i){

		int i_primero = Randint(0,poblacion.size()-1);
		int i_segundo = Randint(0,poblacion.size()-1);
		
		if(i_primero != i_segundo){
			if(costes[i_primero] < costes[i_segundo]){
					poblacion_seleccionada.push_back(poblacion[i_segundo]);
					costes_seleccionada.push_back(costes[i_segundo]);
			}
			else{
					poblacion_seleccionada.push_back(poblacion[i_segundo]);
					costes_seleccionada.push_back(costes[i_segundo]);
			}
		}

		else
			--i;
	}
		
}

/////////////////////////////////////////////////////////////////////////////////
// Función para cruzar los cromosomas
/////////////////////////////////////////////////////////////////////////////////
void cruce(const vector<vector<double>> &v, vector<vector <double> > &vNormalizado, vector<vector<int>>& poblacion, vector<int> & coste_poblacion, int total_cruces, int num_filas, int num_atributos, int &n_eval){

	int punto1, punto2;
	int primero = 0;
	int segundo = 1;

	for(int i=0; i < total_cruces; ++i){

		punto1 = Randint(0, poblacion.size()-1);
		punto2 = Randint(punto1, poblacion.size()-1);

		if(punto1 != punto2) {

			for(int j=punto1; j < punto2; ++j){
				int aux = poblacion[primero][j];
				poblacion[primero][j] = poblacion[segundo][j];
				poblacion[segundo][j] = aux;
			}

			coste_poblacion[primero] = funcionObjetivo(v,vNormalizado,poblacion[primero],num_filas,num_atributos);
			coste_poblacion[segundo] = funcionObjetivo(v,vNormalizado,poblacion[segundo],num_filas,num_atributos);
			n_eval+=2;

			primero+=2;
			segundo+=2;
		}

		else
			--i;
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo AM
/////////////////////////////////////////////////////////////////////////////////
void AM(const vector<vector <double> > &v, vector<vector <double> > &vNormalizado, vector<int> &mascara, int num_filas, int num_atributos, vector<int> &resultado, int tipo){

	vector<vector<int>> poblacion,poblacion_sel;
	vector<int> ev_poblacion,ev_poblacion_sel;	
	int num_cruces = 0.7*10/2;
	int num_eval = 0;
	int contador = 0;

	mascara.clear();
	generarMascaraInicial(num_atributos,mascara);

	//Generamos la población
	generarVecindario(mascara,poblacion,10);

	int num_mutaciones = ceil(poblacion.size() * (num_atributos-1) * 0.001); 

	//Evaluamos la población
	for(int i=0; i < 10; ++i)
		ev_poblacion.push_back(funcionObjetivo(v,vNormalizado,poblacion[i],num_filas,num_atributos));
	
	while(num_eval < 1500){

		ev_poblacion_sel.clear();
		poblacion_sel.clear();

		//Seleccionamos
		torneo(poblacion,ev_poblacion,poblacion_sel,ev_poblacion_sel,10);

		//Cruzamos las parejas
		cruce(v,vNormalizado,poblacion_sel,ev_poblacion_sel,num_cruces,num_filas,num_atributos,num_eval);

		//Mutamos
		for(int i=0; i < num_mutaciones; ++i){

			int pos = Randint(0,poblacion_sel.size()-1);

			mutacion(poblacion_sel[pos],1);

			ev_poblacion_sel[pos] = funcionObjetivo(v,vNormalizado,poblacion_sel[pos],num_filas,num_atributos);
		}

		num_eval+=num_mutaciones;

		//Buscamos el mejor de la población anterior y el peor de la nueva población
		int mejor=1;
		int coste_mejor = ev_poblacion[1];
		for(int i=1; i < ev_poblacion.size(); ++i){
			if(coste_mejor < ev_poblacion[i]){
				coste_mejor = ev_poblacion[i];
				mejor = i;
			}
		}

		int peor=1;
		int coste_peor = ev_poblacion_sel[1];
		for(int i=1; i < ev_poblacion_sel.size(); ++i){
			if(coste_peor > ev_poblacion_sel[i]){
				coste_peor = ev_poblacion_sel[i];
				peor = i;
			}
		}

		//Actualizamos la población
		poblacion_sel[peor] = poblacion[mejor];
		poblacion = poblacion_sel;
		ev_poblacion_sel[peor] = ev_poblacion[mejor];
		ev_poblacion = ev_poblacion_sel;
		contador++;

		//Llamamos a la búsqueda local según el tipo AM
		if(contador == 10){
			switch(tipo){

				case 1: 	
					contador = 0;
					for(int i=0; i < 10; ++i)
						ls(v,vNormalizado,poblacion[i],num_filas,num_atributos,ev_poblacion[i],num_eval);
					break;


				case 2:
					contador = 0;
					ls(v,vNormalizado,poblacion[Randint(0,9)],num_filas,num_atributos,ev_poblacion[Randint(0,9)],num_eval);
					break;

				case 3:
					contador = 0;
					ls(v,vNormalizado,poblacion[indiceMax(ev_poblacion)],num_filas,num_atributos,ev_poblacion[indiceMax(ev_poblacion)],num_eval);
					break;
			}
		}
		
	}

	mascara = poblacion[indiceMax(ev_poblacion)];

}

/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/
int main(){

	int opcionMenu;
	int numSimulacion;

	cout << "-----------------" << endl;
	cout << "Menú:" << endl;
	cout << "Opción 1: ejecutar AM(10,1.0)" << endl;
	cout << "Opción 2: ejecutar AM(10,0.1)" << endl;
	cout << "Opción 3: ejecutar AM(10,0.1M)" << endl;
	cout << "Opción 4: salir" << endl;
	cout << "-----------------" << endl;
	//cin >> opcionMenu;
	opcionMenu = 1;
	numSimulacion = 1;
	Seed = 5;

	while(opcionMenu != 4){
	
		cout << "-----------------" << endl;
		cout << "Introduzca el número de simulación(1-30): " << endl;
		cout << "Simulaciones 1-10: movement_libras" << endl;
		cout << "Simulaciones 11-20: arrhythmia" << endl;
		cout << "Simulaciones 21-30: wdbc" << endl;
		//cin >> numSimulacion;
		cout << "Introduzca el valor de la semilla: " << endl;
		//cin >> Seed;



		vector<vector <double> > vectorEntero;
		vector<vector <double> > mitad1;
		vector<vector <double> > mitad2;

		int num_filasI;
		int num_atributos;
		lectorARFF(numSimulacion,vectorEntero,num_filasI,num_atributos);

		//Separamos los datos en dos vectores, mitad en cada uno
		for(int i=0;i<num_filasI;++i){
			if(i<num_filasI/2)
				mitad1.push_back(vectorEntero[i]);
			else
				mitad2.push_back(vectorEntero[i]);
		}

		//Guardamos el numero de filas de cada una de las matrices
		int filasMitad1 =  num_filasI-num_filasI/2-1;
		int filasMitad2 =  num_filasI-filasMitad1;

		vector<double> minimoColumnas;
		vector<double> maximoColumnas;

		clock_t startTime = clock();
		vector<vector <double> > vNormalizado;
		vector<vector <double> > vNormalizado2;
		calcularMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,num_atributos);
		normalizar(mitad1,vNormalizado,filasMitad1,num_atributos,minimoColumnas,maximoColumnas);
		vector<int> resultado;

		switch(opcionMenu){
			case 1:
				cout << "Ejecutando AM(10,1.0). Simulación: " << numSimulacion << endl;
				AM(mitad1,vNormalizado,resultado,filasMitad1,num_atributos,resultado,1);
				break;

			case 2:
				cout << "Ejecutando AM(10,0.1). Simulación: " << numSimulacion << endl;
				AM(mitad1,vNormalizado,resultado,filasMitad1,num_atributos,resultado,2);
				break;

			case 3:
				cout << "Ejecutando AM(10,0.1M). Simulación: " << numSimulacion << endl;
				AM(mitad1,vNormalizado,resultado,filasMitad1,num_atributos,resultado,3);
				break;

		}
 
 		cout << "-----------------" << endl;
		cout << "Resultado: " << endl;
		for(int i=0;i<num_atributos-1;++i)
			cout << resultado[i] << " ";
		cout << endl;

		int aciertos = 0;
		int fallos = 0;
		calcularMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),num_atributos);
		normalizar(mitad2,vNormalizado2,mitad2.size(),num_atributos,minimoColumnas,maximoColumnas);
		for(int i=0;i<mitad2.size();++i)
			if(knn(mitad2,vNormalizado2,i,resultado,mitad2.size(),num_atributos) == mitad2[i][num_atributos-1])
				aciertos++;
			else
				fallos++;

		int seleccionadas = 0;
		for(int i=0; i < resultado.size(); ++i)
			if(resultado[i] == 1)
				seleccionadas++;

		cout << "Número de aciertos: " << aciertos <<  "\tNúmero de fallos: " << fallos << endl;
		cout << "Porcentaje de acierto: " << (float)aciertos*100/(aciertos+fallos) << "%" << "\tPorcentaje de reducción: " << (float)100*(num_atributos-1-seleccionadas)/(num_atributos-1) << "%" << endl;
		cout << "Se han seleccionado " << seleccionadas << " de " << num_atributos-1 << " características." << endl;
		cout << "Tiempo de ejecución: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " segundos." << endl;
		
		cout << "-----------------" << endl;
		cout << "Menú:" << endl;
		cout << "Opción 1: ejecutar AM(10,1.0)" << endl;
		cout << "Opción 2: ejecutar AM(10,0.1)" << endl;
		cout << "Opción 3: ejecutar AM(10,0.1M)" << endl;
		cout << "Opción 4: salir" << endl;
		cout << "-----------------" << endl;
		//cin >> opcionMenu;
	
		if(numSimulacion==30){
			numSimulacion = 1;
			opcionMenu++;
		}

		else
			numSimulacion++;
	}

	cout << "Saliendo..." << endl;
}
