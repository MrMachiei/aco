#include <iostream>
#include "Graph.hpp"
#include <fstream>
#include <math.h>
#include <thread>
#include <chrono>

#define ALPHA 7
#define BETA 1
#define ANTS 10000
#define Q 10
#define POKOLENIA 10
#define RHO 0.65

int wybierzKrawedz(int V,int pozycja, double **feromony, bool* tabu, int** graf){
    double* szanse = new double[V];
    double maks = 0;
    for(int i = 0; i < V; i++){
        if(graf[pozycja][i] > 0 && tabu[i]){
            szanse[i] = pow(feromony[pozycja][i], ALPHA) * pow(1.0/graf[pozycja][i], BETA);
            maks += szanse[i];
        }else{
            szanse[i] = 0;
        }
    }
    
    double suma = 0;
    double* sumaDo = new double[V];
    for(int i = 0; i < V; i++){
        suma += szanse[i];
        sumaDo[i] = suma;
    }
    
    double x = ((double) rand() /RAND_MAX)*maks ;
    int i = 0;
    /*for(int j = 0; j < V; j++){
        std::cout << "SU: " << sumaDo[j] << "  : "  << i << std::endl;
    }*/
    while(sumaDo[i] <= x){
        if( i == V-1 ) break;
        i++;
    }
    delete[] sumaDo;
    delete[] szanse;
    return i;
}

void wyparowywanie(int x, int V, double **feromony, int** graf){
    for(int i = 0; i < V; i++){
        feromony[x][i] *= RHO;
        //feromony[x][i] += Q;
    }
}

void noweFeromony(int V, int** trasy, int** mrowki, double** feromony, int** graf){
    for(int i = 0; i < V; i++){
        wyparowywanie(i, V, feromony, graf);
    }
    for(int i = 0; i < ANTS; i++){
        int droga = mrowki[i][1];
        for(int j = 0; j < V-1; j++) feromony[trasy[i][j]][trasy[i][j+1]] += Q;
        feromony[trasy[i][V-1]][trasy[i][0]] += Q;
    }
}

void zerujTabu(bool** tabu, int V){
    for(int i = 0; i < ANTS; i++)
    for(int j = 0; j < V; j++) tabu[i][j] = 1;
}

bool** utworzTabu(int V){
    bool** tabu = new bool*[ANTS];
    for(int i = 0; i < ANTS; i++){
        tabu[i] = new bool[V];
        for(int j = 0; j < V; j++) tabu[i][j] = 1;
    }
    return tabu;
}

int** utworzMrowki(int V, bool** t, int** trasy){
    int** mrowki = new int*[ANTS];
    for(int i = 0; i < ANTS; i++){
        mrowki[i] = new int[2];
        mrowki[i][0] = i%V;
        mrowki[i][1] = 0;
        t[i][mrowki[i][0]] = 0;
        trasy[i][0] = mrowki[i][0];
    }
    return mrowki;
}

void mrowkiKrok(int V, int krok, int** mrowki, bool** tabu, double** feromony, int** graf, int** trasy){
    for(int i = 0; i < ANTS; i++){
        int nast = wybierzKrawedz(V, mrowki[i][0], feromony, tabu[i], graf);
        tabu[i][nast] = 0;
        mrowki[i][1] += graf[mrowki[i][0]][nast];
        trasy[i][krok] = nast;
        mrowki[i][0] = nast;
    }
}

void mrowkiKrokWsp(int V, int krok, int** mrowki, int start, int stop, bool** tabu, double** feromony, int** graf, int** trasy){
    for(int i = start; i < stop; i++){
        int nast = wybierzKrawedz(V, mrowki[i][0], feromony, tabu[i], graf);
        tabu[i][nast] = 0;
        mrowki[i][1] += graf[mrowki[i][0]][nast];
        trasy[i][krok] = nast;
        mrowki[i][0] = nast;
    }
}

double** nadajFeromon(int V){
    double** feromon = new double*[V];
    for(int i = 0; i < V; i++){
        feromon[i] = new double[V];
        for(int j = 0; j < V; j++){
            if (i == j) feromon[i][j] = 0;
            else feromon[i][j] = 1.0;
        }
    }
    return feromon;
}

int** utworzTrasy(int V){
    int** trasy = new int*[ANTS];
    for(int i = 0; i < ANTS; i++){
        trasy[i] = new int[V];
        for(int j = 0; j < V; j++) trasy[i][j] = 0;
    }
    return trasy;
}

void startujMrowki(int V, bool** t, int** trasy, int** mrowki){
    
    for(int i = 1; i < ANTS; i++){
        int j = rand() % V;
        mrowki[i][0] = j;
        mrowki[i][1] = 0;
        t[i][mrowki[i][0]] = 0;
        trasy[i][0] = mrowki[i][0];
    }
}

int najkrotszaTrasa(int** mrowki){
    int min = mrowki[0][1];
    for(int i = 1; i < ANTS; i++){
        if(mrowki[i][1] < min) min = mrowki[i][1];
    }
    return min;
}

void zerujTrasy(int V, int** trasy){
    for(int i = 0; i < ANTS; i++) for(int j = 0; j < V; j++) trasy[i][j] = 0;
}

int mrowkiIter(Graph* G, int V){
    bool** tabu = utworzTabu(V);
    int** trasy = utworzTrasy(V);
    int** mrowki = utworzMrowki(V, tabu, trasy);
    double** feromony = nadajFeromon(V);
    int wynik;
    
    for(int p = 0; p < POKOLENIA; p++){
        startujMrowki(V, tabu, trasy, mrowki);
        for(int i = 1; i < V; i++){//pokolenie
            mrowkiKrok(V, i, mrowki, tabu,feromony, G->adj, trasy);
        }
        for(int i = 0; i < ANTS; i++){
            mrowki[i][1] += G->adj[trasy[i][V-1]][trasy[i][0]];
        }
        wynik = najkrotszaTrasa(mrowki);
        noweFeromony(V, trasy, mrowki, feromony, G->adj);
        zerujTabu(tabu, V);
        zerujTrasy(V, trasy);
    }
    
    for(int i = 0; i < ANTS; i++) delete[] tabu[i];
    delete[] tabu;
    for(int i = 0; i < ANTS; i++) delete[] mrowki[i];
    delete[] mrowki;
    for(int i = 0; i < V; i++) delete [] feromony[i];
    delete[] feromony;
    for(int i = 0; i < ANTS; i++) delete [] trasy[i];
    delete[] trasy;
    return  wynik;
}

int mrowkiWsp(Graph* G, int V){
    unsigned int ileProc = std::thread::hardware_concurrency();
    int ileMrowekNaProc = ANTS/ileProc;
    std::thread **threads = new std::thread*[ileProc];
    bool** tabu = utworzTabu(V);
    int** trasy = utworzTrasy(V);
    int** mrowki = utworzMrowki(V, tabu, trasy);
    double** feromony = nadajFeromon(V);
    int wynik;
    
    for(int p = 0; p < POKOLENIA; p++){
        startujMrowki(V, tabu, trasy, mrowki);
        for(int i = 1; i < V; i++){//pokolenie
        //std::thread r1(mrowkiKrokWsp, V, i, mrowki, 0, 500, tabu, feromony, G->adj, trasy);
        //std::thread r2(mrowkiKrokWsp, V, i, mrowki, 500, 1000, tabu, feromony, G->adj, trasy);
            for(int j = 0; j < ileProc; j++){
                threads[j] = new std::thread(mrowkiKrokWsp, V, i, mrowki, j*ileMrowekNaProc, (j+1)*ileMrowekNaProc, tabu, feromony, G->adj, trasy);
            }
            for(int j = 0; j < ileProc; j++){
                threads[j]->join();
            }
            for(int i = 0; i < ileProc; i++) delete threads[i];
        }
        for(int i = 0; i < ANTS; i++){
            mrowki[i][1] += G->adj[trasy[i][V-1]][trasy[i][0]];
        }
        wynik = najkrotszaTrasa(mrowki);
        noweFeromony(V, trasy, mrowki, feromony, G->adj);
        zerujTabu(tabu, V);
        zerujTrasy(V, trasy);
    }
    
    //for(int i = 0; i < ileProc; i++) delete threads[i];
    delete[] threads;
    for(int i = 0; i < ANTS; i++) delete[] tabu[i];
    delete[] tabu;
    for(int i = 0; i < ANTS; i++) delete[] mrowki[i];
    delete[] mrowki;
    for(int i = 0; i < V; i++) delete [] feromony[i];
    delete[] feromony;
    for(int i = 0; i < ANTS; i++) delete [] trasy[i];
    delete[] trasy;

    return  wynik;
}


int main(){
    srand(time(NULL));
    for(int V = 2; V <= 13; V++){
    Graph *test = new Graph(V, 100);
    std::cout << "V = " << V << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Value: " << mrowkiWsp(test, V) << " ";
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); 
    std::cout <<"Synchronus time: "<<  duration.count()/1000.0<< std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::cout << "Value: " << mrowkiIter(test, V) << " ";
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); 
    std::cout <<"Iteration time: "<<  duration.count()/1000.0<< std::endl;

    start = std::chrono::high_resolution_clock::now();
    test->wyczerpujacy_cykl();
    std::cout<< "Value: " <<test->get_opt_c()<< " " ;
    stop = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout <<"Brute-force time: "<<  duration2.count()/1000000.0<< std::endl;

    start = std::chrono::high_resolution_clock::now();
    test->zachlanny_cykl();
    std::cout<< "Value: " << test->get_naj_zachlanny_c()<<  " " ;
    stop = std::chrono::high_resolution_clock::now();
    duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout <<"Nearest neighbour time: "<<  duration2.count()/1000000.0<< std::endl;
    std::cout << std::endl << std::endl;
    delete test;
    }
    return 0;
}