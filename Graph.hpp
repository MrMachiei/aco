/*
#include <iostream>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
using namespace std;
*/
#ifndef GRA_HPP
#define GRA_HPP

class Graph
{
    public: 
    int V;
    int Max;
    int **adj;


    int gen_c_opt;
    int opt_c;
    int naj_zachlanny_c;
  

    Graph(int V, int max);
    ~Graph();
    
    void print();
    void gen_cykl_optimum(); // generuje sciezke opt
    Graph* getGrap();

    int get_gen_c_opt(){return gen_c_opt;} // wynik generatora opt

    void wyczerpujacy_cykl(); // wywoluje void wyczerpujacy_s(int per[], int l, int r, int len);
    void wyczerpujacy_c(int per[], int l, int r, int len);

    int get_opt_c(){return opt_c;}  // wynik wyczerpujacy
    int get_naj_zachlanny_c(){return naj_zachlanny_c;} // wynik najblizszego sasiada


    void zachlanny_cykl(); // najblizszy sasiad

    int** get_graph();
};

#endif //GRA_HPP