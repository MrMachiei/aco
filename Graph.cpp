#include <iostream>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include "Graph.hpp"
using namespace std;

Graph::Graph(int v, int max) 
{ 
    this->V = v;
    this->Max = max;
    adj = new int* [V]; 
    for (int i=0; i<V; i++) 
    {
       adj[i] = new int[V]; 
       //memset(adj[i], false, V*sizeof(int)); 
    }
    //srand((unsigned)time(0)); 
    int Rand;
    for (int i=0; i<V; i++)
    {
        for (int j=i; j<V; j++)
        {
            if(i == j)
            {
                adj[i][j] = -1;
            }
            else
            {   
                Rand = (rand() % (max/2))+1+(max/2); 
                adj[i][j] = Rand;
                adj[j][i] = Rand;
            }
            
        }
    }
}
int** Graph::get_graph(){
    return adj;
}

Graph* Graph::getGrap(){
    return this;
}

Graph::~Graph(){
    for(int i = 0; i < V; i++) delete adj[i];
    delete[] adj;
}

void Graph::print() 
{
    ofstream myfile;
    myfile.open ("matrix.txt");
    myfile << V << " " << V << endl;
   for (int u=0; u<V; u++) 
    { 
        for (int v=0; v<V; v++) 
        {
            myfile << adj[u][v];
            cout << adj[u][v];
            if(v<V-1) {
                myfile << " ";
                cout << " ";
                }
        }
      myfile << endl;
      cout << endl;
    } 
    myfile.close();
}

void Graph::gen_cykl_optimum()
{
    gen_c_opt = 0;
    int *opt_c;
    opt_c = new int [V];
    for(int i = 0; i< V-1; i++){
        opt_c[i]=i+1;
    }
    opt_c[V-1] = 0;

    int a,b,temp;
    for(int i = 0; i < V*2; i++){
        a = (rand() % (V-1))+1;
        b = (rand() % (V-1))+1;
        temp = opt_c[a];
        opt_c[a] = opt_c[b];
        opt_c[b] = temp;
    }

    /*
    for(int i = 0; i< V; i++){
        cout <<opt_c[i]<< " -> ";
    }
    cout<<opt_c[0]<<endl;
    */ 

    int Rand;
    for(int i = 0; i< V-1; i++){
                Rand = (rand() % (Max/2))+1; 
                adj[opt_c[i+1]][opt_c[i]] = Rand;
                adj[opt_c[i]][opt_c[i+1]] = Rand;
                gen_c_opt += Rand;  
    }

    if(V>2){
    Rand = (rand() % (Max/2))+1;
    adj[opt_c[0]][opt_c[V-1]] = Rand;
    adj[opt_c[V-1]][opt_c[0]] = Rand;
    }
    gen_c_opt += Rand;
}

void Graph::wyczerpujacy_c(int per[], int l, int r, int len)  {  

    if (per[0]>0) return;
    if (l == r){
        int naj = 0;
        for (int i = 0; i< V-1; i++){
            naj += adj[per[i]][per[i+1]];          
        }
        naj += adj[per[V-1]][0];

        if (naj < opt_c) opt_c = naj;
    }  
    else
    {  
        int temp;
        // Permutations made  
        for (int i = l; i <= r; i++)  
        {  
  
            // Swapping done  
            temp = per[l];
            per[l] = per[i];
            per[i] = temp;
  
            // Recursion  
            wyczerpujacy_c(per, l+1, r, len);  
  
            // backtrack  
            temp = per[l];
            per[l] = per[i];
            per[i] = temp; 
        }  
    }  
} 

void Graph::wyczerpujacy_cykl(){
    int per[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    opt_c = 9999999;
    wyczerpujacy_c(per,0,V-1,V);
}

void Graph::zachlanny_cykl()
{   
    if(V > 1){
        naj_zachlanny_c = 9999999;
        for(int start = 0; start < V; start++){
            int x = 0;
            int wsk, next;
            int suma = 0;
            int najmniejszy;
            wsk = start;


            int *odwiedzona;
            odwiedzona = new int [V];
            odwiedzona[wsk] = 1;
            int *cykl;
            cykl = new int [V];
            cykl[x] = wsk;
            x++;


            for (x = 1; x<V; x++)
            {   
                najmniejszy = 9999999;
                for(int i = 0; i < V; i++)
                {
                    if((najmniejszy > adj[wsk][i]) && (odwiedzona[i] != 1))
                    {   
                        next = i;
                        najmniejszy = adj[wsk][i];
                    }
                }
                wsk = next;
                odwiedzona[wsk] = 1;
                cykl[x] = wsk;
            }

            /*
            cout << endl;
            cout<<cykl[0]+1;
            for(int i = 1; i < V; i++)
            {
                cout<<" -> "<<cykl[i]+1;
            }
            cout << endl;
            */

            for(int i=0; i < V-1; i++)
            {
                suma = suma + adj[cykl[i]][cykl[i+1]];
            }
            suma = suma + adj[cykl[V-1]][cykl[0]];
            if (naj_zachlanny_c>suma)naj_zachlanny_c=suma;
        }
    }
}

/*
int main(){
    Graph g(3,2);
    g.gen_cykl_optimum();
    g.wyczerpujacy_cykl();
    g.zachlanny_cykl();
    g.print();
    cout<< g.get_gen_c_opt()<<" "<<g.get_opt_c() << " "<< g.get_naj_zachlanny_c();
}
*/