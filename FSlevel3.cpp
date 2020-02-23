#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<set>
#include<algorithm>
#include<map>
using namespace std;
#define edge pair<int, double>

/*
[1] Hon Nian Chua, Kang Ning, Wing-Kin Sung, Hon Wai Leong, and Lim-soon Wong. 
    Using indirect protein–protein interactions for protein complex prediction.
    Journal of bioinformatics and computational biology, 6(03):435–466, 2008.
*/

/*
FS Weighting for the level 3 neighbors of PPINs using new formula
Author: Almie P. Carajay
Date: February 23, 2020
Project: PGC-IMBUE- Internship

*/

int main(int argc, char **argv) {
	if(argc != 5) {
		printf("Usage: <input filename> <output filename> <mlr filename> <threshold>\n");
		return 0;
	}
	double SFSThresh = atof(argv[4]);

	if(SFSThresh == 0.0) {
		printf("Usage: <input filename> <output filename> <mlr filename> <threshold>\n");
		return 0;
	}
	
	FILE* input = fopen(argv[1], "r");    // argv[1] is the input filename to be read

    if(input == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* output = fopen(argv[2], "w+");
    
    if(output == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* mlrformat = fopen(argv[3], "w+");
    
    if(mlrformat == NULL) {
    	perror("");
    	return 0;
    }
	
    int V, i;                            // V is the number of lines in a file
    int Vnew = 0, Enew = 0;
    double navg = 0.0;
    fscanf(input, "%d %*d %*d", &V);    //&V is the address of V which is the number of nodes

    set<int> G[V];           //set of G[V] which contains level-1 neighbors
    vector<edge > Gw[V];     // edge of Gw[V]
    vector<int> active;     //vector
    map<int,int> mapping;   //mapping

    set <int > :: iterator val; 
    //inserting each node 
    for(int c = 0; c < V; c++) {
        while(true) {
            fscanf(input, "%d %*d", &i);
            if(i == -1){ break;}
            else{ G[c].insert(i-1); }         //inserts the level-1 neighors of Np
        }
        navg += (double)(G[c].size());
        G[c].insert(c);
    }

    navg /= (double)(V);   //navg = navg/V

    set<int> thirdpath[V];

    for(int u = 0; u < V; u++) {
        set<int> firstsecondneighbors;
        set<int> thirdneighbors;

		for(set<int>:: iterator j = G[u].begin(); j != G[u].end(); j++) {    //checks level 1
			int k = *j;

			if(k != u) {
                firstsecondneighbors.insert(k);   //insert level-1 neighbors
                for(set<int>:: iterator l = G[k].begin(); l != G[k].end(); l++) { //checks level 2
                    if(*l >= 0 and *l != u and *l != k){
                        firstsecondneighbors.insert(*l);  //insert level-2 neighbors
                        for(set<int>:: iterator m = G[*l].begin(); m != G[*l].end(); m++){ //checks level 3 neighbors
                            if(*m >=0 and *m != u and *m != k and *m != *l){
                                thirdneighbors.insert(*m);   // insert level 3 neighbors
                            }
                        }
                    }
                }                                    
			}
		}

         /*-- level 2 weighting --*/
        for(set<int>::iterator it = firstsecondneighbors.begin(); it != firstsecondneighbors.end(); it++) {
            int v = *it;
            double interuv = 0.0;
            double diffuv = (double)(G[u].size());
            double diffvu = (double)(G[v].size());
            double lambdauv = 0.0;
            double lambdavu = 0.0;
            
            for(set<int>::iterator it = G[u].begin(); it != G[u].end(); it++) {
                if(G[v].find(*it) != G[v].end()){
                    interuv += 1.0;
                } 

            }
            if(interuv == 0.0) continue;


            diffuv -= interuv;
            diffvu -= interuv;

            lambdauv = max(lambdauv, navg - diffuv - interuv);
            lambdavu = max(lambdavu, navg - diffvu - interuv);

            double FSuv = (4.0*interuv*interuv) / ((diffuv+(2.0*interuv)+lambdauv)*(diffvu+(2.0*interuv)+lambdavu));

            if(FSuv >= SFSThresh) {
                Enew++;
                if(Gw[v].empty()) { Vnew++; active.push_back(v+1); }
                Gw[v].push_back(edge(u, FSuv));
            }

        }

        set <int> thirdOnly;
        set_difference(thirdneighbors.begin(), thirdneighbors.end(), firstsecondneighbors.begin(), firstsecondneighbors.end(), inserter(thirdOnly, thirdOnly.end()));  
        
        for (set <int > :: iterator third = thirdOnly.begin(); third != thirdOnly.end(); ++third){            
            for (set <int > :: iterator th = thirdpath[*third].begin(); th != thirdpath[*third].end(); ++th){ 
                int v = *th;
                int w = *third;
                double interuv = 0.0;
                double interuw = 0.0;
                double intervw = 0.0;
                double interuvw = 0.0;
                double diffuvw = (double)(G[u].size());
                double diffvuw = (double)(G[v].size());     
                double diffwuv = (double)(G[w].size()); 
                double lambdauvw = 0.0;
                double lambdavuw = 0.0;
                double lambdawuv = 0.0;

                /*----- intersection of u and v -----*/
                for(set<int>::iterator xt = G[u].begin(); xt != G[u].end(); xt++) {
                    if(G[v].find(*xt) != G[v].end()){
                        interuv += 1.0;
                    } 
                }

                /*----- intersection of v and w -----*/
                for(set<int>::iterator vt = G[v].begin(); vt != G[v].end(); vt++) {
                    if(G[w].find(*vt) != G[w].end()){
                        intervw += 1.0;
                    } 
                }

                /*----- intersection of u and w -----*/
                for(set<int>::iterator ut = G[u].begin(); ut != G[u].end(); ut++) {
                    if(G[w].find(*ut) != G[w].end()){
                        interuw += 1.0;
                    } 
                }

                /*----- intersection of u, v and w ------*/
                for(set<int>::iterator zt = G[u].begin(); zt != G[u].end(); zt++) {
                    if(G[v].find(*zt) != G[v].end()){
                        if(G[w].find(*zt) != G[w].end()){
                            interuvw += 1.0;
                        }
                    } 
                }

                double unionvw = diffvuw + diffwuv - intervw;
                double unionuw = diffuvw + diffwuv - interuw;
                double unionuv = diffuvw + diffvuw - interuv;  
         
                diffuvw -= unionvw;
                diffvuw -= unionuw;
                diffwuv -= unionuv;

                lambdauvw = max(lambdauvw, navg - (abs(diffuvw) + abs(interuv) + abs(interuw) - abs(interuvw)));
                lambdavuw = max(lambdavuw, navg - (abs(diffvuw) + abs(interuv) + abs(intervw) - abs(interuvw)));
                lambdawuv = max(lambdawuv, navg - (abs(diffwuv) + abs(interuw) + abs(intervw) - abs(interuvw)));

                //calculate the weight for each path the level 3 nodes
                double numerator = 3.0 * interuvw;
                double denum1 = (abs(diffuvw) + (2.0 * abs(interuv)) + (2.0 * abs(interuw)) + (3.0 * abs(interuvw)) + lambdauvw);
                double denum2 = (abs(diffvuw) + (2.0 * abs(interuv)) + (2.0 * abs(intervw)) + (3.0 * abs(interuvw)) + lambdavuw);
                double denum3 = (abs(diffwuv) + (2.0 * abs(interuw)) + (2.0 * abs(intervw)) + (3.0 * abs(interuvw)) + lambdawuv);

                long double FSNew = (numerator/denum1)*(numerator/denum2)*(numerator/denum3);


                if(FSNew >= SFSThresh) {
                    Enew++;
                    if(Gw[w].empty()) { Vnew++; active.push_back(w+1); }

                    Gw[w].push_back(edge(u, FSNew));
                }

            } 
        } 

		if((u+1)%500 == 0)
			printf("Done with node %d.\n", u+1);
    }

 


/*------- printing the weight results ------------*/

    sort(active.begin(), active.end());

    for(int c = 0; c < (int)(active.size()); c++){
    	mapping[active[c]] = c+1;
    }
   
    
    fprintf(output, "%d %d %d\n", V, Vnew, Enew);
    fprintf(mlrformat, "%d %d 1\n", Vnew, Enew);
    
    for(int c = 0; c < V; c++) {
    	if(!Gw[c].empty()) {
    		fprintf(output, "%d", c+1);

    		for(int d = 0; d < (int)(Gw[c].size()); d++) {
    			fprintf(output, " %d %lf", Gw[c][d].first+1, Gw[c][d].second);
    			if(d > 0) fprintf(mlrformat, " ");
    			fprintf(mlrformat, "%d %lf", mapping[Gw[c][d].first+1], Gw[c][d].second*100.0);
    		}
    		fprintf(output, "\n");
    		fprintf(mlrformat, "\n");
    	}
    }
    fclose(input);
    fclose(output);
    fclose(mlrformat);
}

