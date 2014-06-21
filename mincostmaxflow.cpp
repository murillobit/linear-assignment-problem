#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<string>
#include<vector>
#include<queue>
#include<set>
#include<algorithm>
#include<utility>
#include<stack>
#include<map>

#define MAX_V 2020
#define MAX_N 1010
#define INF 2e9
#define EPS 1e-7
#define min(a,b) a<b ? a : b

using namespace std;

typedef pair<int, double> id;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<id> vii;

bool lessThan(double a, double b){
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

double c[MAX_N][MAX_N];
int res[MAX_V][MAX_V];

class Maxflow {
public:
    int s, t, n, V, mf, f;;
    vi p;
    vector<vii> AdjList;
    
    Maxflow(){
        scanf("%d",&n);
        V = 2*n + 2;
        
        for(int i=0; i<V; i++){
            for(int j=0; i<V; i++){
                res[i][j] = 0;
            }
        }
        
        s = 0;
        t = V - 1;
        AdjList.assign(V, vii());
        
        for(int i=0; i<n; i++){
            res[0][i+1] = 1;
            AdjList[0].push_back(id(i+1, 0.0));
            AdjList[i+1].push_back(id(0, 0.0));
            
            AdjList[n+i+1].push_back(id(n+n+1, 0.0));
            AdjList[n+n+1].push_back(id(n+i+1, 0.0));
            
            for(int j=0; j<n; j++){
                scanf("%lf",&c[i][j]);
                
                res[i+1][n+j+1] = 1;
                AdjList[i+1].push_back(id(n+j+1, c[i][j]));
                AdjList[n+j+1].push_back(id(i+1, -c[i][j]));
                
                res[n+j+1][n+n+1] = 1;
            }
        }
    }
    
    void augment(int v, int minEdge){
        if(v==s){
            f = minEdge;
            return;
        }
        else if(p[v] != -1){
            augment(p[v], min(minEdge, res[p[v]][v]));
            res[p[v]][v] -= f;
            res[v][p[v]] += f;
        }
    }

    void match(){
        //Edmonds-Karp algorithm using SPFA (Shortest Path Faster Algorithm) to find shortest (min cost) paths
        
        mf = 0;
        
        while(1){
            f = 0;
            vd dist(V, INF);
            dist[s] = 0.0;
            p.assign(V, -1);
            
            queue<int> q; q.push(s);
            vi in_queue(V,0); in_queue[s] = 1;
            
            while(!q.empty()){
                int u = q.front(); q.pop(); in_queue[u] = 0;
                
                for (int j = 0; j < (int)AdjList[u].size(); j++) {
                    id v = AdjList[u][j];
                    
                    if(res[u][v.first] > 0){
                        double a = dist[v.first];
                        double b = (dist[u] + v.second);
                        
                        if(lessThan(b,a)){
                            p[v.first] = u;
                            dist[v.first] = b;

                            if(!in_queue[v.first]){
                                q.push(v.first);
                                in_queue[v.first] = 1;
                            }
                        }
                    }
                }
            }
            
            augment(t, INF);
            
            if(f == 0) break;
            
            mf += f;        
        }
    }
    
    void printM(){
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(res[i+1][n+j+1] == 0){
                    printf("%d %d\n",i,j);
                }
            }
        }
    }
    
    void printC(){
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                printf("%.6lf%s", c[i][j], (j<n-1 ? "\t" : ""));
            }
            printf("\n");
        }
    }
    
   
    double matchingValue(){
        double sum = 0.0;
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(res[i+1][n+j+1] == 0){
                    sum += c[i][j];
                }
            }
        }
        return sum;
    }
    
};


int main(){
    Maxflow f;
    f.match();
    f.printM();
    printf("\n%lf\n",f.matchingValue());
    
    return 0;
}
