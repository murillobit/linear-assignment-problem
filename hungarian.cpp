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

#define MAX_V 1010
#define MAX_U 1010
#define MAX_UV 2020
#define INF 2e9
#define EPS 1e-7
#define min(a,b) a<b ? a : b

using namespace std;

typedef pair< double, pair<int,int> > edge;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;
double c[MAX_V][MAX_U];

bool lessThan(double a, double b){
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}


class Hungarian{
public:
    double a[MAX_V], b[MAX_V], x[MAX_V], minC;
    int M[MAX_UV], parent[MAX_UV], visited[MAX_UV];
    vector<int> eq[MAX_UV];
    set<int> star;
    int n, ci, cj;

    Hungarian(){
        scanf("%d",&n);
        
        b[0] = INF;
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                scanf("%lf",&c[i][j]);

                if(lessThan(c[i][j],b[0])){
                    b[0] = c[i][j];
                }
            }
            
            a[i] = 0;
            eq[i].clear();
            eq[i+n].clear();
            M[i] = M[i+n] = -1;
        }
        
        for(int j=0; j<n; j++) b[j] = b[0];
    }
    
    void classic(void){
        for(int k=1; k<=n; k++){

            bool augment = false;
            int last = -1;
            while(!augment){
                queue<int> trees;
                for(int i=0; i<n; i++) eq[i].clear(), eq[i+n].clear();
                
                for(int i=0; i<n; i++){
                    if(M[i] == -1){
                        //root nodes of hungarian trees
                        trees.push(i);
                    }
                    for(int j=0; j<n; j++){
                        if(fabs(c[i][j] - (a[i]+b[j])) < EPS){
                            //add admissible edges to equality subgraph
                            eq[i].push_back(n+j);
                            eq[n+j].push_back(i);
                        }
                    }
                }

                star.clear();

                while(!trees.empty()){
                    int s = trees.front(); trees.pop();
                    for(int i=0; i<n; i++) { parent[i] = parent[i+n] = -1; visited[i] = visited[i+n] = 0; }
                    queue<int> q;
                    q.push(s);
                    
                    while(!q.empty()){
                        //bfs to find an augmenting path
                        int node = q.front(); q.pop();
                        star.insert(node);
                        
                        if(node >= n && M[node]==-1){
                            augment = true;
                            last = node;
                        }
                        
                        if(!visited[node]){
                            for(int i=0; i<(int)eq[node].size(); i++){
                                int son = eq[node][i];
                                
                                //search through alternating paths
                                if(!visited[son] && parent[node] != son && ((node<n && M[son]!=node)||(node>=n && M[son]==node))){
                                    q.push(son);
                                    parent[son] = node; //record the path
                                }
                            }
                            visited[node]=1;
                        }
                    }
                    if(augment) break;
                }
                
                if(!augment){
                    double lambda=INF;
                    
                    for(int i=0; i<n; i++){
                        for(int j=0; j<n; j++){
                        
                            if(star.count(i) && !star.count(j+n)){
                                double diff = c[i][j] - a[i] - b[j];
                                
                                if(lessThan(diff,lambda)){
                                    lambda = diff;
                                }
                            }
                        }
                    }

                    for(int i=0; i<n; i++){
                        //change dual variables
                        a[i] += (star.count(i))? lambda : 0;
                        b[i] += (star.count(i+n))? -(lambda) : 0;
                    }
                }
                
            }

            int node = last;
            
            while(parent[node] >= 0){
                //augment current matching
                
                if(M[node]==-1){
                    M[node] = parent[node];
                    M[parent[node]] = node;
                }
                else{
                    M[parent[node]] = -1;
                }
                node = parent[node];
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
    
    void printAB(){
        for(int i=0; i<n; i++){
            printf("%lf%s", a[i], (i<n-1 ? "\t" : ""));
        }
        printf("\n");
        
        for(int i=0; i<n; i++){
            printf("%lf%s", b[i], (i<n-1 ? "\t" : ""));
        }
        printf("\n");
    }
    
    void printM(){
        for(int i=0; i<n; i++){
            printf("%d %d\n", i, M[i]-n);
        }
    }
    
    double matchingValue(){
        double sum = 0.0;
        
        for(int i=0; i<n; i++)
            if(M[i] != -1)
                sum += c[i][M[i]-n];
                
        return sum;
    }
    
};

int main(){
    Hungarian h;
    h.classic();
    h.printM();
    printf("\n%lf\n",h.matchingValue());
    
    return 0;
}
