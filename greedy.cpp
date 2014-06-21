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

#define MAX_V 1000
#define MAX_U 1000
#define MAX_UV 2000
#define INF 2e9
#define min(a,b) a<b ? a : b

using namespace std;

typedef pair< double, pair<int,int> > edge;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;

double c[MAX_V][MAX_U];

class Greedy{
public:
    
    int M[MAX_UV], visited[MAX_UV];
    int n;
    vector< edge > q;
    
    Greedy(){
        scanf("%d",&n);
        
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                scanf("%lf",&c[i][j]);
                q.push_back(edge(c[i][j], make_pair(i,j)));
            }
            M[i] = M[i+n] = -1;
        }
    }
    
    void printM(){
        for(int i=0; i<n; i++){
            printf("%d %d\n", i, M[i]-n);
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
        
        for(int i=0; i<n; i++)
            if(M[i] != -1)
                sum += c[i][M[i]-n];
                
        return sum;
    }
    
    void match(void){

        sort(q.begin(), q.end());
        
        vector< edge >::iterator it;
        int m_size;
        
        for(it = q.begin(), m_size=0; it != q.end() && m_size<n; it++){
            
            int v = it->second.first;
            int u = it->second.second + n;
            
            if(M[v] == -1 && M[u] == -1){
                M[v] = u;
                M[u] = v;
                
                m_size++;
            }
        }
    }
};

int main(){
    Greedy g;
    g.match();
    g.printM();
    printf("\n%lf\n",g.matchingValue());
    
    return 0;
}
