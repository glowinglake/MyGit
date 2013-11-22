#include <vector>
#include <string>
#include <iostream>
	
using namespace std;		
class Solution {
public:
    vector<vector<string> > res;
	void NQueens(int n, int N, vector<vector<int> > arr){
		for(int i=0;i<N;i++){
			if(arr[n][i] != 0)
				continue;
			else{
				vector<vector<int> > newArr;
				copyArr(arr,newArr,N);
				colorArr(newArr, n, i, N);
				if(n<N-1)
					NQueens(n+1, N, newArr);
				else
					addConfig(newArr, N);
				}
			}
		}
	void copyArr(vector<vector<int> > arr, vector<vector<int> > newArr, int N){
		for(int i=0;i<N;i++){
			vector<int> row;
			newArr.push_back(row);
				for(int j=0;j<N;j++)
					newArr[i].push_back(arr[i][j]);
			}
		}
	void colorArr(vector<vector<int> > arr, int n, int i, int N){
		arr[n][i]=2;
		for(int x=0;x<N;x++)
			if(arr[n][x]!=2)
				arr[n][x]=1;
		for(int x=0;x<N;x++)
			if(arr[x][i]!=2)
				arr[x][i]=1;
		// diag
		for(int x=n,  y=i;x<N,y<N;x++,y++)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x>=0,y<N;x--,y++)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x>=0,y>=0;x--,y--)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x<N,y>=0;x++,y--)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		}
	void addConfig(vector<vector<int> > arr, int N){
		vector<string> v;
		for(int i=0;i<N;i++){
			string s;
			for(int j=0;j<N;j++){
				//char c=(arr[i][j]==2) ? 'Q':'.';
				s.append(1, (arr[i][j]==2) ? 'Q':'.');
				}
			v.push_back(s);
			cout<<s<<endl;
			}
		cout<<endl;
		res.push_back(v);
	}
    vector<vector<string> > solveNQueens(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<int> > arr;
		for(int i=0;i<n;i++){
			vector<int> row;
			arr.push_back(row);
			for(int j=0;j<n;j++)
				arr[i].push_back(0);
			}
		NQueens(0,n,arr);
		return res;
    }
};
int main(){
	Solution s;
	s.solveNQueens(3);
}