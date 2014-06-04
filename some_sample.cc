// scratch of some coding problem solutions(C++) in bundle
// no separate files for each problem, just scratch
// search for anything interesting

class Solution {
public:
    vector<int> inorderTraversal(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  vector<int> *res=new vector<int>;
	  stack<TreeNode*> stk;
	  TreeNode *cur=root;
	  TreeNode *prev=NULL:
	  stk.push(cur);
	  if(!root)
		return *res;
	  while(!stk.empty()) {
		cur=stk.top();
		if(cur->left) {
		  stk.push(cur->left);
		} else {
		  
		  stk.pop();
		  res->push_back(cur->val);
		  if(cur->right)
			stk.push(cur->right);
		}
	  }
	  return *res;
    }
};
#include <iostream>
using namespace std;
//longest non-decreasing subsequence
int find(int val, int arr[], int n, int s, int e) {
  cout<<"s="<<s<<"      e="<<e<<"val="<<val<<endl;
  if(val >= arr[e])
    return e+1;
  else if(val<arr[s])
    return s;
  while(s<e) {
    int mid=(s+e)/2;  
    if(val >= arr[mid])
      s=mid+1;
    else
      e=mid;
  }
  cout<<"idx="<<s<<endl;
  return s;
}
//longest increasing subsequence
int find2(int val, int arr[], int n, int s, int e) {
  //cout<<"s="<<s<<"      e="<<e<<"val="<<val<<endl;
  if(val > arr[e])
    return e+1;
  else if(val == arr[e])
    return -1;
  else if(val<arr[s])
    return s;
  while(s<e) {
    int mid=(s+e)/2;  
    if(val > arr[mid])
      s=mid+1;
    else if(val==arr[mid])
      return -1;
    else
      e=mid;
  }
  if(arr[s]==val)
    return -1;
  else
    return s;
}
int LIS(int arr[], int n) {
  if(n<=1)
    return n;
  int endNumForEachLen[n+1];
  int curMaxLen=1;
  endNumForEachLen[1]=arr[0];
  for(int i=1; i<n; ++i) {
    int val=arr[i];
    int idx=find2(val, endNumForEachLen, n+1, 1, curMaxLen);
    //always update an element
    if(idx != -1)
      endNumForEachLen[idx]=val;
    if(idx==curMaxLen+1)
      curMaxLen++;
  }
  return curMaxLen;
}
int main() {
  int a[]={1,5,3,2,4,0,6,6,9,9,8,9,1,1,1,1,1,1,1,1,1,1,1};
  cout<<LIS(a, 23)<<endl;
}


void LIS(int arr[], int len) {
  map<int, pair<int, int> > endVal2LenIndex;
  for(int i=0; i<len; i++) {
    auto it = endVal2LenIndex.uppper_bound(arr[i]);
    if(it==endVal2LenIndex.end())
      endVal2LenIndex[arr[i]] = 
#include <map>
#include <list>
#include <iostream>
using namespace std;
template<class T>
class LRU{
  struct node{
    int key;
    int content;
  };
  list<node> List;
  map<int, typename list<node>::iterator> key2Node;
  int maxSize;
public:
  LRU(int mSize):maxSize(mSize) {};
  bool find(int key) {
    return (key2Node.find(key) != key2Node.end());
  }
  T& get(int key) {
    if(key2Node.find(key) != key2Node.end())
      return key2Node.find(key)->content;
  }
  void insert(int key, T content) {
    if(!find(key)) {
      node one;
      one.key=key; one.content=content;
      List.push_front(one);
      key2Node[key] = List.begin();
      if(List.size() > maxSize) {
        
        node back=List.back();
        List.pop_back();
        key2Node.erase(back.key);
        cout<<"remove "<<back.content<<endl;
      }
    } else {//update value only
      List.erase(key2Node[key]);
      key2Node.erase(key);
      node one;
      one.key=key; one.content=content;
      List.push_front(one);
      key2Node[key]=List.begin();
    }
  }
  void print() {
    for(typename list<node>::iterator it=List.begin(); it != List.end(); it++)
      cout<<it->content<<endl;
  }
};

int main(){
  LRU<double> lru(3);
  cout<<lru.find(2)<<endl;
  lru.insert(1,1.1);
  lru.insert(2,2.2);
  lru.insert(3,3.3);
  lru.insert(4,4.4);
  lru.insert(2,2.2);
  lru.print();
  return 0;
}
#include <stdio.h>

int adder(int a, int b) {
  if(a==0)
    return b;
  if(b==0)
    return a;
  int carry=(a&b) << 1;
  int sum=a|b;
  return adder(carry, sum);
}

int multiply(int a, int b) {
  int k=1, sum=0;
  int bb=b;
  while(k) {
    if(a & k) {
      sum=adder(sum, bb);
    }
    k=k<<1;
    bb=bb<<1;
  }
  return sum;
}

int main(){
  int res=
  multiply(1243564543, 35346);
  printf("%d\n", res);
  return 0;
}
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
using namespace std;
class Solution {
public:
    vector<string> anagrams(vector<string> &strs) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<pair<string, string> > sorted2orig;
        set<string> seen;
        set<string> dups;
        vector<string> res;
        for(int i =0; i != strs.size()-1; i++) {
            string s=strs[i];
            string sorted=s;
            sort(sorted.begin(), sorted.end());
            if(seen.find(sorted) != seen.end())
                dups.insert(sorted);
            seen.insert(sorted);
            pair<string, string> p;
            p.first=sorted; p.second=s;
            sorted2orig.push_back(p);
        }
        cout<<sorted2orig.size()<<endl;
        for(int pi=0; pi<sorted2orig.size(); pi++) {
            pair<string, string> p = sorted2orig[pi];
            if(dups.find(p.first) != dups.end())
                res.push_back(p.second);
        }
        return res;
    }
};

int main(){
  Solution s;
  string e;
  vector<string> v;
  v.push_back(e);
  v.push_back(e);
  s.anagrams(v);
  return 0;
}


class AVLTree{
  struct node {
    int height;
    K key;
    T val;
    node *left;
    node *right;
  };
public:
  node* find(K key) {
    node *cur = root;
    while(cur) {
      if(cur->key == key)
        return cur;
      else if(cur->key > key)
        cur=cur->left;
      else
        cur=cur->right;
    }
    return NULL;
  }
  void insert(K key, T val) {
    if(node* n = find(key)) {
      n->val = val;
      return;
    }
    node *n = new node(key, val);
    node *cur=root;
    while(cur) {
      if(cur->key > key)
        cur=cur->left;
      else
        cur=cur->right;
    }


  };
  void remove(K key);

private:
  node* root;
};


#include <iostream>
class a{
  int t;
public:
  void haha() {int i=1; i++;std::cout<<"ahah";
    t++;}
};

int main(){
  a* ii=0;
  ii->haha();
  return 1;
}
class AVL{
  struct Node{
	int val;
	Node *left;
	Node *right;
	Node *par;
  };
public:
  bool insert(int n);
  bool remove(int n);
  Node* find(int n);
  int height(Node *n);
  Node *root;
  int size;
};

Node* AVL::find(int x){
  Node *n=root;
  while(n){
	if(n->val == x)
	  return n;
	else if(n->val > x)
	  n=n->left;
	else
	  n=n->right;
  }
  return n;
}

bool AVL::insert(int x){
  Node *n = new Node;
  n->val=x;
  size++;
  if(!root){
	root=n; return;
  }
  // root not empty
  Node *cur=root;
  Node *prev=NULL;
  while(cur){
	if(cur->val > x){
	  prev=cur; cur=cur->left;
	} else {
	  prev=cur; cur=cur->right;
	}
  }
  if(prev->val > x)//insert on left NULL
	prev->left=n;
  else
	prev->right=n;
  
  while(prev){
	int heightL=prev->left; int heightR=prev->right;
	int heightBalance=heightL-heightR;
	if(heightBalance>1 || heightBalance<-1) {
	  Node *big=heightL>heightR ? prev->left : prev->right; 
	  Node *small=  ;
	  int bigL=height(big->left), bigR=height(big->right);
	  Node *bigBig=bigL>bigR ? big->left : big->right;
	  Node *bigSmall=   ;
	  // determine the case and rotate tree, set prev to new 'prev'
	  
	}
	prev=prev->parent;
  }
}
class Solution {
public:
  TreeNode *build(vector<int> &pre, int a, b, vector<int> &in, int, x, y) {
    TreeNode* n=new TreeNode(pre[0]);
    int ini=index[n->val];
    int leftLen=ini-x;
    int rightLen=y-ini;
    if(leftLen>0)
      n->left=build(pre, a+1, a+leftLen, in, x, x+leftLen-1);
    if(rightLen>0)
      n->right=build(pre, b-rightLen+1, b, in, y-rightLen+1, y);
    return n;
  }
  map<int, int> index;
  TreeNode *buildTree(vector<int> &preorder, vector<int> &inorder) {
    // IMPORTANT: Please reset any member data you declared, as
    // the same Solution instance will be reused for each test case.
    if(preorder.size()==0)
      return NULL;
    index.clear();
    for(int i=0; i<inorder.size(); i++)
      index[inorder[i]]=i;
    return build(preorder, inorder);
  }
};
class Solution {
public:
  void recAdd(vector<vector<int>> &res, int sum, vector<int> &cur, int ci, int target, vector<int> &cand) {
	if(sum==target) {
	  res.push_back(cur);
	  return;
	}
	if(sum>target)
	  return;
	if(i<0)
	  return;
	cur.push_back(cand[i]);
	recAdd(res, sum+cand[i], cur, ci, target, cand);
	cur.pop_back();
	recAdd(res, sum, cur, ci-1, target, cand);
  }
    vector<vector<int> > combinationSum(vector<int> &candidates, int target) {
        // IMPORTANT: Please reset any member data you declared, as
        // the same Solution instance will be reused for each test case.
	  
		
        sort(candidates.begin(), candidates.end());
        vector<vector<int>> res;
		if(candidates.size()==0)
		  return res;
		vector<int> cur;
		recAdd(res, 0, cur, candidates.size()-1, target, candidates);
		return res;
    }
};
#include <map>
#include <iostream>
#include <string>
using namespace std;
void computeDecimal(int a, int b) {
  string res;
  while(a >= b) {
    int c=a/b;
    a = a-c*b;
    char cc=c+'0';
    res=res+cc;
  }
  res+='.';
  std::map<int, int> m;
  while(a) {
    a=a*10;
    if(m.find(a) != m.end()) {
      int idx=m[a];
      res.insert(idx, 1, '(');
      res+=')';
      break;
    }
    m[a]=res.length();
    int c=a/b;
    char cc='0'+c;
    res=res+cc;
    a=a-c*b;
  }
  std::cout<<res<<std::endl;
}

int main() {
  computeDecimal(213,2144);
  computeDecimal(1, 9);
  computeDecimal(0, 7);
  computeDecimal(7, 7);
  computeDecimal(6, 7);
  computeDecimal(2, 7);
  computeDecimal(22, 7);
  computeDecimal(14, 7);
  return 0;
}


class Solution {
public:
  int min(int a, int b){
	return (a<b)?a:b;
  }
  int maxArea(vector<int> &height) {
	  int len=height.size();
	  if(len<2)
		return 0;
	  int max;
	  int left=0;
	  int right=len-1;
	  int prevLH=-1, prevRH=len;
	  while(left<right){
		if(prevLH>-1 && height[left]<=height[prevLH]){
		  left++;
		continue;
		}
		if(prevRH<len && height[right]<=height[prevRH]){
		  right--;
		  continue;
		}
		// compare
		curMax=min(height[left], height[right]) * (right-left);
		if(curMax>max)
		  max=curMax;
		if(height[left]>height[right]){
		  prevRH=right;
		  right--;
		}
		else{
		  prevLH=left;
		  left++;
		}
	  }
	  return max;
    }
};
ListNode* copyListWithRandomPtr(ListNode *head) {
  ListNode *cur=head;
  if(!head)
	return NULL;
  while(cur) {
	ListNode *nextN=cur->next;
	ListNode* inst=new ListNode;
	inst->next=nextN;
	cur->next=inst;
	cur=nextN;
  }
  cur=head;
  ListNode *copy=head->next;
  ListNode *res=copy;
  while(cur) {
	copy->random = cur->random->next;
	cur=copy->next;
	if(cur)
	  copy=cur->next;
  }
  cur=head; copy=head->next;
  while(cur) {
	if(!copy->next) {
	  cur->next=NULL;
	  break;
	} else {
	  cur->next=copy->next;
	  copy->next=cur->next->next;
	  cur=cur->next;
	  copy=copy->next;
	}
  }
  return res;
}
//??getLeftChildNode(TreeNode)???????
//??getRightChildNode(TreeNode)???????
//??isNullNode(TreeNode)????????
int count_complete_binary_tree_nodes(TreeNode root) {
  if(isNullNode(root))
	return 0;
  int leftD=0, rightD=0;
  TreeNode n=root;
  while(1) {
	n=getLeftChildNode(n);
	if(isNullNode(n))
	  break;
	else {
	  leftD++;
	}
  }
  n=root;
  while(1) {
	n=getRightChildNode(n);
	if(isNullNode(n))
	  break;
	else {
	  rightD++;
	}
  }
  if(leftD==rightD)
	return ((1<<(n+1)) -1);
  else
	return count_complete_binary_tree_nodes(getLeftChildNode(root)) + count_complete_binary_tree_nodes(getRightChildNode(root));
}


// toposort
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

void dfs(char c, vector<char>& res, set<char>& visited, map<char, set<char> >& m) {
  if(visited.find(c) != visited.end())
    return;
  for(set<char>::iterator it = m[c].begin(); it != m[c].end(); it++)
    dfs(*it, res, visited, m);
  res.push_back(c);
  visited.insert(c);
}

void dictOrder(vector<string> &v) {
  map<char, set<char> > m;
  for(int i=0; i < v.size(); i++)
    for(int j=0; j < v[i].length(); j++)
      for(int k=j+1; k < v[i].length(); k++) {
        if(m.find(v[i][j])==m.end()) {
          //set<char> s;
          m[v[i][j]] = set<char>();
          //m.insert(pair<char, set<char> >(v[i][j], set<char>())); 
          //m[v[i][j]] = one;
        }
        if(v[i][j] != v[i][k])
          m[v[i][j]].insert(v[i][k]);
      }
  vector<char> res;
  set<char> visited;
  for(map<char, set<char> >::iterator it=m.begin(); it != m.end(); it++)
    dfs(it->first, res, visited, m);
  for(int s=0; s<res.size(); s++)
    cout<<res[s];
}

int main() {
  vector<string> v;
  v.push_back("wrt");
v.push_back("wrf");
v.push_back("er");
v.push_back("ett");
 v.push_back("rftt");
 dictOrder(v);
 return 0;
}
#include <iostream>
using namespace std;



class Solution {
public:
    long divideLong(long dividend, long divisor) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  //cout<<"dividing"<<dividend<<"and"<<divisor<<endl;
      if(divisor<0)
		return -divideLong(dividend, -divisor);
        if(dividend<0)
            return -divideLong(-dividend, divisor);
	  if(divisor==1)
		return dividend;
	  if(dividend<divisor)
		return 0;
	  if(dividend==divisor)
		return 1;
	  int res=1;
	  int temp=divisor;
	  while(temp<=dividend) {
          if(temp<0)
            break;
		temp=temp<<1;
		res=res<<1;
	  }
	  temp=temp>>1;
	  res=res>>1;
	  int newDividend=dividend-temp;
	  return res+divideLong(newDividend, divisor);
	}
    int divide(int dividend, int divisor){
	  if(dividend==0)
		return 0;
	  bool neg = false;
	  if((dividend>0) != (divisor>0))
		neg = true;
	  long divid=(dividend>0 ? (long)dividend:-((long)dividend));
	  long divis=(divisor>0 ? (long)divisor:-((long)divisor));
	  int res = divideLong(divid, divis);
	  if(neg)
		return -res;
	  else
		return res;
    }
};

int main() {
  Solution s;
  //cout<<s.divide(-2147483648,1);
}
class Solution {
public:
    int minDistance(string word1, string word2) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  int len1=word1.length();
	  int len2=word2.length();
	  int dist[len1+1][len2+1];
	  for(int i=0; i<=len1;i++)
		dist[i][0]=i;
	  for(int i=0; i<=len2;i++)
		dist[0][i]=i;
	  for(int i=1;i<=len1;i++)
		for(int j=1;j<=len2;j++) {
		  if(word1[i-1] == word2[j-1]) {
			dist[i][j]=min(dist[i-1][j-1], min(dist[i][j-1], dist[i-1][j])+1);
		  } else {
			dist[i][j]=min(dist[i-1][j-1]+1, min(dist[i][j-1], dist[i-1][j])+1);
		  }
		}
	  return dist[len1][len2];
    }
};
class Solution {
public:
    int numDecodings(string s) {
	  if(s.length()>1 && s[0] == '0')
		return 0;
      if(s.length() == 0)
		return 0;
	  if(s.length() == 1){
		if(s[0]=='0')
		  return 0;
		else
		  return 1;
	  }
	  bool isDigit=false;

	  int digit=(s[0]-'0')*10+s[1]-'0';
	  if(digit>=1 && digit<=26)
		isDigit=true;
	  if(s.length() == 2){
		if(!isDigit){
		  string s4=s;
		  s4.erase(0,1);
		  return numDecodings(s4);
		}
		else {
		  string s5=s;
		  s4.erase(0,1);
		  return 1+numDecodings(s4);
		}
	  }
	  else{
		if(isDigit){
		  string s1=s, s2=s;
		  s1.erase(0,1);
		  s2.erase(0,2);
		  return numDecodings(s1)+numDecodings(s2);
		}else{
		  string s3=s;
		  s3.erase(0,1);
		  return numDecodings(s3);
		}
	  }
	}
};
void escapeSpace(char *str) {
  int space=0;
  int len=strlen(str);
  for(int i=0; i<len; i++)
	if(str[i]==' ')
	  space++;
  int back=i+space*2;
  str[back]=0;
  back--;
  i--;
  while(back>=0 && i>=0) {
	if(str[i]==' '){
	  str[back--]='0';
	  str[back--]='2';
	  str[back--]='%';
	  i--;
	} else
	  str[back--]=str[i--];
  }
  return;
}
int eval(const string & expr, int a, int b) {
  int addSub=-1;
  int mult=-1;
  int num=0;
  bool sub=false;
  for(int i=a;i<=b;i++) {
	if(expr[i]=='+' || expr[i]=='-'){
	  addSub=i;
	  if(expr[i]=='-')
		sub=true;
	} else if(expr[i]=='*') {
	  mult=i;
	} else
	  num=num*10+expr[i]-'0';
  }
  if(addSub != -1 && !sub)
	return eval(expr,a,addSub-1)+eval(expr,addSub+1,b);
  if(addSub != -1 && sub)
	return eval(expr,a,addSub-1)-eval(expr,addSub+1,b);
  if(mult != -1)
	return eval(expr,a,mult-1)*eval(expr,mult+1,b);
  return num;
}
int evaluate(const string& expr) {
  int test[10];
  test[11]=1;
  return eval(expr,0,expr.length()-1);
}
int excelToDec(string excelNum) {
  return e2d(excelNum)+1;
}
//????????excel?
string d2e(int i) {
  int a=i/26;
  int b=i%26;
  string tail;
    if(b==0){
      tail.append(1,'Z'); a--;
    }
    else
      tail.append(1, b+'A'-1);

  string res;
  if(a>0)
	res=d2e(a);
  res.append(tail);
  return res;
}

string decToExcel(int decNum) {
  return d2e(decNum);
  
}
//?excel????????
int e2d(string s) {
  int len=s.length();
  int res=0;
  res+=(s[len-1]-'A')+1;
  if(len>1)
	res+=26*e2d(s.substr(0, len-1));
  return res;
}



int excelToDec(string excelNum) {
  return e2d(excelNum);
}
// using standard exceptions
#include <iostream>
#include <exception>
using namespace std;

class myexception: public exception
{
  virtual const char* what() const throw()
  {
    return "My exception happened";
  }
};

int main () {
  try
  {
    throw  myexception();
  }
  catch (exception& e)
  {
    cout << e.what() << '\n';
  }
  return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
/* Spawn a child process running a new program. PROGRAM is the name
of the program to run; the path will be searched for this program.
ARG_LIST is a NULL-terminated list of character strings to be
passed as the program"s argument list. Returns the process ID of
the spawned process. */
int spawn (char* program, char** arg_list)
{
pid_t child_pid;
/* Duplicate this process. */
child_pid = fork ();
if (child_pid != 0)
/* This is the parent process. */
return child_pid;
else {
/* Now execute PROGRAM, searching for it in the path. */
execvp (program, arg_list);
/* The execvp function returns only if an error occurs. */
fprintf (stderr, "an error occurred in execvp\n");
abort ();
}
}
int main ()
{
/* The argument list to pass to the "ls" command. */
char* arg_list[] = {
"ls", /* argv[0], the name of the program. */
"-l",
"/",
NULL /* The argument list must end with a NULL. */
};
/* Spawn a child process running the "ls" command. Ignore the
returned child process ID. */
spawn ("ls", arg_list);
printf ("done with main program\n");
return 0;
}
#include <iostream>
#include <algorithm>
int find(int arr[], int a, int b, int k) {
  if(a>b || k<0)
    throw;
  int pivot=arr[a], fill=a;
  std::swap(arr[a], arr[b]);
  for(int i=a;i<b;i++) {
    if(arr[i] <= pivot) {
      std::swap(arr[fill++], arr[i]);
    }
  }
  //a ... fill-1 fill ... b
  std::swap(arr[fill], arr[b]);
  if(k == fill-a+1)
    return pivot;
  else if(k > fill-a+1)
    return find(arr, fill+1, b, k-(fill-a+1));
  else
    return find(arr, a, fill-1, k);
}

int main() {
  int a[]={32,73,6l,32,7,3,7};
  std::cout << find(a, 0, 5, 1);
}
bool find(int i, int j, vector<vector<char> > &visited, vector<vector<char> > &grid, string pattern, int idx){
  if(i<0 || i>visited.size() || j<0 || j>visited[0].size())
	return false;
  if(visited[i][j] == 'v')
	return false;
  if(idx == pattern.length()-1 && grid[i][j]==pattern[idx])
	return true;
  if(grid[i][j] != pattern[idx])
	return false;
  // match but not last char
  visited[i][j]='v';
  bool found = find(i+1, j, visited, grid, pattern, idx+1) || find(i-1, j, visited, grid, pattern, idx+1) || find(i, j-1, visited, grid, pattern, idx+1) || find(i, j+1, visited, grid, pattern, idx+1);
  visited[i][j]='n';
  return found;
}

bool exists(vector<vector<char> > &grid, string pattern) {
  vector<vector<char> > v=grid;
  for(int i=0; i<grid.size(); i++)
	for(int j=0; j<grid[0].size(); j++){
	  bool res=find(i,j, v, grid, pattern, 0);
	  if(res)
		return true;
	}
  return false;
}
#include <stdlib.h>
#include <iostream>
int flipCoin() {
  int res=rand()%2;
  return res;
}

int flipUntilHead() {
  int cnt=0;
  while(flipCoin() == false) {
    cnt++;
  }
  return cnt;
}


int doFlipNTimes(int N) {
  if(N<=0)
    return 0;
  int totalTail=0;
  while(N--) {
    totalTail += flipUntilHead();
    std::cout<<totalTail<<std::endl;
  }
  return totalTail;
}

class simulateFlipNTimes {
private:
  simulateFlipNTimes() {};
public:
  static simulateFlipNTimes* sim;
  static void simulate(int N) {
    if(!sim) {
      sim = new simulateFlipNTimes();
    }
    std::cout<<doFlipNTimes(N)<<std::endl;
  }
};

simulateFlipNTimes* simulateFlipNTimes::sim=NULL;

int main() {
  simulateFlipNTimes::simulate(433);
};
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
int main ()
{
pid_t child_pid;
printf ("the main program process ID is %d\n", (int) getpid ());
child_pid = fork ();
if (child_pid != 0) {
printf ("this is the parent process, with id %d\n", (int) getpid ());
printf ("the childs process ID is %d\n", (int) child_pid);
}
else
printf ("this is the child process, with id %d\n", (int) getpid ());
return 0;
}
#include <iostream>
#include <stack>
using namespace std;

void solve(int s, int m, int d, int num) {
  if(num==1) {
    cout<<"move disk "<<num<<" from tower "<<s<<" to tower "<<d<<endl;
  } else {
    solve(s, d, m, num-1);
    cout<<"move disk "<<num<<" from tower "<<s<<" to tower "<<d<<endl;
    solve(m, s, d, num-1);
  }
}

struct hanoiInstance {
  int s;
  int m;
  int d;
  int num;
  int top;
  hanoiInstance(int src, int mid, int dest, int cnt, int t) : s(src), m(mid), d(dest), num(cnt), top(t) {};
};

void solveIterative(int s, int m, int d, int num) {
  stack<hanoiInstance*> stk;
  stk.push(new hanoiInstance(s, m, d, num, num));
  while(!stk.empty()) {
    hanoiInstance* inst = stk.top(); stk.pop();
    if(inst->num == 1)
      cout<<"move disk "<<inst->top<<" from tower "<<inst->s<<" to tower "<<inst->d<<endl;
    else {
      stk.push(new hanoiInstance(inst->m, inst->s, inst->d, inst->num - 1, inst->top-1));     
      stk.push(new hanoiInstance(inst->s, inst->m, inst->d, 1, inst->top));
      stk.push(new hanoiInstance(inst->s, inst->d, inst->m, inst->num - 1, inst->top-1));
    }
    delete inst;
  }
}

int main(){
  solveIterative(1,2,3, 4);
  return 0;
}
#include <vector>
#include <iostream>
using namespace std;
// max heap 
void heapify(vector<int>& heap) {
  int sz=heap.size();
  for(int i=sz-1; i>=1; i--) {
    int cur=i;
    while(2*cur < sz) {
      if(2*cur+1 < sz) {
        int larger=heap[2*cur] > heap[2*cur+1] ? 2*cur : 2*cur+1;
        if(heap[larger] > heap[cur]) {
          swap(heap[larger], heap[cur]);
          cur=larger;
        } else
          break;
      } else {
        if(heap[2*cur] > heap[cur]) {
          swap(heap[2*cur], heap[cur]);
          cur=2*cur;
        } else
          break;
      } 
    }
  }
}

int main() {
  vector<int> arr;
  for(int i=0; i<10; i++)
    arr.push_back(i);
  heapify(arr);
  for(int i=0; i<arr.size(); i++)
    cout<<arr[i]<<"  ";
}








class List {
  bool isElement;
  vector<List> list;
};

class HierList {
  stack<pair<*List, int>> level;
  bool hasNext() {
    if(level.empty())
      return false;
    pair<*List, int> cur = level.top();
    if(cur.second() == cur.first().size()-1) { // last element
      level.pop();
      return hasNext();
    }
    
#include <vector>
#include <stack>
#include <iostream>
using namespace std;

class Solution {
public:
  int areaM;
  void update(int a){
    areaM=max(areaM, a);
  }
    int largestRectangleArea(vector<int> &height) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(height.size()==0)
            return 0;
	  stack<int> st;
	  areaM = -1;
	  for(int i=0;i<height.size();i++) {
		if(st.empty())
		  st.push(i);
		else if(height[i] >= height[st.top()])
		  st.push(i);
		else {
		  while(!st.empty() && height[st.top()] > height[i]){
			int topIdx=st.top();
			st.pop();
            int left;
			if(st.empty())
			  left=0;
			else
			  left=st.top()+1;
			int area = height[topIdx] * (i - left);
			cout<<"removing...left="<<left<<"  i="<<i<<"  area="<<area<<endl;
			update(area);
		  }
		  st.push(i);
		}
	  }
	  while(!st.empty()) {
		int idx=st.top();
		st.pop();
		int left;
		if(st.empty())
		  left=0;
		else
		  left=st.top()+1;
		int area=height[idx]*(height.size() - left);
		update(area);
	  }
	  return areaM;
    }
};


int main(){
  int temp[10];
  temp[11]=1;
  cout<<temp[11]<<endl;
  vector<int> v;
  v.push_back(2);
  v.push_back(1);
  v.push_back(2);
  Solution s;
  cout<<s.largestRectangleArea(v)<<endl;
}
#include <iostream>
#include <vector>
#include <set>
using namespace std;

set<vector<vector<int> > > tried;

bool inRange(int a, int b) {
  return (a>=0 && a<3) && (b>=0 && b<3);
}
void compute(vector<vector<int> > table, int x, int y, vector<vector<int> >&target) {
  if(tried.find(table) != tried.end()) {
    cout<<"visited, reject"<<endl;
    return;
  }
  tried.insert(table);
  if(table == target) {
    cout<<"found"<<endl;
    return;
  }
  cout<<endl<<"try "<<x<<" "<<y<<endl;
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) 
      cout<<table[i][j]<<" ";
    cout<<endl;
  }
  int xx=x-1, yy=y;
  if(inRange(xx, yy)) {
    vector<vector<int> > t2 = table;
    swap(t2[x][y], t2[xx][yy]);
    compute(t2, xx, yy, target);
  }
  xx=x+1; yy=y;
 if(inRange(xx, yy)) {
    vector<vector<int> > t2 = table;
    swap(t2[x][y], t2[xx][yy]);
    compute(t2, xx, yy, target);
  }
 xx=x; yy=y-1;
 if(inRange(xx, yy)) {
   vector<vector<int> > t2 = table;
   swap(t2[x][y], t2[xx][yy]);
   compute(t2, xx, yy, target);
  }
 xx=x; yy=y+1;
 if(inRange(xx, yy)) {
   vector<vector<int> > t2 = table;
   swap(t2[x][y], t2[xx][yy]);
   compute(t2, xx, yy, target);
 }
}

int main() {
  int row1[]={1,2,3};
  int row2[]={4,5,6};
  int row3[]={0,7,8};
  vector<int> v1(row1, row1+3);
  vector<int> v2(row2, row2+3);
  vector<int> v3(row3, row3+3);
  vector<vector<int> > table;
  table.push_back(v1);
  table.push_back(v2);
  table.push_back(v3);
  int rowa[]={1,2,3};
  int rowb[]={4,5,6};
  int rowc[]={7,8,0};
  vector<int> va(rowa, rowa+3);
  vector<int> vb(rowb, rowb+3);
  vector<int> vc(rowc, rowc+3);
  vector<vector<int> > target;
  target.push_back(va);
  target.push_back(vb);
  target.push_back(vc);
  compute(table, 2, 1, target);
  return 0;
}
#include <iostream>
#include <stack>
using namespace std;
struct node{
  int val;
  node* left;
  node* right;
  node(int v) {left=right=NULL; val=v;}
};
void traverse(node *n) {
  node *cur=n;
  stack<node*> s;
  while(1) {
	while(cur){
	  s.push(cur);
	  cur=cur->left;
	}
	cout<<s.top()->val<<endl;
	cur=s.top()->right;
	s.pop();
	if(!cur && s.empty())
	  break;
  }
}

void postOrder(node *n) {
  node *cur=n;
  node *prev=NULL;
  stack<node *> s;
  s.push(n);
  while(!s.empty()) {
	cur=s.top();
	if((!cur->left && !cur->right)||(cur->right && prev==cur->right)||(cur->left && prev==cur->left)) {
	  cout<<cur->val<<endl; prev=cur; s.pop();
	} else {
	  if(cur->right)
		s.push(cur->right);
	  if(cur->left)
		s.push(cur->left);
	}
  }
}

int main() {
  node a(1), b(2), c(3), d(4), e(5);
  a.left=&b;
  a.right=&d;
  b.left=&c;
  b.right=&e;
  postOrder(&a);
}
class Solution {
public:
    bool isInterleave(string s1, string s2, string s3) {
        int len1=s1.length();
        int len2=s2.length();
        int len3=s3.length();
        if(len1+len2 != len3)
            return false;
        bool isInter[len1+1][len2+1];
        for(int i=0; i<len1+1; i++)
            for(int j=0; j<len2+1; j++)
                isInter[i][j]=false;
        isInter[0][0]=true;
        for(int i=1; i< len1+1; i++){
            if(s1[i-1]==s3[i-1])
                isInter[i][0]=true;
            else
                break;
        }
        for(int i=1; i< len2+1; i++){
            if(s2[i-1]==s3[i-1])
                isInter[0][i]=true;
            else
                break;
        }
        for(int i=1; i<len1+1; i++)
            for(int j=1; j<len2+1; j++){
                if((isInter[i][j-1] && s3[i+j-1]==s2[j-1]) || (isInter[i-1][j] && s3[i+j-1]==s1[i-1]))
                    isInter[i][j]=true;
            }
        return isInter[len1][len2];
    }
};
#include <set>
#include <string>
using namespace std;

class Solution {
public:
    int lengthOfLongestSubstring(string s) {
        // IMPORTANT: Please reset any member data you declared, as
        // the same Solution instance will be reused for each test case.
        set<char> cnt;
        if(s.length()<=0) return 0;
        int a=0, b=0;
		int gLen=0;
        while(b<s.length()-1) {
		  gLen=max(gLen, b-a+1);
		  b++;
		  char c=s[b];
		  if(cnt.find(c) == cnt.end()) {
			cnt.insert(c);
			continue;
		  } else {
			while(cnt.find(c)!=cnt.end()) {
			  cnt.erase(s[a]); a++;
			}
			cnt.insert(c);
		  }
		}
		gLen=max(gLen, b-a+1);
		return gLen;
    }
};

int main() {
  string s="hchzvfrkmlnozjk";
  Solution sol;
  sol.lengthOfLongestSubstring(s);
}
int maxConsSum2(const vector<int> &arr) {
  int gM=0;
  int curM=0;
  int sz=arr.size();
  for(int i=0;i<sz;i++) {
	if(curM>0)
	  curM=curM+arr[i];
	else
	  curM=arr[i];
	gM=max(gM, curM);
  }
  int curMi=0, gMi=0;
  for(int i=0; i<sz; i++) {
	if(curMi<0)
	  curMi = curMi+arr[i];
	else
	  curMi = arr[i];
	gMi=min(gMi, curMi);
  }
  int sum=0;
  for(int i=0; i<sz; i++)
	sum+=arr[i];
  int cirMax=sum-gMi;
  return max(cirMax, gM);
}
/*????
struct Point {double x; double y}
*/
bool eq(double x, double y) {
  return ((x-y)<0.0000001 && (x-y)>-0.0000001) ? true:false;
}

bool cmp(double x, double y){
  return (x<y) ? true:false;
}

int maxPointsOnLine(vector<Point> &points) {
  int sz=points.size();
  if(sz==0)
	return 0;
  if(sz==1)
	return 1;
  int globalMax=0;
  for(int i=0; i<sz; i++) {
	vector<double> slopes;
	int infiSlope=0;
	int maxC=0;
	for(int j=0; j<sz; j++) {
	  if(i==j)
		continue;
	  if(eq(points[j].x, points[i].x)){
		infiSlope++;
	  } else {
		double slope=(points[j].y-points[i].y)/(points[j].x-points[i].x);
		slopes.push_back(slope);
	  }
	}
	sort(slopes.begin(), slopes.end(), cmp);
	int b=0; int cur=1;
	for(b=1;b<slopes.size();b++) {
	  if(eq(slopes[b], slopes[b-1]))
		cur++;
	  else{
		maxC=max(maxC, cur);
		cur=1;
	  }
	}
	maxC=max(maxC, cur);
	maxC=max(maxC, infiSlope);
	globalMax=max(globalMax, maxC);
  }
  return globalMax;
}
int maxRectSum(vector<vector<int> > &matrix) {
  int row=matrix.size();
  int col=matrix[0].size();
  //int subSum[col][col][row];

  vector<vector<vector<int> > > subSum(col, vector<vector<int> >(col, vector<int>(row, 0)));

  for(int i=0; i<row; i++)
	for(int j=0;j<col; j++)
	  for(int k=j; k<col; k++) {
		if(j==k)
		  subSum[j][k][i] = matrix[i][j];
		else if(j==0)
		  subSum[j][k][i] = subSum[j][k-1][i] + matrix[i][k];
		else
		  subSum[j][k][i] = subSum[0][k][i] - subSum[0][j-1][i];
	  }
  int gMax=0;
  for(int c1=0;c1<col; c1++)
	for(int c2=c1; c2<col;c2++) {
	  int curMax = subSum[c1][c2][0];
	  gMax=max(gMax, curMax);
	  for(int r=1; r<row; r++) {
		if(curMax > 0)
		  curMax = curMax + subSum[c1][c2][r];
		else
		  curMax = subSum[c1][c2][r];
		gMax=max(gMax, curMax);
	  }
	}
  return gMax;
}
class Solution {
public:
  double getMed(int A[], int a, int b){
	if((b-a)%2==0)
	  return A[(a+b)/2];
	else
	  return (A[(a+b)/2]+A[(a+b)/2+1])/2.0;
  }
  double getMed2(int arr[], int a, int b, int val) {
	if((b-a)%2==0) {
	  int x=arr[(a+b)/2-1];
	  int y=arr[(a+b)/2];
	  int z=arr[(a+b)/2+1];
	  if(val <= x)
		return (x+y)/2;
	  else if(val >= z)
		return (y+z)/2;
	  else
		return (val+y)/2;
	} else {
	  int x=arr[(a+b)/2];
	  int y=arr[(a+b)/2+1];
	  if(val <= x)
		return x;
	  else if(val >= y)
		return y;
	  else
		return val;
	}
  }
  double find(int A[], int a1,int a2, int B[], int b1,int b2) {
	if(a1==a2 && b1==b2)
	  return (A[a1]+B[b1])/2;
	if(a2==a1) {
	  int aval=A[a1];
	  if(aval <= B[b1]) 
		return getMed(B, b1, b2-1);
	  else if(aval >= B[b2])
		return getMed(B, b1+1, b2);
	  else {//aval within B
		return getMed2(B, b1, b2, aval);
	  }
	} else if(b1==b2) {
	  int bval=B[b1];
	  if(bval <= A[a1])
		return getMed(A, a1, a2-1);
	  else if(bval >= A[a2])
		return getMed(A, a1+1, a2);
	  else {//aval within B
		return getMed2(A, a1, a2, bval);
	  }
	}
	int am=(a1+a2)/2;
	int bm=(b1+b2)/2;
	int cuta=a2-am;
	int cutb=b2-bm;
	if(A[am]>B[bm]) {
	  if(cuta>cutb)
		return find(A, a1, a2-cutb, B, b1+cutb, b2);
	  else
		return find(A, a1, a2-cuta, B, b1+cuta, b2);
	} else {
	  if(cuta>cutb)
		return find(A, a1+cutb, a2, B, b1, b2-cutb);
	  else
		return find(A, a1+cuta, a2, B, b1, b2-cuta);
	}
  }
    double findMedianSortedArrays(int A[], int m, int B[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  if(m==0)
		return getMed(B, 0, n-1);
	  else if(n==0)
		return getMed(A, 0, m-1);
	  return find(A, 0, m-1, B, 0, n-1);
    }
};









///////////////////////////
class Solution {
public:
  double getMed(int A[], int a, int b){
    if((b-a)%2==0)
	  return A[(a+b)/2];
	else
	  return (A[(a+b)/2]+A[(a+b)/2+1])/2.0;
  }
  double getMed2(int arr[], int a, int b, int val) {
	if((b-a)%2==0) {
	  int x=arr[(a+b)/2-1];
	  int y=arr[(a+b)/2];
	  int z=arr[(a+b)/2+1];
	  if(val <= x)
		return (x+y)/(2.0);
	  else if(val >= z)
		return (y+z)/(2.0);
	  else
		return (val+y)/(2.0);
	} else {
	  int x=arr[(a+b)/2];
	  int y=arr[(a+b)/2+1];
	  if(val <= x)
		return x;
	  else if(val >= y)
		return y;
	  else
		return val;
	}
  }
  double find(int A[], int a1,int a2, int B[], int b1,int b2) {
	if(a1==a2 && b1==b2)
	  return (A[a1]+B[b1])/(2.0);
	if(a2==a1) {
	  int aval=A[a1];
	  if(aval <= B[b1]) 
		return getMed(B, b1, b2-1);
	  else if(aval >= B[b2])
		return getMed(B, b1+1, b2);
	  else {//aval within B
		return getMed2(B, b1, b2, aval);
	  }
	} else if(b1==b2) {
	  int bval=B[b1];
	  if(bval <= A[a1])
		return getMed(A, a1, a2-1);
	  else if(bval >= A[a2])
		return getMed(A, a1+1, a2);
	  else {//aval within B
		return getMed2(A, a1, a2, bval);
	  }
	}
	int am=(a1+a2)/2;
	int bm=(b1+b2)/2;
	int cuta=a2-am;
	int cutb=b2-bm;
	if(A[am]>B[bm]) {
	  if(cuta>cutb)
		return find(A, a1, a2-cutb, B, b1+cutb, b2);
	  else
		return find(A, a1, a2-cuta, B, b1+cuta, b2);
	} else if(A[am]<B[bm]) {
	  if(cuta>cutb)
		return find(A, a1+cutb, a2, B, b1, b2-cutb);
	  else
		return find(A, a1+cuta, a2, B, b1, b2-cuta);
	} else { //
        if((a2-a1)%2 == 0 || (b2-b1)%2 == 0)
            return A[am];
        else {
                //both even length
                int x=A[am];
                int y=B[bm];
                int p=A[am+1];
                int q=B[bm+1];
                int z=(p>q) ? q:p;
                return (x+z)/2.0;
        }
	}
  }
    double findMedianSortedArrays(int A[], int m, int B[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  if(m==0)
		return getMed(B, 0, n-1);
	  else if(n==0)
		return getMed(A, 0, m-1);
	  return find(A, 0, m-1, B, 0, n-1);
    }
};
/**
 * Definition for binary tree with next pointer.
 * struct TreeLinkNode {
 *  int val;
 *  TreeLinkNode *left, *right, *next;
 *  TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
 * };
 */
class Solution {
public:
    void connect(TreeLinkNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  queue<TreeLinknode*> q;
	  if(!root)
		return;
	  q.push(root); q.push(NULL);
	  TreeLinkNode *prev=NULL;
	  while(!(q.size()==1 && q.back()==NULL)) {
		TreeLinkNode *out=q.front(); q.pop();
		if(prev)
		  prev->next=out;
		prev=out;
		if(out==NULL){
		  q.push(NULL);
		}
		else{
		  if(out->left)
			q.push(out->left);
		  if(out->right)
			q.push(out->right);
		}
	  }
	}
};
#include <iostream>
#include <string>
#include <stdio.h>
using namespace std;
string foo(string oct) {
  unsigned int num=0;
  for(int i=0; i<oct.length(); i++) {
    num = num<<3;
    printf("after shift %x\n", num);
    num = num | (oct[i]-'0');
    cout<<num<<endl;
    printf("%x\n", num);
  }
  cout<<num<<endl;
  string hex;
  while(num) {
    int mask=num & 15;
    char c;
    if(mask>9)
      c=mask-10+'a';
    else
      c=mask+'0';
    cout<<c<<endl;
    hex=c+hex;
    num = num>>4;
  }
  return hex;
}

int main() {
  string oct="321";
  std::cout<<foo(oct);
  return 0;
}




static const int signature = 0xDEADBEEF;


void *operator new(std::size_t size) throw(std::bad_alloc) {
  using namespace std;
  size_t realSize=size+2*sizeof(int);
  void* pMem=malloc(realSize);
  *((int*)pMem)=signature;
  *((int*)(((char*)pMem)+size+sizeof(int)))=signature;
  return ((int*)(((char*)pMem)+sizeof(int)));
}

int main() {
  char *c=new(3);
  return 0;
}
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
int main(){
  int fd;
  abort();
  assert(0);
}







#include <iostream>
#include <stdio.h>
int findLeftBound(int arr[], int a, int b, int t) {
  int candidate=-1;
  while(a<=b) {
    int mid=a+(b-a)/2;
    if(arr[mid]==t) {
      candidate=mid;
      b=mid-1;
    } else if(arr[mid] > t) {
      b=mid-1;
    } else {
      a=mid+1;
    }
  }
  return candidate;
}

int main() {
  int a[]={1,2,2,2,3,4,4,4,5,6,8,99};
  printf("%d\n", findLeftBound(a, 0, 11, 99));
 printf("%d\n", findLeftBound(a, 0, 11, 991));
 printf("%d\n", findLeftBound(a, 0, 11, 1));
 printf("%d\n", findLeftBound(a, 0, 11, 2));
 printf("%d\n", findLeftBound(a, 0, 11, 7));
  return 1;
}
class Solution {
public:
    ListNode *partition(ListNode *head, int x) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  ListNode *head2 = NULL;
	  ListNode *head1 = NULL;
	  ListNode *cur2=NULL;
	  ListNode *cur1=NULL;
	  ListNode *cur=head;
	  while(cur) {
		if(cur->val >= x) {
		  if(!head2){
			head2=cur;
			cur2=head2;
		  } else {
			cur2->next=cur;
			cur2=cur2->next;
		  }
		} else {
		  if(!head1) {
			head1=cur;
			cur1=head1;
		  } else {
			cur1->next=cur;
			cur1=cur1->next;
		  }
		}
		cur=cur->next;
	  }
	  if(cur1)
		cur1->next=head2;
	  if(cur2)
		cur2->next=NULL;
	  return cur1 ? head1:head2;
    }
};
#include <vector>
#include <iostream>
using namespace std;

class Solution {
public:
  void perm(vector<vector<int> > &res, vector<int> &cur, vector<bool> &v, vector<int> &num, int i){
	if(i==num.size())
	  res.push_back(cur);
	else {
	  for(int j=0; j<num.size(); j++) {
		if(v[j])
		  continue;
		if(j!=0 && num[j-1]==num[j] && !v[j-1])
		  continue;
		v[j]=true;
		cur[i]=num[j];
		perm(res, cur, v, num, i+1);
		v[j]=false;
	  }
	}
  }
    vector<vector<int> > permuteUnique(vector<int> &num) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  vector<bool> v(num.size(), false);
	  vector<int> cur=num;
	  vector<vector<int> > res;
	  perm(res, cur, v, num, 0);
	  return res;
    }
};




int main() {
  Solution s;
  vector<int> in;
  in.push_back(2);
  in.push_back(2);
  in.push_back(2);
  in.push_back(2);
  vector<vector<int> > res=s.permuteUnique(in);
  cout<<res.size()<<"  "<<res[0].size()<<endl;
}
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
/* Write COUNT copies of MESSAGE to STREAM, pausing for a second
between each. */
void writer (const char* message, int count, FILE* stream)
{
for (; count > 0; --count) {
/* Write the message to the stream, and send it off immediately. */
fprintf (stream, "%s\n", message);
/* fflush (stream); */
/* Snooze a while. */
sleep (1);
}
 fflush(stream);
}
/* Read random strings from the stream as long as possible. */
void reader (FILE* stream)
{
char buffer[1024];
/* Read until we hit the end of the stream. fgets reads until
either a newline or the end-of-file. */
while (!feof (stream)
&& !ferror (stream)
&& fgets (buffer, sizeof (buffer), stream) != NULL)
fputs (buffer, stdout);
}
int main ()
{
int fds[2];
pid_t pid;
/* Create a pipe. File descriptors for the two ends of the pipe are
placed in fds. */
pipe (fds);
/* Fork a child process. */
pid = fork ();
if (pid == (pid_t) 0) {
FILE* stream;
/* This is the child process. Close our copy of the write end of
the file descriptor. */
close (fds[1]);
/* Convert the read file descriptor to a FILE object, and read
from it. */
stream = fdopen (fds[0], "r");
reader (stream);
close (fds[0]);
}
else {
/* This is the parent process. */
FILE* stream;
/* Close our copy of the read end of the file descriptor. */
close (fds[0]);
/* Convert the write file descriptor to a FILE object, and write
to it. */
stream = fdopen (fds[1], "w");
writer ("Hello, world.", 5, stream);
close (fds[1]);
}
return 0;
}
#include <pthread.h>
#include <stdio.h>
/* Compute successive prime numbers (very inefficiently). Return the
Nth prime number, where N is the value pointed to by *ARG. */
int compute_prime (void* arg)
{
int candidate = 2;
int n = *((int*) arg);
while (1) {
int factor;
int is_prime = 1;
/* Test primality by successive division. */
for (factor = 2; factor < candidate; ++factor)
if (candidate % factor == 0) {
is_prime = 0;
break;
}
/* Is this the prime number we?re looking for? */
if (is_prime) {
if (--n == 0)
/* Return the desired prime number as the thread return value. */
return candidate;
}
++candidate;
}
return 0;
}
int main ()
{
pthread_t thread;
int which_prime = 5000;
int prime;
/* Start the computing thread, up to the 5,000th prime number. */
pthread_create (&thread, NULL, &compute_prime, &which_prime);
/* Do some other work here... */
/* Wait for the prime number thread to complete, and get the result. */
pthread_join (thread, prime);
/* Print the largest prime it computed. */
printf("The %dth prime number is %d.\n", which_prime, prime);
return 0;
}
#include <stdio.h>
void putlong(long num) {
  int sz=sizeof(num);
  char charArr[sz];
  long mask=0xff;
  for(int i=0; i < sz; i++) {
    char oneChar = (char)(num & mask);
    printf("%d\n", oneChar);
    num >>= 8;
    charArr[sz-i-1] = oneChar;
  }
  for(int i=sz-1; i >= 0; i--) {
    putchar((int)charArr[i]);
  }
  putchar('f');
}
int main() {
  putchar('a');
  long n=3243;
  putlong(n);
}
  
map<int, int> nextMap;
map<int, int> prevMap;
void init(int N) {
  nextMap.clear(); prevMap.clear();
  for(int i=0;i<N-1;i++)
	nextMap[i]=i+1;
  nextMap[N]=-1;
  for(int i=1; i<N; i++)
	prevMap[i]=i-1;
  prev[0]=-1;
}
void removeNum(int x) {
  if(nextMap.find(x) == nextMap.end())
	return;
  else if(nextMap[x]==-1 && prevMap[x]==-1) {
	nextMap.erase(x); prevMap.erase(x);
  } else if(nextMap[x]==-1) {
	int prevMapNum=prevMap[x];
	prevMap.erase(x);
	nextMap.erase(x);
	nextMap[prevMapNum]=-1;
  } else if(prevMap[x]==-1) {
	int nextMapNum=nextMap[x];
	prevMap.erase(x);
	nextMap.erase(x);
	prevMap[nextMapNum]=-1;
  } else {
	int prevMapNum=prevMap[x];
	int nextMapNum=nextMap[x];
	nextMap[prevMapNum]=nextMapNum;
	prevMap[nextMapNum]=prevMapNum;
	nextMap.erase(x);
	prevMap.erase(x);
  }
}
int query(int x) {
  if(nextMap.find(x) != nextMap.end())
	return nextMap[x];
  else
	return -1;
}




set<int> num;
void init(int N) {
  for(int i=0; i<N; i++)
	num.insert(N);
}


void removeNum(int x){
  if(set.find(x) != set.end())
	set.erase();
}
		  
		  


int query(int x){
  if(set.find(x) != set.end())
	return x;
  set<int>::iterator it = set.lower_bound(x);
  if(it != set.end())
	return *it;
  else
	return -1;
}
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;
bool comp(const pair<int,int> &a, const pair<int,int> &b) {
  return a.first > b.first;
}
vector<int> rearrange(vector<pair<int,int> > in){
  sort(in.begin(), in.end(), comp);
  cout<<"sort ok"<<endl;
  vector<int> res(in.size(), -1);
  for(int i=0; i<in.size(); i++) {
    int shortThanI = in[i].second;
    cout<<"this is "<<in[i].first<<" num shorter than me is "<<shortThanI<<endl;
    int slot=0;
    
    while(shortThanI) {
      if(res[slot] == -1) {
        slot++;
        shortThanI--;
      } else
        slot++;
    }
    while(res[slot]!=-1)
      slot++;
    cout<<"slot "<<slot<<" is "<<in[i].first<<endl;
    res[slot]=in[i].first;
  }
  return res;
}

int main() {
  vector<pair<int,int> > in;
  pair<int, int> a,b,c;
  a.first=3; a.second=0;
  in.push_back(a);
  b.first=2; b.second=0;
  in.push_back(b);
  c.first=1; c.second=0;
  in.push_back(c);


  rearrange(in);
  return 0;
}
    
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    TreeNode *f;
    TreeNode *s;
    TreeNode *fNext;
    bool fSet;
    int cur;
    TreeNode *prev;
    void tra(TreeNode *n) {
        if(n->left)
            tra(n->left);
        if(cur > n->val){
            if(!fSet){
                f=prev;
                fNext=n;
                fSet=true;
            } else
                s=n;
        }
        cur=n->val;
        prev=n;
        if(n->right)
            tra(n->right);
    }
    void recoverTree(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        cur=-99999;
        f=NULL;
        s=NULL;
        fNext=NULL;
        fSet=false;
        if(root)
            tra(root);
        if(s)
            swap(f->val, s->val);
        else
            swap(f->val, fNext->val);
    }
};
#include <vector>
struct A{
  std::vector<A> v;
  A() {v.resize(0);}
};
int main() {
  A a;
}
#include <iostream>
using namespace std;

 
 struct ListNode {
      int val;
      ListNode *next;
      ListNode(int x) : val(x), next(NULL) {}
  };
 
class Solution {
public:
    ListNode *reverse(ListNode *n, ListNode *& tail) {
        if(!n){
            tail=NULL; return NULL;
        }
        if(!n->next) {
            tail=n; return n;
        }
        ListNode *t;
        ListNode *head=reverse(n->next, t);
        t->next=n;
        tail=n;
        return head;
    }
    ListNode* merge(ListNode *a, ListNode *b) {
        if(!a && !b)
            return NULL;
        if(!b)
            return a;
        if(!a)
            return b;
        //a && b
        ListNode *sub=merge(a->next, b->next);
        a->next=b;
        b->next=sub;
        return a;
    }
    void reorderList(ListNode *head) {
        // IMPORTANT: Please reset any member data you declared, as
        // the same Solution instance will be reused for each test case.
        if(head==NULL || head->next==NULL)
            return;
        ListNode* slow=head, *fast=head, *prev=NULL;
        while(fast && fast->next) {
            fast=fast->next->next;
            prev=slow;
            slow=slow->next;
        }
        // now slow is the 2nd list's head
        if(fast) {
            prev=slow; slow=slow->next; prev->next=NULL;
        } else {
            prev->next=NULL;
        }
        ListNode *tail;
        
        ListNode *second = reverse(slow, tail);
        merge(head, second);
        
    }
};

int main() {
  ListNode a(1), b(2), c(3), d(4);
  a.next=&b;
  b.next=&c;
  c.next=&d;
  Solution s;
  s.reorderList(&a);
  ListNode *cur=&a;
  while(cur) {
	cout<<cur->val<<" ";
	cur=cur->next;
  }
}
#include <stdio.h>


unsigned reverse(unsigned x) {
  unsigned res=0;
  printf("size %d\n", sizeof(x));
  for(int i=0; i<sizeof(x) * 8; i++) {
    res = res << 1;
    res = res | (x & 1);
    x = x >> 1;
  }
  return res;
}

int main() {
  unsigned r = reverse(0x1323);
  printf("%x\n", r);
}
#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
// 3xabc 1x3xabc
// 332x4abc 2x31x2x4abc
char numstr[10];
bool digit(char c) {
  return c<='9' && c>='0';
}
string encode(string in) {
  cout<<endl<<in<<endl;
  int cur=1;
  string res;
  int cnt=1;
  in.insert(in.begin(), '.');
  while(cur < in.length()) {
    bool done=false;
    int start = cur;
    while(start < in.length() && digit(in[start]))
      start++;
    if(start < in.length() && in[start] == 'x') {
      done=true; // num+x case
      // first write the previous chars before nums
      if(cnt==1)
        res += in[cur-1];
      else {
        res += std::to_string(cnt);
        res += in[cur-1];
      }
      cnt=1;
      // now handle number
      cur = cur+1;
      while(cur <= start) {
        if(in[cur] == in[cur-1]) {
          cnt++;
        } else {
          res += to_string(cnt);
          res += 'x';
          res +=in[cur-1];
          cnt=1;
        }
        cur++;
      }
      cur=start; // start from 1st char after num
    }
    if(!done) {
      if(in[cur] == in[cur-1]) {
        cnt++;
      } else {
        if(cnt==1)
          res += in[cur-1];
        else {
          res += to_string(cnt);
          res += in[cur-1];
        }
        cnt=1;
      }
    }
    cur++;
  }
  // write the final chars
   if(cnt==1)
     res += in[cur-1];
   else {
     res += to_string(cnt);
     res += in[cur-1];
   }
   cout << res << endl;
   return res;
}
    
  
int main() {
  encode("abc");
  encode("1xb");
  encode("123xb");
  encode("113xb");
  encode("133xb");
  encode("13dsf3x2");
  encode("xb32b43xxx24");
  return 1;
}
#include <queue>
#include <iostream>
#include <exception>
class myex: public std::exception{
  char* what() {
    return "myEX";
  }
};
class runningAvg {
  std::queue<double> window;
  int size;
  double sum;
public:
  runningAvg(int sz) {
    if(sz<=0)
      throw;
    size=sz;
    sum=0;
  }
  double next(double in) {
    if(window.size() < size) {
      window.push(in);
      sum+=in;
      return sum/window.size();
    } else {
      double removed=window.front();
      window.pop();
      sum-=removed;
      window.push(in);
      sum+=in;
      return sum/size;
    }
  }
};

int main() {
  try{
    throw myex;
  } catch(exeption& e) {
    std::cout<<e.what();
  }
  runningAvg(0);
  runningAvg(1);
  return 1;
}
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    bool isSameTree(TreeNode *p, TreeNode *q) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  if(p==NULL && q==NULL)
		return true;
	  else if(p==NULL || q==NULL)
		return false;
	  else if(p->val != q->val)
		return false;
	  else
		return isSameTree(p->left, q->left) && isSameTree(p->right, q->right);
    }
};
#include <iostream>
using namespace std;
struct node{
  char item;
  node* left;
  node* right;
};
node* deserialBT(int &idx, char out[]) {
  char val;
  val=out[idx++];
  if(val=='#')
	return NULL;
  else {
	node *n=new node;
	n->item=val;
	n->left=deserialBT(idx, out);
	n->right=deserialBT(idx, out);
	return n;
  }
}

void print(node *n) {
  if(!n)
	cout<<"# ";
  else{
	cout<<n->item<<" ";
	print(n->left);
	print(n->right);
  }
}

int main() {
  char out[]={'1','2','3','#','#','#','4','5','#','#','6','#','#'};
  int idx=0;
  node* res = deserialBT(idx, out);
  print(res);
}




int foo(vector<char> in, set<char> target) {
  map<char, pair<int, int> > lenIdx; // pair<index, length>
  map<char, int> order; // fill this from target in order
  int minCover = -1;
  for(int i=0; i<in.length(); i++) {
    char c=in[i];
    if(target.find(c)) {
      if(order[c]==0) {
        lenIdx[c]={i, 1};
      } else {
        char prevChar = (map.lower_bound(c)-1)->first;
        if(order.find[prevChar]) {
          int prevIdx=lenIdx[prevChar].first;
          int prevLen=lenIdx[prevCHar].second;
          lenIdx[c]={i, i-prevIdx+prevLen};
          if(c == target[target.length()-1]) {
            minCover = min(minCover, lenIdx[c]);
          }
        }
      }
    }
  }
  return minCover;
}
#include <vector>

using namespace std;

class Solution {
public:
    int singleNumber(int A[], int n) {
        vector<int> cnt(32, 0);
        for(int i=0; i<n; i++) {
            for(int d=0; d<32; d++) {
                int mask=(1<<d);
                if(mask & A[i])
                    cnt[d]++;
            }
        }
        int res = 0;
        for(int s=0; s<32; s++) {
            if(cnt[s]%3 != 0)
                res = res | (1<<s);
        }
        return res;
    }
};

int main() {
  Solution s;
  int arr[]={1};
  s.singleNumber(arr, 1);
  return 0;
};
#include <stdio.h>
#include <iostream>
using namespace std;

template <typename OBJ>
class singleton {
private:
  singleton() {};
  //void operator=(S const&);
  //~singleton();
  static OBJ* instance;
  static const int a=23;
public:
  static OBJ* get() {
	if(instance)
	  return instance;
	else {
	  instance = new OBJ;
	  return instance;
	}
  }
  static OBJ& getRef() {
	if(!instance)
	  instance = new OBJ;
	return *instance;
  }
};

template<typename OBJ>
OBJ* singleton<OBJ>::instance = new OBJ;

class mySing{
public:
  int a;
  int p(){printf("fdsf\n");return 11;}
};

int main() {
  cout<<singleton<mySing>::get()->p();
  singleton<mySing>::get()->a=3;
  cout<<singleton<mySing>::get()->a;
  //cout<<sgt.get()<<endl;
  //cout<<sgt.getRef()<<endl;
  //printf("%p,   \n", sgt.get());
}
class Solution {
public:
    void sortColors(int A[], int n) {
	  int one=0;
	  int three=n-1;
	  int i=0;
	  while(i<=three) {
		if(A[i] == 2){
		  i++;
		  continue;
		} else if(A[i] == 1) {
		  swap(A[i], A[one++]);
		} else {
		  swap(A[i], A[three--]);
		}
	  }
    }
};
bool exists(vector<vector<int> > &matrix, int target) {
  int row=matrix.size();
  int col=matrix[0].size();
  int i=row-1, j=0;
  while(1) {
	if(i<0 || j>=col)
	  return false;
	if(matrix[i][j] == target)
	  return true;
	eles if(matrix[i][j] > target)
	  i--;
	else
	  j++;
  }
}
class Solution {
public:
    int sqrt(int x) {
	  if(x==0)
		return 0;
	  if(x==1)
		return 1;
	  int a=1, b=x;
	  long long mid=(b-a)/2+a;
	  while(!(mid*mid <= x && (mid+1)*(mid+1)>x)) {
		if(mid*mid > x){
		  b=mid;
		} else {
		  a=mid+1;
		}
		mid=(a+b)/2;
	  }
	  return mid;
    }
	  
};
class A{
public:
  struct B{
    int x,y;
  };
  double d;
};

main(){
  B b;
}
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
class Solution {
public:
    vector<int> con;
    vector<vector<int> > res;
    vector<int> concat(vector<int> a, vector<int> b) {
      cout<<"merging ";
	  for(int p=0;p<a.size();p++)
		cout<<a[p]<<"  ";
	  cout<<endl<<"and ";
	  for(int q=0;q<b.size();q++)
		cout<<b[q]<<"  ";
	  cout<<endl;
	  con=a;
        for(int i=0;i<b.size();i++)
            con.push_back(b[i]);
        return con;
    }
  void printVec(vector<int> v){
	for(int i=0;i<v.size();i++)
	  cout<<v[i]<<" ";
	cout<<endl;
  }
    void vecMult(vector<vector<int> > &a, vector<vector<int> > &b) {
	  cout<<"vector mult, a size="<<a.size()<<"    "<<endl;
	  for(int x=0; x<a.size(); x++)
		printVec(a[x]);
	  cout<<"vector mult, b size="<<b.size()<<"    "<<endl;
	  for(int x=0; x<b.size(); x++)
		printVec(b[x]);
        for(int i=0; i<a.size(); i++)
            for(int j=0; j<b.size(); j++)
                res.push_back(concat(a[i], b[j]));
    }
    int dupCount(vector<int> &S, int cur, int idx) {
        int count = 0;
        for(int t=idx; t<S.size(); t++){
            if(S[t] == cur)
                count++;
            else
                break;
        }
        return count;
    }
    void generateVec(int cur, int cnt, vector<vector<int> > &single) {
        vector<int> each;
        single.push_back(each);
        for(int i=1; i<=cnt; i++) {
            each.push_back(cur);
            single.push_back(each);
        }
    }
    vector<vector<int> > subsetsWithDup(vector<int> &S) {
        res.clear();
        vector<int> empty;
        res.push_back(empty);
        sort(S.begin(), S.end());
		for(int t=0;t<S.size();t++)
		  cout<<S[t]<<" ";
		cout<<endl;
        vector<int>::iterator Si=S.begin();
        int cur;
        int len = S.size();
        int i=0;
        int prev=S[0]-1;
        while( i != len ) {
            cur = S[i];
            i++;
            if(cur == prev)
                continue;
            else{
                vector<vector<int> > single;
                int cnt = dupCount(S, cur, i-1);
				cout<<"current val="<<cur<<" cnt="<<cnt<<endl;
                generateVec(cur, cnt, single);
				cout<<"single's size="<<single.size()<<endl;
                vector<vector<int> > a=res;
				res.clear();
				cout<<"a's size="<<a.size()<<endl;
                vecMult(a, single);
                prev=cur;
            }
        }
        return res;
    }
};

int main() {
  Solution s;
  vector<int> v;
  v.push_back(1);v.push_back(2);
  s.subsetsWithDup(v);
}
// empty cells are indicated by '.'
class Solution{
public:
  int getBox(int x, int y) {
    return x/3*3 + y/3 + 1;
  }
  bool num(char c) {
    return c>='0' && c<='9';
  }
  void next(int &x, int &y) {
    if(y==8)
      x++;
    else
      y++;
  }
  bool recFill(int x, int y, vector<vector<bool>> &row, vector<vector<bool>> &col, vector<vector<bool>> &box, vector<vector<char>> &board) {
    if(x==8 && y>8)
      return true;
    if(num(board[x][y])) {
      next(x, y);
      return recFill(x, y, row, col, box, board);
    }
    for(int t=1; t<=9; t++) {
      if(!row[x][t] && !col[y][t] && !box[getBox(x, y)][t]) {
        row[x][t] = true; col[y][t] = true; box[getBox(x, y)][t]=true;
        int xx=x, yy=y; next(xx, yy);
        bool res = recFill(xx, yy, row, col, box, board);
        if(res)
          return true;
        row[x][t] = false; col[y][t] = false; box[getBox(x, y)][t]=false;
      }
    }
  }
  void solveSudoku(vector<vector<char>> &board) {
    if(board.size() != 9 || board[0].size() != 9)
      return;
    vector<vector<bool>> row(9, vector(10, false)); // 1 ... 9
    vector<vector<bool>> col(9, vector(10, false));
    vector<vector<bool>> box(9, vector(10, false));
    for(int i=0; i<9; i++)
      for(int j=0; j<9; j++) {
        if(num(board[i][j])){
          row[i][board[i][j]-'0']=true;
          col[j][board[i][j]-'0']=true;
          box[getBox(i,j)][board[i][j]-'0']=true;
        }
      }
    recFill(0, 0, row, col, box, board);
  }
};
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
  int total;
  void sum(TreeNode *n, int cur) {
	cur=cur*10 + n->val;
	if(!n->left && !n->right)
	  total+=cur;
	else{
	  if(n->left)
		sum(n->left, cur);
	  if(n->right)
		sum(n->right, cur);
	}
  }
    int sumNumbers(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  total=0;
	  if(!root)
		return 0;
	  sum(root, 0);
      return total;
	}
};
class Solution {
public:
  bool isRegion(int x, int y, vector<vector<char> > &b) {
	//out of bound
	if(x==-1 || y==-1 || x==b.size() || y==b[0].size())
	  return false;
	if(b[x][y]=='C')
	  return true;
	if(b[x][y]=='X')
	  return true;
	else {
	  b[x][y]='C';
	  return isRegion(x-1,y,b) && isRegion(x+1,y,b) && isRegion(x,y-1,b) && isRegion(x,y+1,b);
	}
  }
  void uncolor(int x, int y, vector<vector<char> > &b) {
	if(x==-1 || y==-1 || x==b.size() || y==b[0].size())
	  return;
	if(b[x][y]=='X' || b[x][y]=='O')
	  return;
	if(b[x][y]=='C') {
	  b[x][y]='O';
	  uncolor(x-1,y,b);
	  uncolor(x+1,y,b);
	  uncolor(x,y-1,b);
	  uncolor(x,y+1,b);
	}
  }

  void find(int &x, int &y, vector<vector<char> > &b) {
	for(int i=0; i<b.size(); i++)
	  for(int j=0; j<b[0].size(); j++) {
		if(b[i][j]=='O')
		  if(isRegion(i,j,b)) {
			x=i; y=j;
			uncolor(i,j,b);
			return;
		  }
		  else
			uncolor(i,j,b);
	  }
	x=-1; y=-1;
	return;
  }

  void fill(int x, int y, vector<vector<char> > &b) {
	if(x==-1 || y==-1 || x==b.size() || y==b[0].size())
	  return;
	if(b[x][y]=='O') {
	  b[x][y]='X';
	  fill(x-1,y,b);
	  fill(x+1,y,b);
	  fill(x,y-1,b);
	  fill(x,y+1,b);
	}
  }



    void solve(vector<vector<char>> &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  int &x, &y;
	  x=0;y=0;
	  while(1) {
		find(x,y,board);
		if(x==-1 && y==-1)
		  return;
		else
		  fill(x,y,board);
	  }
    }
};


//  faster

class Solution {
public:
void color(int x, int y, vector<vector<char> > &b) {
	if(x==-1 || y==-1 || x==b.size() || y==b[0].size())
	  return;
	if(b[x][y]=='X' || b[x][y]=='C')
	  return;
	if(b[x][y]=='O') {
	  b[x][y]='C';
	  uncolor(x-1,y,b);
	  uncolor(x+1,y,b);
	  uncolor(x,y-1,b);
	  uncolor(x,y+1,b);
	}
  }

   void solve(vector<vector<char>> &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	 int i, j;
	 for(i=0;i<board.size();i++)
	   color(i,0,board);
	 for(i=0;i<board.size();i++)
	   color(i,board[0].size()-1,board);
	 for(i=0;i<board[0].size();i++)
	   color(0,i,board);
	 for(i=0;i<board[0].size();i++)
	   color(board.size()-1,i,board);


	 for(i=0;i<board.size();i++)
	   for(j=0;j<board[0].size();j++)
		 if(board[i][j]=='O')
		   board[i][j]='X';

	 for(i=0;i<board.size();i++)
	   for(j=0;j<board[0].size();j++)
		 if(board[i][j]=='C')
		   board[i][j]='O';
	 


    }
};
#include <stdio.h>
#include <iostream>
using namespace std;
void sw(int &a, int &b) {
  int temp=a;
  a=b;
  b=temp;
}

template<class T>
void SW(T &a, T&b) {
  T temp=a;
  a=b;
  b=temp;
}

class test{
  int p;
  void get(){cout<<p;}
};

int main() {
  int x=1;
  int y=2;
  SW(x,y);
  cout<<"x"<<x<<"Y"<<y<<endl;
  long t=324;
  long s=234;
  SW(t,s);
  cout<<t<<" "<<s<<endl;
  test t1;
  test t2;
  SW(t1,t2);
};
class Solution {
public:
  bool isSym(TreeNode *p, TreeNode *q) {
	if(p==NULL && q==NULL)
	  return true;
	  else if(p==NULL || q==NULL)
		return false;
	  else if(p->val != q->val)
		return false;
	  else
		return isSym(p->left, q->right) && isSym(p->right, q->left);  
  }
  bool isSymmetric(TreeNode *root) {
	// Start typing your C/C++ solution below
	// DO NOT write int main() function
	if(!root)
	  return true;
	return isSym(root->left, root->right);
  }
};



// iterative
class Solution {
public:
  bool isSym(TreeNode *root) {
	queue<TreeNode*> ps, qs;
	ps.push(root->left);
	qs.push(root->right);
	while(!ps.empty() && !qs.empty()){
	  TreeNode *p=ps.front();
	  ps.pop();
	  TreeNode *q=qs.front();
	  qs.pop();
	  if(p==NULL && q==NULL)
		continue;
	  else if(p==NULL || q==NULL)
		return false;
	  else if(p->val != q->val)
		return false;
	  else{
		ps.push(p->left);
		ps.push(p->right);
		qs.push(q->right);
		qs.push(q->left);
	  }
	}
  }
  bool isSymmetric(TreeNode *root) {
	// Start typing your C/C++ solution below
	// DO NOT write int main() function
	if(!root)
	  return true;
	isSym(root);
  }
};

#include <string>
#include <vector>
#include <iostream>
using namespace std;


class Solution {
public:
    vector<string> fullJustify(vector<string> &words, int L) {
        // Note: The Solution object is instantiated only once and is reused by each test case.
	  vector<string> curLine;
	  vector<string> res;
	  int curLen=0;
	  for(int i=0; i<words.size(); i++) {
		if(curLen+words[i].length() + curLine.size() <= L) {
		  curLine.push_back(words[i]);
		  curLen+=words[i].length();
		} else {//not enough L
		  i--;
		  int spaces=L-curLen;
		  int wordCnt=curLine.size();
		  string line;
		  if(wordCnt==1) {
			line+=curLine[0];
			for(int cc=0; cc<spaces; cc++)
			  line+=' ';
		  } else {
			int gap=spaces/(wordCnt-1);
			int mod=spaces % (wordCnt-1);
			line=line+curLine[0];
			for(int s=0; s<gap+mod; s++)
			  line.append(1,' ');
			for(int j=1; j<curLine.size()-1; j++) {
			  line=line+curLine[j];
			  for(int ss=0; ss<gap; ss++)
				line=line+' ';
			}
			line=line+curLine[curLine.size()-1];
		  }
		  res.push_back(line);
		  curLine.clear();
		  curLen=0;
		}
	  }
	  //  last line
	  string lastLine;
	  for(int c=0; c<curLine.size()-1; c++)
		lastLine=lastLine + curLine[c] + ' ';
	  lastLine+=curLine[curLine.size()-1];
	  for(int c=0; c< L - curLen -(curLine.size()-1); c++)
		lastLine+=' ';
	  res.push_back(lastLine);
	  return res;
    }
};


int main() {
  Solution s;
  vector<string> words;
  words.push_back(string("Don't"));
 words.push_back(string("go"));
 words.push_back(string("around"));
 words.push_back(string("saying"));
 words.push_back(string("the"));
 words.push_back(string("world"));
 words.push_back(string("owes"));  
 words.push_back(string("you"));
 words.push_back(string("a"));
 words.push_back(string("living;"));
 words.push_back(string("the"));
 words.push_back(string("world"));
 words.push_back(string("owes"));
 words.push_back(string("you"));
 words.push_back(string("nothing;"));
 words.push_back(string("it"));
 words.push_back(string("was"));
 words.push_back(string("here"));
 words.push_back(string("first."));




vector<string> res = s.fullJustify(words, 30);
  for(int i=0; i<res.size(); i++)
	cout<<res[i]<<"--END"<<endl;
}
				  
#include <pthread.h>
#include <stdio.h>
/* Prints x?s to stderr. The parameter is unused. Does not return. */
pthread_mutex_t m;
pthread_cond_t cond;
bool var=false;
void* print_xs (void* unused)
{
  printf("thread\n");
  pthread_mutex_lock(&m);
  var=false;
  pthread_mutex_unlock(&m);
  printf("processing\n");

  pthread_mutex_lock(&m);
  var=true;  
  pthread_cond_signal(&cond);
  pthread_mutex_unlock(&m);
  printf("finished signal\n");
  return NULL;
}
/* The main program. */
int main ()
{
  pthread_t thread_id;
  pthread_mutex_destroy(&m);
  perror("destroy");
  if(pthread_mutex_lock(&m) != 0)
    perror("aa");
  pthread_mutex_init(&m, NULL);
  pthread_cond_init(&cond, NULL);
  /* Create a new thread. The new thread will run the print_xs
     function. */

  pthread_create(&thread_id, NULL, &print_xs, NULL);
  /* Print o?s continuously to stderr. */
  pthread_mutex_lock(&m);
  if(pthread_mutex_trylock(&m) != 0) {
    perror("error");
  }
  printf("waiting...\n");
  if(var==false)
    pthread_cond_wait(&cond, &m);
  printf("woke up\n");
  pthread_mutex_unlock(&m);
  

  printf("done\n");
  return 0;
}
#include <iostream>
#include <map>
using namespace std;

template <class Obj>
class Union {
public:
  // return group leader
  Obj* find(Obj* o) {
	if(parent.find(o) != parent.end()){
	  while(parent[o] != nullObj)
		o=parent[o];
	  return o;
	} else {
	  return nullObj;
	}
  }
  
  //return new group leader
  bool unite(Obj* a, Obj* b) {
	if(parent.find(a)==parent.end() || parent.find(b)==parent.end())
	  return false;
	if(find(a) == find(b))
	  return false;
	parent[find(a)]=find(b);
	unionCnt--;
	return true;
  }
  
  bool makeGroup(Obj* o) {
	if(parent.find(o) != parent.end()){
	  return false; // already exist
	} else {
	  parent[o]=nullObj;
	  unionCnt++;
	  return true;
	}
  }
  
  Union () {
	unionCnt=0;parent.clear(); nullObj=NULL;
  }

  int count() { return unionCnt;}
private:
  int unionCnt;
  map<Obj*, Obj*> parent;
  Obj* nullObj;
};

class node {
  int x;
};

int main() {
  Union<node> u;
  node a,b;
  u.makeGroup(&a);
  u.makeGroup(&b);
  u.find(&a);
  u.unite(&a,&b);
  u.unite(&a,&b);
  cout<<u.count();
  return 0;
}
class Solution {
public:
    int numTrees(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(n==0 || n==1)
            return 1;
        int count[n+1];
        count[0]=count[1]=1;
        for(int i=2;i<=n;i++){
            int sum=0;
            for(int j=1; j<=i; j++){
                sum+=count[j-1] * count[i-j];
            }
            count[i]=sum;
        }
        return count[n];
        
    }
};
class Solution {
public:
  vector<TreeNode *>* genRange(int a, int b) {
	vector<TreeNode *> *res=new vector<TreeNode *>;
	if(a > b){
	  res->push_back(NULL);
	  return res;
	} else if(a == b){
	  TreeNode *tn = new TreeNode(a);
	  res->push_back(tn);
	  return res;
	} else {
	  for(int i=a; i<=b; i++){
		vector<TreeNode*> *left=genRange(a, i-1);
		vector<TreeNode*> *right=genRange(i+1, b);
		for(vector<TreeNode*>::iterator li=left->begin(); li != left->end(); li++)
		  for(vector<TreeNode*>::iterator ri=right->begin(); ri != right->end(); ri++){
			TreeNode *tn = new TreeNode(i);
			tn->left=(*li);
			tn->right=(*ri);
			res->push_back(tn);
		  }
	  }
	  return res;
	}
  }

    vector<TreeNode *> generateTrees(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  return *genRange(1, n);
    }
};
#include <iostream>
using namespace std;
int print(int i){
  return i+2;
}
struct p{
  int val;
};
int main(){
  p* pp;
  print(pp->val);
  pp=new p;
  print(pp->val);
  delete pp;
  print(pp->val);
}
class Solution {
public:
  bool valid(TreeNode *n, int &mi, int &ma) {
	int curVal = n->val;
	int leftMax, leftMin, rightMin, rightMax;
	leftMax=curVal;
	leftMin=curVal;
	rightMax=curVal;
	rightMin=curVal;
	if(n->left) {
	  if(!valid(n->left, leftMin, leftMax))
		return false;
	  if(leftMax > curVal)
		return false;
	}
	if(n->right) {
	  if(!valid(n->right, rightMin, rightMax))
		return false;
	  if(rightMin < curVal)
		return false;
	}
	mi=leftMin;
	ma=rightMax;
	return true;
  }
	  
    bool isValidBST(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  int mi,ma;
	  if(!root)
		return true;
	  return valid(root, mi, ma);
    }
};
#include <iostream>
#include <vector>
using namespace std;
#if 0
class Solution {
public:
    bool isValidSudoku(vector<vector<char> > &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<bool> > rows(9, vector<bool>(9, false));
        vector<vector<bool> > cols(9, vector<bool>(9, false));
        vector<vector<bool> > blocks(9, vector<bool>(9, false));

        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                if (board[i][j] == '.') continue;
                int c = board[i][j] - '1';
                if (rows[i][c] || cols[j][c] || blocks[i - i % 3 + j / 3][c])
                    return false;
                rows[i][c] = cols[j][c] = blocks[i - i % 3 + j / 3][c] = true;
            }
        }
        return true;
    }
};
#endif


class Solution {
public:
  void next(int x, int y, int &xx, int &yy){
	if(x==8 && y==8){
	  xx=-1; yy=-1;
	} else if(y==8) {
	  xx=x+1; yy=0;
	} else {
	  xx=x; yy=y+1;
	}
  }
  bool valid(int x, int y, vector<vector<char> > &board, bool row[9][9], bool col[9][9], bool box[9][9]) {
	if(x==-1 && y==-1)
	  return true;
	bool thisRes=false;
	if(board[x][y] != '.'){
	    int xx, yy;
	  next(x, y, xx, yy);
	  thisRes=valid(xx, yy, board, row, col, box);
	} else { // to be filled
	  for(int i=0; i<9; i++) {
		if(row[x][i]==false && col[y][i]==false && box[x-x%3 + y/3][i]==false){
		  row[x][i] = col[y][i] = box[x-x%3 + y/3][i] = true;
		  int xx, yy;
		  next(x,y,xx,yy);
		  thisRes=valid(xx, yy, board, row, col, box);
		  if(thisRes)   return true;
		  row[x][i] = col[y][i] = box[x-x%3 + y/3][i] = false;
		}
	  }
	}
	return thisRes;

  }
    bool isValidSudoku(vector<vector<char> > &board) {
	  vector<vector<char> > myBoard=board;
	  // fill the hashes
	  bool row[9][9], col[9][9], box[9][9];
	  for(int i=0; i<9; i++)
		for(int j=0; j<9; j++) {
		  row[i][j]=col[i][j]=box[i][j]=false;
		}
	  for(int i=0; i<9; i++)
		for(int j=0; j<9; j++) {
		  if(board[i][j] != '.') {
			int val=board[i][j]-'1';
			row[i][val]=col[j][val]=box[i-i%3+j/3][val]=true;
		  }
		  return valid(0,0, board, row, col, box);
		}
	}
};


int main() {
  vector<vector<char> > b(9, vector<char>(9, '.'));
  b[0][0]=b[0][1]='1';
  //b[8][0]='9';
  Solution s;
  bool res = s.isValidSudoku(b);
  if(res)
	cout<<"valid"<<endl;
  return 0;
}
class Solution {
public:
  unordered_set<string> visited;
  unordered_set<string> bound;
  unordered_set<string> tempBound;
  int wordSize;
  
  void getNeigh(string s, unordered_set<string> &dict) {
	visited.insert(s);
	for(int i=0; i < wordSize; i++)
	  for(char c='a'; c<='z'; c++){
		string oneMove=s;
		oneMove[i]=c;
		if(dict.find(oneMove) != dict.end() && visited.find(oneMove) == visited.end() && oneMove != s)
		  tempBound.insert(oneMove);
	  }
  }

  int explore(string end, int cnt, unordered_set<string> &dict) {
	for(unordered_set<string>::iterator it=bound.begin(); it != bound.end(); it++) {
	  //put neighbour into tempBound;
	  if(*it == end)
		return 0;
	  getNeigh(*it, dict);
	}
	if(tempBound.size() == 0)
	  return -1;
	bound.clear();
	bound=tempBound;
	tempBound.clear();
	return 1;
  }
  int ladderLength(string start, string end, unordered_set<string> &dict) {
	// Start typing your C/C++ solution below
	// DO NOT write int main() function
	wordSize=end.length();
	visited.clear();
	bound.clear();
	tempBound.clear();
	bound.insert(start);
	int cnt=1;
	while(1) {
	  
	  int found=explore(end, cnt, dict);
	  if(found==0)
		break;
	  else if(found==-1)
		return 0;
	  cnt++;
	  
	}
	return cnt;
  }
};




class Solution {
public:
  void rec(vector<vector<string> > &res, int wid, int step, string end, vector<string> vec, unordered_set<string> &dict){
	if(step==0){
	  if(vec.back()==end)
		res.push_back(vec);
	  return;
	}
	for(int i=0; i<wid; i++){
	  for(int j=0; j<26; j++){
		char c='a'+j;
		string t=vec.back();
		t[i]=c;
		if(dict.find(t) != dict.end()){
		  vec.push_back(t);
		  rec(res, wid, step-1, end, vec, dict);
		  vec.pop_back();
		}
	  }
	}
  }
    int findLadderss(string start, string end, unordered_set<string> &dict) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
	  int step=1;
	  queue<string> q;
	  unordered_set<string> visited;
	  int level=1;
	  int wid=start.length();
	  q.push(start);
	  bool found=false;
	  q.push(string(""));
	  while(!q.empty()){
		string s=q.front();
		q.pop();
		if(s == "") {
		  if(q.empty())
			break;
		  q.push(string(""));
		  level++;
		  continue;
		}
		visited.insert(s);
		for(int i=0; i<wid; i++){
		  for(int j=0; j<26; j++){
			char c='a'+j;
			string t=s;
			t[i]=c;
			if(t==end) {
			  found=true;
			  return level;
			}
			else if(dict.find(t) != dict.end() && visited.find(t) == visited.end() && s != t)
			  q.push(t);
		  }
		}
	  }
	  if(!found)
		return -1;
	}
  vector<vector<string>> findLadders(string start, string end, unordered_set<string> &dict) {
	int level=findLadderss( start,  end,dict); 
	  vector<vector<string> > res;
	  if(level==-1)
		return res;
	  vector<string> vec;
	  vec.push_back(start);
	  rec(res, end.length(), level, end, vec, dict);
	  return res;
	  
    }
};
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
int main ()
{
pid_t child_pid;
/* Create a child process. */
child_pid = fork ();
if (child_pid > 0) {
/* This is the parent process. Sleep for a minute. */
  exit(0);
}
else {
/* This is the child process. Exit immediately. */
sleep (60);
exit (0);
}
return 0;
}
#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
using namespace std;


class Solution {
public:
    vector<vector<int> > threeSum(vector<int> &num) {
        vector<vector<int> > res;
        set<vector<int> > found;
        if(num.size() < 3) return res;
        sort(num.begin(), num.end());
        for(int i=0; i<num.size()-3; i++) {
            if(i>0 && num[i]==num[i-1])
                continue;
            int j=i+1, k=num.size()-1;
            int target=0-num[i];
            cout<<"start  "<<j<<"  "<<k<<endl;
            while(j<k) {
                if(num[j] + num[k] < target)
                    j++;
                else if(num[j] + num[k] > target)
                    k--;
                else {
                  cout<<j<<"  "<<k<<endl;
                    vector<int> ans;
                    ans.push_back(num[i]);
                    ans.push_back(num[j]);
                    ans.push_back(num[k]);
                    if(found.find(ans)==found.end()) {
                        found.insert(ans);
                        res.push_back(ans);
                    }
                    j++; k--;
                }
            }
        }
        return res;
    }
};

int main() {
  Solution s;
  vector<int> v;
  v.push_back(0);
  v.push_back(0);
  v.push_back(0);
  s.threeSum(v);
  return 0;
}
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
	void copyArr(vector<vector<int> > & arr, vector<vector<int> > & newArr, int N){
		for(int i=0;i<N;i++){
			vector<int> row;
			newArr.push_back(row);
				for(int j=0;j<N;j++)
					newArr[i].push_back(arr[i][j]);
			}
		}
	void colorArr(vector<vector<int> > & arr, int n, int i, int N){
		arr[n][i]=2;
		for(int x=0;x<N;x++)
			if(arr[n][x]!=2)
				arr[n][x]=1;
		for(int x=0;x<N;x++)
			if(arr[x][i]!=2)
				arr[x][i]=1;
		// diag
		for(int x=n,  y=i;x<N && y<N;x++,y++)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x>=0 && y<N;x--,y++)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x>=0 && y>=0;x--,y--)
			if(arr[x][y]!=2)
				arr[x][y]=1;
		for(int x=n,  y=i;x<N && y>=0;x++,y--)
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
	cout<<"input the number of n for N Queen problem"<<endl;
	int x;
	cin>>x;
	s.solveNQueens(x);
}


// getFriendsListForUser
// getPurchasesForUser

struct productCnt {
	string product;
	int cnt;
};
bool comp(productCnt &a, product &b) {
	return a.cnt > b.cnt ;
}
list<string> generateRank(string user) {
	list<string> friends = getFriendsListForUser(user);
	list<string> recommendation;
	if(friends.size() == 0) {
		// no friends, no rec
		return recommendation;
	}
	map<string, int> productCount;
	for(int i=0; i<friends.size(); i++) {
		list<string> purchaseList = getPurchasesForUser(friends[i]);
		for(int p=0; p<purchaseList; p++) {
			if(map.find(purchaseList[p]) == map.end())
				productCount[purchaseList[p]] = 1;
			else
				productCount[purchaseList[p]]++;
		}
	}
	vector<productCnt> v; // used to store purchase-count pair and do sorting
	list<string> myList = getPurchasesForUser(user);
	set<string> myOwnPurchase;
	for(auto ii = myList.begin();ii != myList.end(); ii++) {
		myOwnPurchase.insert(*ii);
	}
	for(auto mi = productCount.begin(); mi != productCount.end(); mi++) {
		if(myOwnPurchase.find(mi->first) != myOwnPurchase.end())
			// user already bought this
			continue;
		productCnt pc;
		pc.product = mi->first;
		pc.cnt = mi->second;
		v.push_back(pc);
	}
	// sort in reverse order with comp()
	sort(v.begin(), v.end(), comp);
	for(auto vi = v.begin(); vi != v.end(); vi++) {
		recommendation.push_back(vi->first);
	}
	return recommendation;
}














Rocket Fuel
code challenge: auto racer
电面1:
第一题：贪心
Given a number, can you remove k digits from the number so that the new
formatted number is smallest possible.
input: n = 1432219, k = 3
output: 1219

vector<int> removeKDigit(vector<int> num, int k) {
	if(k==0)
		return num;
	int m=INT_MAX;
	int mi=0;
	for(int i=0; i<=k; i++) {
		if(num[i]<m) {
			mi=i;
			m=num[i];
		}
	}
	for(int i=0; i<=mi; i++)
		num.erase(num.begin());
	vector<int> subRes=removeKDigit(num, k-mi);
	subRes.insert(subRes.begin(), m);
	return subRes;
}
		
		
第二题：DP
BT(binary tree), want to find the LIS(largest independent set) of the BT
LIS: if the current node is in the set, then its chilren should not be in
the set. So that the set has the largest number of nodes.
map<node*, int> size;
int LIS(node* n) {
	if(!n) return 0;
	if(size.find(n) != size.end())
		return size[n];
	int useCur=1;
	useCur+=(n->left ? LIS(n->left->left)+LIS(n->left->right) : 0) + (n->right ? LIS(n->right->left)+LIS(n->right->right)) : 0);
	int dontUseCur=0;
	dontUseCur+=LIS(n->left)+LIS(n->right);
	int res=max(useCur, dontUseCur);
	size[n]=res;
	return res;
}

电面2：
第一题：Median of Two Sorted Arrays
第二题：DP
一个二维数组，元素是0或1，找出最大的由1构成的"X"形状

onsite:
1. print all subsets
   system design(N topics, publishers, subscribers, scalability, distributed)
   the most frequent urls in the past minute, hour, day
2. manager interview
   code review
3. shortest path between two nodes of a tree(no parent pointer)
4. machine learning(不懂)
5. machine learning(不懂)

Rocket Fuel是自己投的，因为在网上看到code challenge挺有意思。onsite的时候了
解到他家最近要搬进新楼里，应该招人很多，大家可以试一试，题目不简单

Google:
电面：
remove duplicate lines of a file(what if the file is very large which could
not be held in the main memory)
开关灯问题
Trapping Rain Water(leetcode)
sometimes a program works, sometimes it does not. Possible reasons
multithreading, multiprocess, different input, different memory access(cache hit/miss), rely on random numbers


onsite:
1. clone directed graph(recursive, non-recursive)
   longest common suffix of two linked list
//first push both to two stacks, then pop and compare.
vector<int> longestCommon(node* a, node* b) {
	stack<node> sa, sb;
	while(a) {
		sa.push(a);
		a=a->next;
	}
	while(b) {
		sb.push(b);
		b=b->next;
	}
	vector<int> res;
	while(!sa.empty() && !sb.empty()) {
		if(sa.top()->val==sb.top()->val){
			res.push_back(sa.top());
			sa.pop();
			sb.pop();
		} else
			break;
	}
	reverse(res.begin(), res.end());
	return res;
}


   data structure design
2. how many (m, n) pairs such that m*m+n*n<N
int count(int N) {
	if(N==0) return 0;
	int m=0, n=sqrt(N);
	if(n*n==N)
		n--;
	int cnt=0;
	while(n>=0 && m*m+n*n < N) {
		cnt=cnt+n+1;
		cout<<n<<"  "<<m<<endl;
		m++;
		bool negN=false;
		while(m*m+n*n >= N) {
			n--;
			if(n<0) {
				negN=true;
				break;
			}
		}
		if(negN)
			break;
	}
	return cnt;
}

// solution 2
int count(int N) {
	int m=0, n=sqrt(N);
	int cnt=0;
	while( m*m < N ) {
		while(m*m + n*n >= N)
			n--;
		cnt+=n+1;
		cout<<m<<" "<<n<<endl;
		m++;
	}
}

   线索化二叉树
3. 判断一个点是否在一个凸多边形内, O(n), O(logn)
4. group items(BFS)
   MapReduce(filter a collection of documents, the words which occur more
than 5000 times)

google面的不好，因为实在是太累了，幸运的是还是给offer了。

linkedin
电面1：
第一题：给一个words list, 输入两个单词，找出这两个单词在list中的最近距离(先
写了一个没有预处理的，又写了一个预处理建index的)
['green', 'blue', 'orange', 'purple', 'green']  f.distance(list, 'blue', '
green') # output 1
第二题：类似binary search的一题，要注意boundary case

电面2：
binary tree level order traversal, 写了三种方法。。。(BFS用arraylist，类似
DFS，BFS用queue)

onsite:
1. romanToInt, intToRoman,
   N points, m nearst ones
2. 双向链表，每个node可能都有父节点和子节点，每个父子节点又是一个链表。把它
拍扁，顺序随意，O(1)空间复杂度
//assume each parent/child node could only be the head of a sub list
void toEnd(node* &end) {
	while(end->next)
		end=end->next;
}
node* flatList(node* head) {
	if(!head) return NULL;
	node *end=head, *cur=head;
	toEnd(end);
	while(cur) {
		if(cur->parent) {
			end->next=cur->parent;
			cur->parent->child=NULL;
			cur->parent=NULL;
			toEnd(end);
		} else if(cur->child) {
			end->next=cur->child;
			cur->child->parent=NULL;
			cur->child=NULL;
			toEnd(end);
		} else {
			cur=cur->next;
		}
	}
	return head;
}

   edit distance
3. system deisign: design amazon product page
product id, name, manufacture, price, inventory.
description, review.

4. project presentation
5. group fit


FLAG干货

Linkedin

phone1：烙印
lowest common ancestor w/ and w/o parent pointer

//w/o
node *anc = NULL;
bool find(node* n, int v1, int v2) {
	if(!n) return false;
	bool left=find(n->left, v1, v2);
	bool right=find(n->right, v1, v2);
	if(left && right) {
		anc=n;
		return true;
	} else if(left || right) {
		if(n->val == v1 || n->val == v2) {
			anc=n;
			return true;
		} else
			return true;
	} else if(n->val == v1 || n->val == v2) 
		return true;
	else
		return false;
}

phone2：国人
search in rotated sorted array
leetcode has both version

onsite:
1.两个国人
implement addInterval(start, end) and getcoverage(),
2.两个国人
talk projects and some behavior question
3.烙印
lunch, talk about technologies interest
4.亚裔，不确定是否国人
Manager, talked a lot of behavior questions, interest and projects
5.烙印
Design: tinyurl
6.烙印+小白
1.exclusive array, give an arr1, return a new arr2, arr2[i] is the
multiplication of all elements in arr1 except arr1[i]
2.boolean isMirrorTrees(tree1, tree2)/inplace convert a tree to its mirror
tree/create a new mirror tree
3.find the intersection of two linked list(do not use hashmap)


Amazon

phone1: 烙印
Given a list of test results (each with a test date, Student ID, and the
student’s Score), return the Final Score for each student. A student’s
Final Score is calculated as the average of his/her 5 highest test scores.
You can assume each student has at least 5 test scores.

phone2:白男
1.大数plusOne
2.给你一个按字母顺序排好的字典（但你不知道字母顺序,非英语），要求找出字母顺序
例：
单词顺序：
    wrt
    wrf
    er
    ett
    rftt
字母顺序：
    w,e,r,t,f

topo sort
	
onsite:
1. 白男
class MagicNumber{
boolean isMagicNum(long num);
long nextMagic(long num){
    while(!isMagicNum(num)){
    num++;
}
return num;
}
}
consider a data structure to improve the nextMagic(long num)
2. 烙印
behavior questions and text editor design(insert, add, search, cut, paste)
3. 白男
大数加法 (int数组表示大数，每一个元素代表一个2^31进制数字)
4. 日本人manager
lunch interview:
    4.1 describe a time you are stressful to meet a deadline
    4.2 describe a time you feel most proud in your professional career
    4.3 what would you change in your past project if you have a chance
5. 白男
give API: List<Movie> getMovies(Actor a);  List<Actor> getActors(Movie m);
implement: int findDistance(Actor a, Actor b)
6. 白男
System design, open question, give your solution, describe pros and cons


Google
phone： 白男
1. remove duplicates of the array in place
2. 一道BFS题。具体是什么记不清了

on-site:

1. 白男
count islands in a m*n grid （一个联通的值为1的区域被视为一个island）
例：
    0011010
    0010010
    1000110
    0000001

    4 islands found in above grid
Design： copy and shuffle lines in a 8 GB file, memory limit 1 GB (you are
given multiple machines)
2.国人
void minMSwap(int[] num, int m), return the min array after m swaps， each
swap happens only between two adjacent elements([4,2,1,3], 2 return [1,4,2,3
] )
4,2,1,3
4,1,2,3
1,4,2,3
design a protocol to syncing gmail messages among different client apps

3.小白
give a list of <id, parent id, weight>, build the tree(not limited to a
binary tree), then update each node’s sum value(sum is the sum of all its
descendents’ weights)
int[] num incremental(大数加一)
design interface for memory cache

4. 国女
find a intersection to build office so that the sum of all employees’
commute distances is minimum. （the map is represented as a m*n grid, you
are given each employee’s coordination, they can only move in up-down and
left-right directions）

5. 白男manager
How to find median of unsorted integers in linear time
Design the system architecture(FE and BE) for above service in a distributed
system (find optimal office location).


FB
phone: 小白
word break, suffixtree

Onsite:
1. 越南人？
Talked the resume, project and behavior questions. lowest common ancestor
with parent pointer.
2. 白男
is valid binary search tree (handle edge case), if the tree size exist
memory limit, how to handle?
3. 白男
Design question, FB search
4. 国人
give a time, search in a log file. 需要自己提问需求，并考虑边界情况。
    00:23 *****
    00:24 *****
    00:56 *****
    01:02 *****

how to distribute the work to 10 servers?
5. 白男
Celebrity Problem



发信人: seecloud (seecloud), 信区: JobHunting
标  题: 转载的A家onsite 面经
发信站: BBS 未名空间站 (Sat Feb 22 05:56:36 2014, 美东)

想知道这属于什么难度的 是不是不同的组难度有差别？

Round 1 (Written)
1. Given an array, output an array where every index conains nearest
greatest element to that element on right side.

void nearestPeak(int arr[], int len, int out[]) {
	if(len < 2) return;
	int great = arr[len-1];
	out[len-2] = great;
	for(int i=len-3; i>=0; i--) {
		if(arr[i+1]>arr[i+2])
			great=arr[i+1];
		out[i]=great;
	}
}

2. Program to convert sorted array to Binary Search Tree
3. Find first non-repeating character in String
ex: geeksforgeeks: f
geeksforgeeksFirst:o

Round 2 (F2F)
1. Given linked list as a-x-b-y-c-z
output it as a-b-c-z-y-x
that is reverse alternate element and append to end of list

void foo(node* n, node* left, node* right) {
	if(!n) return;
	left->next=n;
	right->next=n->next;
	n=n->next->next;
	left=left->next; right=right->next;
	foo(n, left, right);
}
void do(){
	foo(n, dummy1, dummy2);
	reverse(dummy2);
	append dummy2 to dummy1
}	

2. Output nearest number greater than given number such that output is
palindrome
ex: 121:131
900:909
99:101
vector<int> nearestPal(vector<int> n, int num) {
	int len=n.size();
	vector<int> res;
	IF(len==1) {res.push_back(n[0]+1); return res;}
	if(len%2 == 1) {
		int a=len/2-1; int b=len/2+1;
		if(n[a]>=n[b]) ...


Round 3 (F2F)
1. Vertical Sum in Tree( I told him I know the solution, he proceeded
further)
2. Given stream of Strings find top 5 words with maximum frequency or count
3. Given 2 nodes in Binary Tree find distance between them

Round 4 (F2F with hiring manager)
1. Projects done so far, HR questions
2. Design Linkedin and find till 2nd level connections and path between 2
connection
for ex: if A is friend of B which is friend of C
print between A and C A-B-C
3. Programming language: Java
About synchronisation, serialization, transient and volatile keyword,
Singleton Class

Round 5 (Bar Raiser)
1. Count Inversion in array that is if i a[j]
Told the solution nlogn of divide and conquer. He asked another solution,
then told by inserting in BST and whenever node goes to left side then
adding 1 and number of children on right side . We have to keep track of
count of right subtree in every node

Round 6 (F2F)
1. HR questions (Why leaving company, projects, SWOT)
2. Program to check for mirror tree
3. Data Structure so that push, pop, getmin, getmax O(1) (using 3 stacks)
4. Data Structure so that push, pop, pop min, pop max
Told Solution till O(logn) by using min heap, max heap with pointers to
doubly linked list nodes



// count and evict
#include <iostream>
#include <vector>
using namespace std;
int count(int num, int div) {
	vector<int> v(num, 0);
	for(int i=0; i<num; i++)
		v[i]=i+1;
	int pos=0, cur=1;
	while(v.size()>1) {
		pos++;
		if(pos>=v.size())
			pos=0;
		cur++;
		if(cur%div == 0) {
			v.erase(v.begin()+pos);
			pos--;
			if(pos<0) pos=0;
		}
	}
	return v[0];
}

int main() {
	cout<<count(6,3);
	return 0;
}


// count in an array where all element appear 3 times, except one, which appears 1 time
#include <iostream>
using namespace std;
int count(int A[], int len) {
	int ones=0, twos=0;//, zeros=0xffffffff;
	for(int i=0; i<len; i++) {
		int newOnes=ones^(~twos & A[i]);
		int newTwos=(ones&A[i]) | (twos & ~A[i]);
		//int newZeros=twos&A[i];
		ones=newOnes;
		twos=newTwos;
		
	}
	return ones;
	
}

int main() {
	int A[]={
		2,3,2,4,4,2,2,2,4,2
	};
	cout<<count(A, 10);
}

// build bst from linked list
node* build(node * &n, int len) {
	if(len<=0) {
		return NULL;
	}
	int leftLen=len/2;
	node* leftChild=build(n, leftLen-1);
	node* parent=n;
	parent->left=leftChild;
	n=n->next;
	node* rightChild=build(n, len-leftLen);
	parent->right=rightChild;
	n=n->next;
	return parent;
}
	
// BST to sorted double LL
void convert(node* n, node* &left, node* &right) {
	node* leftSubLeft=NULL;
	node* leftSubRight=NULL;
	node* rightSubLeft=NULL;
	node* rightSubRight=NULL;
	if(n->left) {
		convert(n->left, leftSubLeft, leftSubRight);
		n->left=leftSubRight;
		leftSubRight->right=n;
		left=leftSubLeft;
	} else {
		left=n;
	}
	if(n->right) {
		convert(n->right, rightSubLeft, rightSubRight);
		n->right=rightSubLeft;
		rightSubLeft->left=n;
		right=rightSubRight;
	} else {
		right=n;
	}
}
	
	
#include <iostream>
#include <vector>
using namespace std;
//  compute coin
void genRec(int t, int cand[], int cans, int clen, vector<int>& cur) {
	if(t<0) return;
	if(t==0) {
		for(int i=0; i<cur.size(); i++)
			cout<<cur[i];
		cout<<endl;
		return;
	}
	if(cans < clen) {
		cur.push_back(cand[cans]);
		genRec(t-cand[cans], cand, cans, clen, cur);
		cur.pop_back();
		genRec(t, cand, cans+1, clen, cur);
	}
}
int main() {
	int candidate[]={5,2,1};
	vector<int> current;
	genRec(17, candidate, 0, 3, current);
}

// is graph a tree?
#include <vector>
#include <map>
#include <set>
#include <iostream>
using namespace std;

struct node{
	int val;
	vector<node*> nodes;
};
void BFS(node& n, set<int>& vNums) {
	if(vNums.find(n.val) != vNums.end())
		return;
	vNums.insert(n.val);
	for(vector<node*>::iterator ni = n.nodes.begin();  ni != n.nodes.end(); ++ni) {
		BFS(*(*ni), vNums);
	}
};
	
bool isTree(vector<pair<int,int> >& edges) {
	map<int, node> visited;
	for(vector<pair<int, int> >::iterator it=edges.begin(); it!=edges.end(); ++it) {
		cout<<"iteration"<<endl;
		if(visited.find((*it).first) ==visited.end()) {//not found
			node one;
			one.val=(*it).first;
			visited[(*it).first]=one;
			cout<<"insert "<<one.val<<endl;
		}
		if(visited.find((*it).second) ==visited.end()) {//not found
			node one;
			one.val=(*it).second;
			visited[(*it).second]=one;
			cout<<"insert "<<one.val<<endl;
		}
		visited[(*it).first].nodes.push_back(&visited[(*it).second]);
		visited[(*it).second].nodes.push_back(&visited[(*it).first]);
	}
	if(visited.size() != edges.size()+1) {
		cout<<visited.size()<<" edge size "<<edges.size()<<endl;
		return false;
	} 
	set<int> visitedNums;
	BFS(visited.begin()->second, visitedNums);
	return visitedNums.size() == visited.size();
};

int main() {
	vector<pair<int,int> > edges;
	pair<int, int> e1, e2, e3;
	e1.first=1; e1.second=2;
	e2.first=2; e2.second=3;
	e3.first=5; e3.second=4;
	
	edges.push_back(e1);
	edges.push_back(e2);
	edges.push_back(e3);
	
	cout<<isTree(edges);
	
}
	
	
// inorder
void foo(node* n) {
	if(n) {
		foo(n->left);
		(*func)(n);
		foo(n->right);
	}
}
void foo2(node* root) {
	stack<node*> st;
	if(!root) return;
	//st.push(root);
	node* cur=root;
	while(1) {
		while(cur) {
			st.push(cur);
			cur=cur->left;
		}
		if(st.empty()) break;
		cur=st.top();
		st.pop();
		(*func)(cur);
		cur=cur->right;
	}
}
// preorder
void preOrder(node* root) {
	stack<node*> st;
	if(!root) return;
	st.push(root);
	while(st.empty()==false) {
		node* n=st.top();
		st.pop();
		cout<<n->val<<endl;
		if(n->right) st.push(n->right);
		if(n->left) st.push(n->left);
	}
}
// postOrder
void postOrder(node* root) {
	stack<node*> st;
	if(!root) return;
	node* previous=NULL;
	while(!st.empty()) {
		node* cur=st.top();
		if(cur->right && previous==cur->right || !cur->right && cur->left && previous==cur->left || !cur->right && !cur->left) {
			cout<<cur->val<<endl;
			st.pop();
			previous=cur;
		} else {
			if(cur->right) st.push(cur->right);
			if(cur->left) st.push(cur->left);
		}
	}
}
	
	
	
// Trie implementation	
#include <string>
#include <vector>
#include <iostream>
using namespace std;

#define CHARSIZE 26
class Trie {
public:
struct node {
	vector<node*> child;
	bool isWord;
	node() : isWord(false), child(vector<node*>(CHARSIZE, false)) {};
};
void insert(string s) {
	insertChar(s, &root);
}
void insertChar(string s, node* n) {
	if(s.empty() == true) {
		n->isWord = true;
		return;
	}
	char c = s[0];
	if(!n->child[c-'a'])
		n->child[c-'a'] = new node;
	insertChar(s.substr(1), n->child[c-'a']);
}
bool find(string s) {
	return findChar(s, &root);
}
bool findChar(string s, node* n) {
	if(s.empty())
		return n->isWord;
	char c=s[0];
	if(!n->child[c-'a'])
		return false;
	else
		return findChar(s.substr(1), n->child[c-'a']);
}
void remove(node* n) {
	if(!n) return;
	for(int i=0; i<CHARSIZE; i++) {
		remove(n->child[i]);
	}
	cout<<"removing a node"<<endl;
	delete n;
}
~Trie() {
	remove(&root);
}
private:
	node root;
};

int main() {
	Trie trie;
	trie.insert("abc");
	trie.insert("abcd");
	trie.insert("ks");
	trie.insert("phd");
	trie.insert("wow");
	cout<<trie.find("abc")<<endl;
	cout<<trie.find("ab")<<endl;
	cout<<trie.find("abcd")<<endl;
	cout<<trie.find("abdc")<<endl;
	cout<<trie.find("")<<endl;
	trie.insert("");
	cout<<trie.find("")<<endl;
	return 0;
}

// end Trie	
	
	
//UNION FIND
#include <iostream>
#include <map>
#include <set>
using namespace std;

template <class T>
class UnionFind {
	public:
		bool Uni(T a, T b) {
			if(member.find(a)==member.end() || member.find(b)==member.end())
				return false;
			leader[Find(b)] = Find(a);
			return true;
		}
		T Find(T a) {
			if(leader[a] != leader[leader[a]])
				leader[a] = Find(leader[a]);
			return leader[a];
		}
		bool add(T a) {
			if(member.find(a) == member.end()) {
				member.insert(a);
				leader[a] = a;
				return true;
			}
			std::cout<<"already member"<<endl;
			return false;
		}
	private:
		map<T, T> leader;
		set<T> member;
};
int main() {
	UnionFind<int*> uf;
	int arr[]={
		1,2,3,4,5,6,7
	};
	for(int i=0;i<7; i++)
		uf.add(&arr[i]);
	if(uf.Find(&arr[0]) == uf.Find(&arr[1]))
		cout<<"0 same 1"<<endl;
	else
		cout<<"not same"<<endl;
	uf.Uni(&arr[0], &arr[1]);
	if(uf.Find(&arr[0]) == uf.Find(&arr[1]))
		cout<<"0 same 1"<<endl;
	else
		cout<<"not same"<<endl;
	uf.Uni(&arr[1], &arr[6]);
	uf.Uni(&arr[6], &arr[3]);
	if(uf.Find(&arr[0]) == uf.Find(&arr[3]))
		cout<<"0 same 3"<<endl;
	else
		cout<<"not same 0 3"<<endl;
	if(uf.Find(&arr[0]) == uf.Find(&arr[2]))
		cout<<"0 same 2"<<endl;
	else
		cout<<"not same 0 2"<<endl;
	return 0;
}
// End Union Find



// Calendar Class
struct eventNode {
	string content;
	clock_t start;
	clock_t end;	
};
	
class comp {
	bool operator<(eventNode &e1, eventNode &e2) {
		return e1.start < e2.start;
	}
};

set<eventNode, comp> q;

void notify(eventNode e) {
	cout<<e.content<<endl;
}

bool fin;
// main thread, periodically wake up and check / dispatch calendar event
void processQueue() {
	while (fin) {
	if(q.empty())
		sleep(5);
	else {
		event=*q.begin();
		clock_t curTime = clock();
		if(equal(curTime, event.start)) {
			// handle event
			thread t(notify, event);
		} else if(curTime < event.start) {
			std::this_thread::sleep_for(event.start - curTime);
		}
	}
}

bool insertEvnet(eventNode e) {
	q.insert(e);
	wakeUpProcessQueue(); // wake up main thread
	return true;
}

bool removeEvent(eventNode e) {
	q.erase(e);
	wakeUpProcessQueue();
	return true;
}
void finishCalendar() {
	fin = true;
}

void run() {
	fin = false;
	thread processThread(processQueue);
}

// End Calendar class


// getMinNumber with m move
int getMinNumber(vector<int> arr, int m) {
	set<pair<int,int> > val2Idx;
	for(int i=0; i<arr.size(); i++) {
		pair<int, int> item={arr[i], i};
		set.insert(item);
	}
	int curFillIndex=0;
	while(m) {
		pair<int, int> top = *val2Idx.begin();
		val2Idx.erase(val2Idx.begin());
		if(top.second - curFillIndex <= m) {
			cout<<top.first;
			m = m - (top.second - curFillIndex);
			curFillIndex++;
		} 
	}
	return 0;
}

// Euler cycle
set<Edge*> usedEdge;
doubleLL<Node*> findACycle(node *n, set<Node*>& haveUnusedEdges) {
	doubleLL<Node*> res;
	Edge* edge = getAnUnusedEdge(n);
	if(!edge)
		continue;
	haveUnusedEdges.push(n);
	usedEdge.insert(edge);
	res.push_back(n);
	node *cur=edge.otherSide(n);
	while(cur != n) {
		Edge * e = getAnUnusedEdge(cur);
		haveUnusedEdges.insert(cur);
		usedEdge.insert(e);
		res.push_back(cur);
		cur = e.otherSide(cur);
	}
	return res;
}

void buildDoubleLinkedList(Graph *g) {
	doubleLL<Node*> nodes;
	set<Node*> haveUnusedEdges;
	haveUnusedEdges.insert(g->nodes().begin());
	while(haveUnusedEdges.empty() == false) {
		Node* one = haveUnusedEdges.begin();
		haveUnusedEdges.erase(haveUnusedEdges.begin());
		doubleLL<Node*> onePath = findACycle(one, haveUnusedEdges);
		if(onePath.size() > 0) {//append onePath to existing result, nodes
			doubleLL.appendAfterNode(one, onePath);
		}
	}
	return nodes;
}
		
// naive hash table with chaining for collision handling
template<class K, V>
class Hash {
	public:
	struct node{
		K key;
		V val;
		node* prev;
		node* next;
		node(K k, V v) : key(key), V(v), next(NULL), prev(NULL) {};
	};
	int size;
	node* table;
	Hash(int tableSize) : size(tableSize) {
		table = new node*[size];
		for(int i=0; i<size; i++)
			table[i] = NULL;
	};
	~Hash() {
		for(int i=0; i<size; i++)
			remove(table[i]); // remove linkedlist
		delete table;
	}
	int slot(K key) {
		return hash(K) % size;
	}
	void insert(K key, V val) {
		int s = slot(key);
		if(!table[s]) {
			node* n=new(key, val);
			table[s]=n;
		} else { // collision
			node* cur=table[s];
			while(cur && cur->key != key) {
				cur = cur->next;
			}
			if(!cur) {// no same key in table
				node *n=new(key, val);
				n->next = table[s];
				table[s]->prev = n;
				table[s] = n;
			} else {
				cur->val = val;
			}
		}
	}
	bool get(K key) {};
	V& get(K key) {};
	bool remove(K key) {
		if(get(key) == false)
			return false;
		int s = slot(key);
		node* cur=table[s];
		while(cur && cur->key != key) {
			cur = cur->next;
		}
		// must have found cur == key
		node *prevNode = cur->prev;
		node *nextNode = cur->next;
		if(prevNode)
			prevNode->next = nextNode;
		else
			table[s] = nextNode;
		if(nextNode)
			nextNode->prev = prevNode;
		delete cur;
	}
};
			
			
// Dijkstra
// DijkNode is a symbol for representing node in the original graph
class DijkNode{
	vector<DijkNode*> neighbours;
	int dis; // distance at current stage
	bool operator<(DijkNode& n1) {
		return this->dis > ni.dis;
	}
};
unordered_map<DijkNode*, int> shortest; //decided distance
unordered_map<pair<DijkNode*, DijkNode*>, int> edges;
priority_queue<DijkNode*> pq;
int shortestPath(DijkNode* src, DijkNode* dest, vector<DijkNode*> nodes) {
	// set all distance to infinity
	for(auto ni=nodes.begin(); ni != nodes.end(); ni++)
		*ni->dis = INT_MAX;
	src->dis = 0;
	pq.push(src);
	while(!pq.empty()) {
		DijkNode* top = pq.top(); pq.pop();
		if(shortest.find(top) != shortest.end())
			continue;
		shortest[top] = top->dis;
		if(top == dest)
			break;
		// find all children, update their dis at current stage
		for(auto ni = top->neighbours.begin(); ni != top->neighbours.end(); ni++) {
			if(shortest.find(*ni) != shortest.end()) // this neighbour's distance has been decided
				continue;
			int newDistance = top->dis + edges[make_pair(top, *ni)];
			if(newDistance < *ni->dis) {
				DijkNode* newNode = new DijkNode(*ni); // make of copy of this node
				newNode->dis = newDistance;
				pq.push(newNode);
			}
		}
	}
	return shortest[dest]; 
}


// DFS check graph has cycle
set<node*> visited;
bool DFS(node *n, node* prev) {
	if(visited.find(n))
		return true;
	else
		visited.insert(n);
	for(int i=0; i<n->neighbour.size(); i++) {
		if(n->neighbour[i] != prev) {
			bool res = DFS(n->neighbour[i], n);
			if(res)
				return true;
		}
	}
	return false;
}

bool wrapper(node *n) {
	return DFS(n, NULL);
}


// given a double between 0 - 1, print the binary representation
void foo(double d) {
	for(int i=0; i<10; i++) {
		d=2.0*d;
		if(d>=1) {
			cout<<"1"; d=d-1;
		} else {
			cout<<"0";
		}
	}
	cout<<endl;
}


// merge sort
void mergeSort(int arr[], int low, int high) {
	if(low>=high)
		return;
	int mid = low + ((high - low) >> 1);
	mergeSort(arr, low, mid);
	mergeSort(arr, mid+1, high);
	// merging
	int temp[high-low+1];
	int i=0, lowPtr=low, highPtr=mid+1;
	while(lowPtr <= mid && highPtr <= high) {
		if(arr[lowPtr] <= arr[highPtr]) {
			temp[i++]=arr[lowPtr++];
		} else {
			temp[i++]=arr[highPtr++];
		}
	}
	while(lowPtr <= mid)
		temp[i++]=arr[lowPtr++];
	while(highPtr <= high)
		temp[i++]=arr[highPtr++];
	for(i=0; i<high-low+1; i++)
		arr[low+i] = temp[i];
}

int main() {
	int arr[] = {3,2,1,2,1,5,3,-9,7,0,-1,-2};
	mergeSort(arr, 0, 11);
	for(int i=0; i<12; i++)
		cout<<arr[i]<<" ";
	return 0;
}
			
// quick sort			
void QSort(int arr[], int low, int high) {
	if(low>=high)
		return;
	int pivot=arr[low];
	int s=low+1, e=high;
	while(s <= e) {
		if(arr[s] <= pivot) {
			s++;
		} else {
			swap(arr[s], arr[e--]);
		}
	}
	swap(arr[low], arr[e]);
	QSort(arr, low, e-1);
	QSort(arr, e+1, high);
}
	
// quick select
int QSelect(int arr[], int low, int high, int k) {	
	if(low>=high)
		return arr[low];
	int pivot=arr[low];
	int s=low+1, e=high;
	while(s <= e) {
		if(arr[s] <= pivot) {
			s++;
		} else {
			swap(arr[s], arr[e--]);
		}
	}
	swap(arr[low], arr[e]);
	if(k == e-low+1)
		return arr[e];
	else if(k < e-low+1)
		return QSelect(arr, low, e-1, k);
	else
		return QSelect(arr, e+1, high, k-(e-low+1));
}
	
// avl tree
class node{
	public:
	int val;
	int height;
	node* left;
	node* right;
	node(int v):val(v), left(NULL), right(NULL), height(1) {
		
	};
};
class AVL{
	public:
	node* root;
	AVL() : root(NULL) {};
int getH(node *n) {
	return n ? n->height : 0;
}
node* rotateLeft(node* n) {
	cout<<"rotate left "<<n->val<<endl;
	node* rightChild=n->right;
	n->right=rightChild->left;
	rightChild->left=n;
	n->height=max(getH(n->left), getH(n->right)) + 1;
	rightChild->height=max(getH(rightChild->left), getH(rightChild->right)) + 1;
	return rightChild;
}
node* rotateRight(node* n) {
	cout<<"rotate right "<<n->val<<endl;
	node* leftChild=n->left;
	n->left=leftChild->right;
	leftChild->right=n;
	n->height=max(getH(n->left), getH(n->right)) + 1;
	leftChild->height=1 + max(getH(leftChild->left), getH(leftChild->right));
	return leftChild;
}
int balance(node *n) {
	int t = getH(n->left) - getH(n->right);
	return t;
}
node* insert(node* n, int x) {
	if(!n) {
		cout<<"insert node "<<x<<endl;
		return new node(x);
	}
	if(n->val >= x)
		n->left = insert(n->left, x);
	else
		n->right = insert(n->right, x);
	n->height = 1 + max(getH(n->left), getH(n->right));
	int t = balance(n);
	if(t<=1 && t>=-1)
		return n;
	else if(t < -1 && x > n->right->val) { // right-right
		return rotateLeft(n);
	} else if(t > 1 && x <= n->left->val) {
		return rotateRight(n);
	} else if(t < -1) {
		n->right = rotateRight(n->right);
		return rotateLeft(n);
	} else if( t > 1) {
		n->left = rotateLeft(n->left);
		return rotateRight(n);
	}
}
void insertWrapper(int x) {
	cout<<"root="<<root<<endl;
	root = insert(root, x);
}
};
int main() {
	AVL avl;
	avl.insertWrapper(1);
	avl.insertWrapper(2);
	avl.insertWrapper(1);
	avl.insertWrapper(2);
	avl.insertWrapper(4);
	avl.insertWrapper(-3);
	avl.insertWrapper(6);
	avl.insertWrapper(5);
	avl.insertWrapper(7);
	avl.insertWrapper(8);
	return 0;
}
	

// generate 1 0 string
void gen(string in) {
	vector<string> res;
	stack<string> stk;
	if(in[0]=='1')
		stk.push("1");
	else if(in[0]=='0');
		stk.push("0");
	else {
		stk.push("1"); stk.push("0");
	}
	while(!stk.empty()) {
		string one=stk.top(); stk.pop();
		int idx=one.length();
		if(idx == in.length()) {
			res.push_back(one);
			cout<<one<<endl;
		} else if(in[0]=='1') {
			one+='1';
			stk.push(one);
		} else if(in[0]=='0') {
			one+='0';
			stk.push(one);
		} else {
			stk.push(one+'0'); stk.push(one+'1');
		}
	}
}
	
	
enum cardName {Kill, Miss, Peach, Wine};
enum cardType {Basic, Trick, Suite};
enum cardSuite {Heart, Diamond, Spade, Club};
enum General {guanyu, zhangfei, liubei};
enum Identity {Host, Anti, Evil, Loyal};
class Card{
	int value;	
	cardType type;
	cardName name;
};
class Deck{
	shuffle(List);
	distribute();
	int totalPlayer;
	List<Card> cards;
	List<Card> discarded;
};
class Wear{
	Card* weapon;
	Card* armor;
	Card* offenseHorse;
	Card* defenseHorse;
};

class Player{
	String name;
	General general;
	Identity identity;
	int tableID;
	Wear wear;
	// reset the hasKilled var, etc. 
	Prelude();
	bool hasKilled;
	getCard(Deck);
	// could be multiple cards, with multiple players
	action(List<Card> cardSet, List<Player> playerSet);
	discardCard();
	List<Card> checkRegion;
	List<Card> inHand;
	
	enum Status {flip, passive, prelude, check, ppickCards, act, finish};
	Status status;
	int HPRange;
	int HP;
}


//balanced binary tree
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    int height(TreeNode *node, bool &isFlat, bool &isBalanced){
		if(!node){
			isFlat=true;
			isBalanced=true;
			return 0;
		}
		bool leftFlat, rightFlat, leftBalanced, rightBalanced;
		int leftH=height(node->left, leftFlat, leftBalanced);
		int rightH=height(node->right, rightFlat, rightBalanced);
		if(!leftBalanced || !rightBalanced){
			isBalanced=false;
			return -1;
		}
		//balanced children
		if(leftH==rightH){
			if(leftFlat && rightFlat){
				isFlat=true;
				isBalanced=true;
				return leftH+1;
			}
			else{
				isFlat=false;
				isBalanced=true;
				return leftH+1;
			}
		}
		else if(leftH==rightH+1){

				isFlat=false;
				isBalanced=true;
				return leftH+1;
			
		}
		else if(rightH==leftH+1){
		
				isFlat=false;
				isBalanced=true;
				return rightH+1;
			
		}
		else{
			isBalanced=false;
			return -1;
		}
	}
	bool isBalanced(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        bool balanced, flat;
		int h;
		h=height(root, flat, balanced);
		return balanced;
	
    }
};


//length of last word
class Solution {
public:
    bool isWord(const char *s){
		return (*s != ' ');
		}
    int lengthOfLastWord(const char *s) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        bool inWord=false;
		const char * t=s;
		if(!t)
			return 0;
		int len=0;
		int lastLen=0;
		while(*t != '\0'){
			if(isWord(t)){
				inWord=true;
				len++;
				lastLen=len;
			}
			else{
				inWord=false;
				len=0;
			}
		t++;
		}
		return lastLen;
    }
};



//minimun path sum
class Solution {
public:
    int minPathSum(vector<vector<int> > &grid) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int width=grid.size();
		int height=grid[0].size();
		int *resArr=new int[width];
		resArr[0]=grid[0][0];
		for(int i=1;i<width;i++)
			resArr[i]=resArr[i-1]+grid[i][0];
		for(int j=1;j<height;j++){
			resArr[0]=resArr[0]+grid[0][j];
			for(int k=1;k<width;k++){
				resArr[k]=grid[k][j]+((resArr[k-1]<resArr[k]) ? resArr[k-1]:resArr[k]);
			}
		}
		return resArr[width-1];
    }
};


// stock  1
class Solution {
public:
    int maxProfit(vector<int> &prices) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		int len = prices.size();
		if(len<2)
			return 0;
		int curMax=0;
		int cur, prev=prices[len-1];
		// search backward for a increasing set of prices
		for(int i=len-2;i>=0;i--){
			cur=prices[i];
			if(cur<=prev){
				curMax = curMax>(prev-cur) ? curMax:(prev-cur);
			}else{
				prev=cur;
			}
		}
		return curMax;
    }
};

// stock 2
// greedy
class Solution {
public:
    int maxProfit(vector<int> &prices) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int len=prices.size();
		bool hold=false;
		int profit=0;
		int a;
		for(int i=0;i<len-1;i++){
			if(!hold){
				if(prices[i+1]>prices[i]){
					a=prices[i];
					hold=true;
				}
				else
					continue;
			}
			else{
				if(prices[i+1]>prices[i])
					continue;
				else{
					profit+=prices[i]-a;
					hold=false;
				}
			}
		}
		if(hold)
			profit+=prices[len-1]-a;
		return profit;
    }
};

// stock 3
class Solution {
public:
    int maxProfit(vector<int> &prices) {
		int len = prices.size();
		if(len<2)
			return 0;
		int curMax=0;
		int backMax[len];
		backMax[len-1]=0;
		int cur, prev=prices[len-1];
		// search backward for a increasing set of prices
		for(int i=len-2;i>=0;i--){
			cur=prices[j];
			if(cur<=prev){
				curMax = max(curMax, prev-cur); // curMax>(prev-cur) ? curMax:(prev-cur);
				backMax[i]=curMax;
			}else{
				prev=cur;
				backMax[i]=backMax[i+1];
			}
		}
		// search forward
		int forwardMax[len];
		forwardMax[0]=0;
		int ma=prices[0], mi=prices[0];
		for(int j=1; j<len; j++) {
			cur=prices[i];
			if(cur<mi){
				mi=cur;
				ma=cur;
				forwardMax[j]=forwardMax[j-1];
			} else if(cur>ma) {
				ma=cur;
				int profit=ma-mi;
				forwardMax[j]=max(forwardMax[j-1], profit);
			} else
				forwardMax[j]=forwardMax[j-1];
		}
		int sum=0;
		for(int i=1;i<len-2;i++)
			sum=max(sum, forwardMax[i]+backMax[i+1]);
		sum=max(sum, backMax[0]);
		return sum;
    }
};

// bt -> linked list
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    TreeNode *flat(TreeNode *n){
		TreeNode *temp=n->right;
		TreeNode *original=n;
		if(n->left)
			n->right=flat(n->left);
		else
			n->right=NULL;
		n->left=NULL;
		if(!temp){
			return original;
		} else {
			while(n->right){
				n=n->right;
			}
			n->right=flat(temp);
			return original;
		}
	}
		
    void flatten(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		if(root)
			flat(root);
    }
};

// n queen 
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
			}
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

// partition list
/**
 * Definition for singly-linked list
 * struct ListNode {
 *     int val;
 *     ListNode *next;
 *     ListNode(int x) : val(x), next(NULL) {}
 * };
 */
class Solution {
public:
    ListNode *partition(ListNode *head, int x) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        ListNode* lastSmall=NULL;
		ListNode* firstBig=NULL;
		ListNode* lastBig=NULL;
		if(head==NULL)
			return NULL;
		ListNode* t=head;
		ListNode* firstSmall=NULL;
		while(t){
			ListNode* temp=t->next;
			if(t->val < x){
				if(!lastSmall){
					firstSmall=t;
					lastSmall=t;
					lastSmall->next=firstBig;
				} else {
					lastSmall->next=t;
					t->next=firstBig;
					lastSmall=t;
				}
			} else {
				if(!firstBig){
					firstBig=t;
					lastBig=t;
					firstBig->next=NULL;
					if(lastSmall)
						lastSmall->next=firstBig;
				} else {
					lastBig->next=t;
					lastBig=t;
					lastBig->next=NULL;
				}
			}
			t=temp;
		}
		if(firstSmall)
			return firstSmall;
		else if(firstBig)
			return firstBig;
    }
};



// pascal's triangle II
class Solution {
public:
    vector<int> getRow(int rowIndex) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<int> res;
		res.push_back(1);
		for(int i=0;i<rowIndex;i++){
			res.push_back(1);
			for(int j=i;j>0;j--){
				res[j]=res[j-1]+res[j];
			}
		}
		return res;
    }
};

//triangle 
class Solution {
public:
	int min(int a,int b){
		return (a<b) ? a:b;
	}
    int minimumTotal(vector<vector<int> > &triangle) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<int> > res=triangle;
        int level=triangle.size();
		for(int i=1;i<level;i++){
			int len=triangle[i].size();
			res[i][0]=res[i-1][0]+res[i][0];
			res[i][len-1]=res[i][len-1]+res[i-1][len-2];
			for(int j=1;j<len-1;j++)
				res[i][j]=res[i][j]+min(res[i-1][j-1], res[i-1][j]);
			}
		int minPath=res[level-1][0];
		for(int i=1;i<res[level-1].size();i++)
			minPath=min(minPath,res[level-1][i]);
		return minPath;
    }
};

//is BST
class Solution {
public:
    int getBiggestNodeVal(TreeNode *n){
		while(n->right){
			n=n->right;
		}
		return n->val;
	}
	int getSmallestNodeVal(TreeNode *n){
		while(n->left){
			n=n->left;
		}
		return n->val;
	}
	bool compare(TreeNode *n){
		if(!n->left && !n->right)
			return true;
		if(!n->left && n->right)
			return n->val < getSmallestNodeVal(n->right);
		if(n->left && !n->right)
			return n->val > getBiggestNodeVal(n->left);
		if(n->left && n->right)
			return (n->val < getSmallestNodeVal(n->right)) && (n->val > getBiggestNodeVal(n->left));
	}
	bool isValidBST(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(!root)
			return true;
		if(!compare(root))
			return false;
		return isValidBST(root->left) && isValidBST(root->right);
    }
};

// binary tree maximum path
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
	int max(int a, int b){
		return (a>b) ? a:b;
	}
	int globalMax;
    int maxSinglePath(TreeNode *n){
		// n is not NULL
		if(!n->left && !n->right){
			return n->val;
		}else if(n->left && !n->right){
			return max(n->val, n->val + maxSinglePath(n->left));
		} else if (!n->left && n->right){
			return max(n->val, n->val + maxSinglePath(n->right));
		} else {
			return max(n->val, n->val + max(maxSinglePath(n->left), maxSinglePath(n->right)));
		}
	}
	int maxPathSum(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int thisAsMid=root->val;
		if(!root->left && !root->right){
			return root->val;
		} else if(!root->left && root->right){
			return max(max(root->val, root->val+maxSinglePath(root->right)), maxPathSum(root->right));
		} else if(root->left && !root->right){
			return max(max(root->val, root->val+maxSinglePath(root->left)), maxPathSum(root->left));
		} else{
			int leftS=maxSinglePath(root->left);
			int rightS=maxSinglePath(root->right);
			int maxCur=root->val;
			if(leftS>0)
				maxCur+=leftS;
			if(rightS>0)
				maxCur+=rightS;
			int leftM=maxPathSum(root->left);
			int rightM=maxPathSum(root->right);
			return max(max(leftM,rightM), maxCur);
		}
    }
};

class Solution {
public:
	int max(int a, int b){
		return (a>b) ? a:b;
	}
	int globalMax;
    void updateGlobalMax(int a){
		globalMax = max(globalMax, a);
	}
	int maxSinglePath(TreeNode *n){
		// n is not NULL
		if(!n->left && !n->right){
			updateGlobalMax(n->val);
			return n->val;
		}else if(n->left && !n->right){
			int curMax=max(n->val, n->val + maxSinglePath(n->left));
			updateGlobalMax(curMax);
			return curMax;
		} else if (!n->left && n->right){
			int curMax=max(n->val, n->val + maxSinglePath(n->right));
			updateGlobalMax(curMax);
			return curMax;
		} else {
			int leftMax=maxSinglePath(n->left);
			int rightMax=maxSinglePath(n->right);
			int curMax=max(n->val, n->val + max(leftMax, rightMax));
			int updateVal=n->val+max(leftMax, 0)+max(rightMax, 0);
			updateGlobalMax(updateVal);
			return curMax;
		}
	}
	int maxPathSum(TreeNode *root) {
       globalMax=-99999;
	   maxSinglePath(root);
	   return globalMax;
	  }
};


//  binary tree zigzag traversal
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    vector<vector<int> > zigzagLevelOrder(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<int> > res;
		stack<TreeNode *> even, odd;
		bool oddTurn=true;
		if(!root)
			return NULL;
		odd.push(root);
		while(1){
			if(oddTurn){
				if(odd.empty()){
					break;
				} else {
					vector<int> level;
					while(!odd.empty()){
						TreeNode* n=odd.top();
						odd.pop();

						level.push_back(n->val);
						if(n->left)
							even.push(n->left);
						if(n->right)
							even.push(n->right);
						}
					res.push_back(level);
					oddTurn=false;
				}
			} else {
				if(even.empty())
					break;
				else {
					vector<int> level;
					while(!even.empty()){
						TreeNode* n=even.top();
						even.pop();
						level.push_back(n->val);
						if(n->right)
							odd.push(n->right);
						if(n->left)
							odd.push(n->left);
						
						}
					res.push_back(level);
					oddTurn=true;
				}
			}
		}
		return res;
    }
};

//interleave string
class Solution {
public:
    bool isInterleave(string s1, string s2, string s3) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function    
		if(s1.length()+s2.length() != s3.length())
			return false;
		if(s1.length()==0 && s2.length()==0 && s3.length()==0)
			return true;
		bool match1=false, match2=false;
		if(!s1.empty() && s1[0]==s3[0]){
			string s1new=s1, s3new=s3;
			s1new.erase(0,1);
			s3new.erase(0,1);
			match1=isInterleave(s1new,s2,s3new);
		}
		if(!s2.empty() && s2[0]==s3[0]){
			string s2new=s2, s3new2=s3;
			s2new.erase(0,1);
			s3new2.erase(0,1);
			match2=isInterleave(s1,s2new,s3new2);
		}
		return match1 || match2;
    }
};

//solution 2 dp
class Solution {
public:
	// F(m+1,n+1)=F(m,n+1)&(s3[n+1]==s2[n+1] + F(m+1,n)&(==)
    bool isInterleave(string s1, string s2, string s3) {
		if(s1.length()==0)
			return s2==s3;
		if(s2.length()==0)
			return s1==s3;
		if(s1.length()+s2.length() != s3.length())
			return false;
		bool canInterleave[s1.length()+1][s2.length()+1];
		canInterleave[0][0]=true;
		for(int i=1;i<s1.length()+1;i++)
			canInterleave[i][0]=canInterleave[i-1][0]&&(s1[i-1]==s3[i-1]);
		for(int i=1;i<s2.length()+1;i++)
			canInterleave[0][i]=canInterleave[0][i-1]&&(s2[i-1]==s3[i-1]);
		for(int i=1;i<s1.length()+1;i++){
			int x=i, y=1;
			while(x>0 && y<s2.length()+1){
				canInterleave[x][y]=(canInterleave[x-1][y] && (s1[x-1]==s3[x+y-1])) || (canInterleave[x][y-1] && (s2[y-1]==s3[x+y-1]));
				x--;
				y++;
			}
		}
		for(int i=2;i<s2.length()+1;i++){
			int y=i, x=s1.length();
			while(x>0 && y<s2.length()+1){
				canInterleave[x][y]=(canInterleave[x-1][y] && (s1[x-1]==s3[x+y-1])) || (canInterleave[x][y-1] && (s2[y-1]==s3[x+y-1]));
				x--;
				y++;
			}
		}
		return canInterleave[s1.length()][s2.length()];
    }
};

//longest common prefix
class Solution {
public:
    string longestCommonPrefix(vector<string> &strs) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int len=strs.size();

        string res;
        if(len==1)
            return strs[0];
        if(len==0)
            return res;
        char c;
        int t=0, shortestStr=strs[0].size();
        for(int i=1;i<len;i++){
            if(strs[i].size()<shortestStr)
                shortestStr=strs[i].size();
        }
        bool end=false;
        while(t<shortestStr){
        for(int i=0;i<len;i++){
            if(i==0)
                c=strs[i][t];
            else if(strs[i][t] != c){
                end=true;
                break;
            }
        }
        if(end)
            break;
        else{
            res.append(1,c);
            t++;
        }
        }
        return res;
    }
};

// populate next right pointer and II
class Solution {
public:
    void connect(TreeLinkNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        queue<TreeLinkNode *> q;
		if(!root)
			return;
		q.push(root);
		q.push(NULL);
		TreeLinkNode* prev=NULL;
		while(!q.empty()){
			TreeLinkNode* node=q.front();
			q.pop();
			if(!node){
				if(q.empty())
					break;
				q.push(NULL);
				prev=NULL;
				continue;
			}
			if(node->left)
				q.push(node->left);
			if(node->right)
				q.push(node->right);
			if(prev)
				prev->next=node;
			node->next=NULL;
			prev=node;
		}
    }
};
// constant space solution 
class Solution {
public:
    void connect(TreeLinkNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		if(!root)
			return;
		TreeLinkNode* begin=NULL, prev=NULL, cur=root;
		bool processLeft=false; // processed cur's left child?
		while(1){
			if(!begin){
				while(cur){
					if(cur->left){
						begin=cur->left;
						prev=begin;
						processLeft=true;
						break;
					}else if(cur->right){
						begin=cur->right;
						cur=cur->next;
						prev=begin;
						processLeft=false;
						break;
					}
					cur=cur->next;
				}
				if(!begin)
					//end of loop
					break;
			}else{
				while(cur){
					if(processLeft){//process right only
						if(cur->right){
							prev->next=cur->right;
							prev=cur->right;
						}
						processLeft=false;
						cur=cur->next;
					} else {//process left if any
						if(cur->left){
							prev->next=cur->left;
							prev=cur->left;
							processLeft=true;
						}else{//no left
							processLeft=true;
						}
					}
				}
				// end of a level
				cur=begin;
				begin=NULL; prev=NULL;
			}
		}
	}
};

//  exchange wrong elements in bst
/* Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
	TreeNode* arr[2];
	TreeNode* wrongNode[2];
	void check(TreeNode* n){
		if(arr[1] != NULL){
			arr[0]=arr[1];
			arr[1]=n;
			//arr[2]=n;
			if(arr[0]->val > arr[1]->val)
				if(wrongNode[0])
					wrongNode[1]=arr[1];
				else{
					wrongNode[0]=arr[0];
					wrongNode[1]=arr[1];
				}
		}else
			arr[1]=n;
	}	
	void traverse(TreeNode* n){
		if(!n)
			return;
		traverse(n->left);
		check(n);
		traverse(n->right);
	}
    void recoverTree(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        arr[0]=NULL;
		arr[1]=NULL;
		wrongNode[0]=NULL;
		wrongNode[1]=NULL;
		//arr[2]=NULL;
		traverse(root);
		if(!wrongNode[0] || !wrongNode[1])
			return;
		int temp=wrongNode[0]->val;
		wrongNode[0]->val=wrongNode[1]->val;
		wrongNode[1]->val=temp;
    }
};

// restore ip address
class Solution {
public:
	vector<string> result;
	bool isValid(string s, int a, int b){
		int len=s.length();
		if(a<0 || a>=len)
			return false;
		if(b<0 || b>=len)
			return false;
		if(a>b)
			return false;
		if(s[a]=='0' && a!=b)
			return false;
		int num=0;
		for(int i=a;i<=b;i++)
			num=num*10+(s[i]-'0');
		if(num>=0 && num<256)
			return true;
		else
			return false;
	}
    void restore(string s, int seg, string addr){
		if(s.length()<seg || s.length()>seg*3)
			return;
		if(seg==1){
			if(isValid(s,0,s.length()-1)){
				addr.append(".");
				addr.append(s);
				if(addr[0]=='.')
					addr.erase(0,1);
				result.push_back(addr);
			}
			}else{
			for(int len=1;len<4;len++){
				if(isValid(s,0,len-1)){
					string newAddr=addr;
					newAddr.append(".");
					newAddr.append(s,0,len);
					string newS=s;
					newS.erase(0,len);
					restore(newS, seg-1, newAddr);
				}
			}
		}
		}
	vector<string> restoreIpAddresses(string s) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		string addr;
		result.clear();
		restore(s, 4, addr);
		return result;
    }
};


//  two sum
class Solution {
public:
    vector<int> twoSum(vector<int> &numbers, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int len=numbers.size();
		vector<int> oldNumbers=numbers;
		vector<int> res;
		for(int i=0;i<len;i++){
			int min=numbers[i];
			for(int j=i;j<len;j++)
				if(numbers[j]<min){
					numbers[i]=numbers[j];
					numbers[j]=min;
					min=numbers[i];
				}
		}
		int a=0,aa,bb;
		int b=len-1;
		while(a<b){
			if(numbers[a]+numbers[b]<target){
				a++;
			}
			else if(numbers[a]+numbers[b]>target)
				b--;
			else
				break;
		}
		for(int i=0;i<len;i++){
			if(oldNumbers[i]==numbers[a])
				aa=i;
			if(oldNumbers[i]==numbers[b])
				bb=i;
			}
		int t;
		if(aa>bb){
			t=aa;
			aa=bb;
			bb=t;
		}
		aa++;bb++;
		res.push_back(aa);
		res.push_back(bb);
		return res;
	}
};

//unique binary tree I
class Solution {
public:
    int numTrees(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int res=0;
        if(n==0)
            return 1;
        if(n==1)
            return 1;
        if(n==2)
            return 2;
        for(int i=0;i<n;i++){
            res+=numTrees(i)*numTrees(n-i-1);
        }
        return res;
    }
};

//unique binary tree II
class Solution {
public:
	vector<TreeNode* > generate(int a, int b){
		//start at a, end at b (non inclusive)
		vector<TreeNode* > res;
		if(a==b-1){
			TreeNode* node=new TreeNode(a);
			res.push_back(node);
			return res;
		} else if(a>=b){
			// empty node;
			TreeNode* node=NULL;
			res.push_back(node);
			return res;
		} else {
			for(int i=a;i<b;i++){
				vector<TreeNode* > leftSubs=generate(a,i);
				vector<TreeNode* > rightSubs=generate(i+1, b);
				for(int p=0;p<leftSubs.size();p++){
					for(int q=0;q<rightSubs.size();q++){
						TreeNode* node=new TreeNode(i);
						node->left=leftSubs[p];
						node->right=rightSubs[q];
						res.push_back(node);
					}
				}
			}
			return res;
		}
	}
    vector<TreeNode *> generateTrees(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        return generate(1, n+1);
    }
};

// sudoku solver
class Solution {
public:
    void clearReg(bool reg[9]){
		for(int i=0;i<9;i++)
			reg[i]=false;
		}
	int isNum(char c){
		if(c>='1' && c<='9')
			return (c-'0');
		else
			return 0;
	}
	bool registered[9];
	bool isValid(vector<vector<char> > &board, int x, int y){
		clearReg(registered);
		for(int i=0;i<9;i++){
			if(!isNum(board[x][i]))
				continue;
			else if(registered[isNum(board[x][i])-1] == false)
				registered[isNum(board[x][i])-1] = true;
			else
				return false;
		}
		clearReg(registered);
		for(int i=0;i<9;i++){
			if(!isNum(board[i][y]))
				continue;
			else if(registered[isNum(board[i][y])-1] == false)
				registered[isNum(board[i][y])-1] = true;
			else
				return false;
		}
		clearReg(registered);
		int xS=x-x%3;
		int yS=y-y%3;
		for(int i=xS;i<xS+3;i++)
			for(int j=yS;j<yS+3;j++){
				if(!isNum(board[i][j]))
					continue;
				else if(registered[isNum(board[i][j])-1] == false)
					registered[isNum(board[i][j])-1] = true;
				else
					return false;
			}
		return true;
	}
	vector<vector<char> > tempBoard;
	void storeBoard(vector<vector<char> > &board){
		tempBoard=board;
	}
	void solve(vector<vector<char> > &board, int x, int y){
    	bool finalCube=false;

        if(x==9 && y==0)
			finalCube=true;
		if(finalCube){
			storeBoard(board);
			return;
		}
		int xNext,yNext;
		if(y==8){
			xNext=x+1; yNext=0;
		}else{
			xNext=x; yNext=y+1;
		}
		
		if(board[x][y] != '.')
			solve(board, xNext, yNext);
		else{
			for(int i=1;i<=9;i++){
				board[x][y]='0'+i;
				if(isValid(board, x, y))
					solve(board, xNext, yNext);
			}
			// revert the current cube
			board[x][y]='.';
		}
	}
    void solveSudoku(vector<vector<char> > &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        solve(board, 0, 0);
		for(int i=0;i<9;i++)
			for(int j=0;j<9;j++)
				board[i][j]=tempBoard[i][j];
    }
};

// valid sudokuclass Solution 
class Solution {
public:
    void clearReg(bool reg[9]){
    	for(int i=0;i<9;i++)
			reg[i]=false;
		}
	int isNum(char c){
		if(c>='1' && c<='9')
			return (c-'0');
		else
			return 0;
	}
	bool registered[9];
	bool validSudoku;
	bool isValid(vector<vector<char> > &board, int x, int y){
		clearReg(registered);
		for(int i=0;i<9;i++){
			if(!isNum(board[x][i]))
				continue;
			else if(registered[isNum(board[x][i])-1] == false)
				registered[isNum(board[x][i])-1] = true;
			else
				return false;
		}
		clearReg(registered);
		for(int i=0;i<9;i++){
			if(!isNum(board[i][y]))
				continue;
			else if(registered[isNum(board[i][y])-1] == false)
				registered[isNum(board[i][y])-1] = true;
			else
				return false;
		}
		clearReg(registered);
		int xS=x-x%3;
		int yS=y-y%3;
		for(int i=xS;i<xS+3;i++)
			for(int j=yS;j<yS+3;j++){
				if(!isNum(board[i][j]))
					continue;
				else if(registered[isNum(board[i][j])-1] == false)
					registered[isNum(board[i][j])-1] = true;
				else
					return false;
			}
		return true;
	}
	vector<vector<char> > tempBoard;
	void storeBoard(vector<vector<char> > &board){
		tempBoard=board;
		validSudoku=true;
	}
	void solve(vector<vector<char> > &board, int x, int y){
		if(validSudoku)
			return;
    	bool finalCube=false;

        if(x==9 && y==0)
			finalCube=true;
		if(finalCube){
			storeBoard(board);
			return;
		}
		int xNext,yNext;
		if(y==8){
			xNext=x+1; yNext=0;
		}else{
			xNext=x; yNext=y+1;
		}
		
		if(board[x][y] != '.')
			solve(board, xNext, yNext);
		else{
			for(int i=1;i<=9;i++){
				board[x][y]='0'+i;
				if(isValid(board, x, y))
					solve(board, xNext, yNext);
			}
			// revert the current cube
			board[x][y]='.';
		}
	}
    bool isValidSudoku(vector<vector<char> > &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        validSudoku=false;
		solve(board, 0, 0);
		if(validSudoku)
			return true;
		else
			return false;
		for(int i=0;i<9;i++)
			for(int j=0;j<9;j++)
				board[i][j]=tempBoard[i][j];
    }
};

// distinct subsequence
// DP solution, the countSubSeq() routine is not used.
class Solution {
public:
	int count;
	void countSubSeq(string& s, int si, string& t, int ti) {
		if(si>=s.length() || ti>=t.length())
			return;
		if(t.length()==ti+1){
		// matching last char
			for(int i=si;i<s.length();i++)
				if(s[i]==t[ti])
					count++;
		}else{
			for(int i=si;i<s.length();i++)
				if(s[i]==t[ti])
					countSubSeq(s, i+1, t, ti+1);
		}
	}
	// subCount[t][s]=sigma(subCount[t+1][ss]) where S[ss-1]=T[t]
    int numDistinct(string S, string T) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        count=0;
		int tchar=T.length();
		int schar=S.length();
		if(tchar==0 || schar==0)
			return 0;
		int subCount[tchar][schar];
		int initCount=0;
		for(int i=schar-1;i>=0;i--){
			if(S[i]==T[tchar-1])
				initCount++;
			subCount[tchar-1][i]=initCount;
		}
		for(int i=0;i<tchar-1;i++)
			subCount[i][schar-1]=0;
		for(int i=tchar-2;i>=0;i--)
			for(int j=schar-2;j>=0;j--){
				if(T[i]!=S[j])
					subCount[i][j]=subCount[i][j+1];
				else{
					subCount[i][j]=subCount[i][j+1]+subCount[i+1][j+1];
				}
			}
		return subCount[0][0];
		//countSubSeq(S, 0, T, 0);
		//return count;
    }
};

// maximal Rectangle
class Solution {
public:
	int maxSpace;
	int checkOnes(int i, int j, int c, vector<vector<char> > &matrix){
		int row=matrix.size();
		int col=matrix[0].size();
		if(c>=col)
			return 0;
		int cOrig=c;
		int cDest=c;
		int width=j-1+1;
		bool hasZero=false;
		while(c<col && !hasZero){
			for(int x=i;x<=j;x++)
				if(matrix[x][c]=='0'){
					cDest=c-1;
					hasZero=true;
					break;
				}
			c++;
		}
		if(!hasZero)
			cDest=col-1;
		return width*(cDest-cOrig+1);
	}
    int maximalRectangle(vector<vector<char> > &matrix) {
		maxSpace=0;
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int row=matrix.size();
		int col=matrix[0].size();
		for(int c=0;c<col;c++){
			for(int i=0;i<row;i++){
				if(matrix[i][c]=='0')
					continue;
				for(int j=i;j<row;j++){
					if(matrix[j][c]=='0')
						break;
					int space=checkOnes(i,j,c,matrix);
					updateSpace(space);
				}
			}
		}
					
    }
};

// valid palindrome
class Solution {
public:
    bool nonChar(char c){
            if(c>='0' && c<='9')
                return false;
            if(c>='a' && c<='z')
                return false;
            if(c>='A' && c<='Z')
                return false;
            return true;
    }
    bool isP(char a, char b){
            if(a==b)
                return true;
            if(a==b-'A'+'a')
                return true;
            if(a==b-'a'+'A')
                return true;
            return false;
    }
    bool isPalindrome(string s) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(s.length()==0)
            return true;
        int a, b=s.length()-1;
        while(a<b){
            if(nonChar(s[a]))
                a++;
            else if(nonChar(s[b]))
                b--;
            else if(!isP(s[a], s[b])
                break;
            else{
                a++;
                b--;
            }
        }
        if(a>=b)
            return true;
        else
            return false;
    }
};



//surrounded regions

class Solution {
public:
    int len, height;
    void recoverOrColor(vector<vector<char>> &b, int x,int y, bool color) {
            if(x<0 ||x> len-1 ||y<0 || y>height-1 )
                return;
            else if(b[x][y]=='X')
                return;
            else if(b[x][y]=='R') {
                if(color)
                    b[x][y]=='X';
                else
                    b[x][y]=='O';
                recoverOrColor(b, x-1,y,color);
                recoverOrColor(b, x+1,y,color);
                recoverOrColor(b, x,y-1,color);
                recoverOrColor(b, x,y+1,color);
            }
			else
				return;
    }
    bool isRegion(vector<vector<char>> &b, int x, int y) {
        if(b[x][y]=='X')
            return true;
        else if(x==0 || x==len-1 || y==0 || y==height-1)
            return false;
        
        else if(b[x][y]=='O') {
                b[x][y]='R';
                bool left=(b[x-1][y]=='R') ? true:isRegion(b, x-1,y);
                bool right=(b[x+1][y]=='R') ? true:isRegion(b, x+1,y);
                bool up=(b[x][y-1]=='R') ? true:isRegion(b, x,y-1);
                bool down=(b[x][y+1]=='R') ? true:isRegion(b, x,y+1);
                if(left && right && up && down)
                    return true;
                else
                    return false;
        }
		else
			return false;
    }
    void solve(vector<vector<char>> &board) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        len=board.size();//x
        if(len==0)
            return;
        height=board[0].size();//y
        if(len ==0 || height== 0)
            return;
        for(int i=0;i<len;i++)
            for(int j=0;j<height;j++) {
                if(board[i][j]=='O') {
                    if(isRegion(board, i, j))
                        recoverOrColor(board, i, j, true);
                    else
                        recoverOrColor(board, i, j, false);
                }
            }
    }
};

//longest consecutive sequence

class Solution {
public:
   set<int> m;
    int findAndDelete(int x) {
        if(m.find(x)==m.end())
            return 0;
        //found, now extend left, right
        int t=x-1;
        int ret=1;
        m.erase(x);
        while(m.find(t) != m.end()) {
            ret++;
            m.erase(t--);
        }
        t=x+1;
        while(m.find(t) != m.end()) {
            ret++;
            m.erase(t++);
        }
        return ret;
    }
    int longestConsecutive(vector<int> &num) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        for(int i=0; i<num.size(); ++i)
            m.insert(num[i]);
        int maxlen=0;
        for(int i=0; i<num.size(); ++i) {
            int len=findAndDelete(num[i]);
            if(len>maxlen)
                maxlen=len;
        }
        return maxlen;
    }
};



//Palindrome Partitioning
class Solution {
public:
    set< pair<int,int> > pa;
    vector<vector<string> > res;
	void find(vector<string> curV, int startPos, string &s) {
		if(startPos >= s.length()){
			res.push_back(curV);
			return;
		}
		for(int ep=startPos; ep<s.length(); ep++) {
			if(pa.find(pair<int,int>(startPos, ep)) != pa.end()) {
				vector<string> newV=curV;
				newV.push_back(s.substr(startPos, ep-startPos+1));
				find(newV, ep+1, s);
			}
		}
	}
	vector<vector<string> > partition(string s) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		res.clear();
		pa.clear();
        for(int i=0; i<s.length(); i++)
			pa.insert(pair<int,int>(i,i));
		for(int i=0; i<s.length()-1; i++)
			if(s[i]==s[i+1])
				pa.insert(pair<int,int>(i,i+1));
		for(int i=3;i<=s.length();i++) {
			for(int start=0; start<s.length()-(i-1); start++) {
				if(s[start]==s[start+i-1] && pa.find(pair<int,int>(start+1, start+i-2)) != pa.end())
					pa.insert(pair<int,int>(start, start+i-1));
				}
			}
		vector<string> vec;
		find(vec,0,s);
		return res;
    }
};



//Palindrome Partitioning II
class Solution {
public:
    set< pair<int,int> > pa;
    vector<vector<string> > res;
    int minCuts;
	void update(int cuts) {
		if(cuts<minCuts)
			minCuts=cuts;
	}
	void find(vector<string> curV, int startPos, string &s, int cuts) {
		if(startPos >= s.length()){
			res.push_back(curV);
			update(cuts);
			return;
		}
		for(int ep=startPos; ep<s.length(); ep++) {
			if(pa.find(pair<int,int>(startPos, ep)) != pa.end()) {
				vector<string> newV=curV;
				newV.push_back(s.substr(startPos, ep-startPos+1));
				find(newV, ep+1, s, cuts+1);
			}
		}
	}
	int minCut(string s) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
		res.clear();
		pa.clear();
		minCuts=s.length();
        for(int i=0; i<s.length(); i++)
			pa.insert(pair<int,int>(i,i));
		for(int i=0; i<s.length()-1; i++)
			if(s[i]==s[i+1])
				pa.insert(pair<int,int>(i,i+1));
		for(int i=3;i<=s.length();i++) {
			for(int start=0; start<s.length()-(i-1); start++) {
				if(s[start]==s[start+i-1] && pa.find(pair<int,int>(start+1, start+i-2)) != pa.end())
					pa.insert(pair<int,int>(start, start+i-1));
				}
			}
		map<pair<int,int>, int> minCutMap;
		for(int i=0; i<s.length(); i++)
			minCutMap[pair<int,int>(i,i)]=0;
		for(int i=0; i<s.length()-1; i++)
			if(pa.find(pair<int,int>(i,i+1)) != pa.end())
				minCutMap[pair<int,int>(i,i+1)]=0;
			else
				minCutMap[pair<int,int>(i,i+1)]=1;
		for(int i=3; i<=s.length(); i++) {
			//for each interval length, compute mincuts
			for(int start=0; start<s.length()-(i-1); start++) {
				if(pa.find(pair<int,int>(start, start+i-1)) != pa.end()) {
					minCutMap[pair<int,int>(start, start+i-1)] = 0;
				} else {// not a palindrome
					int min = minCutMap[pair<int,int>(start,start)] + minCutMap[pair<int,int>(start+1,start+i-1)]; //+1 later
					for(int c=start+2; c<=start+i-1; c++) {
						int newCut = minCutMap[pair<int,int>(start,c-1)] + minCutMap[pair<int,int>(c,start+i-1)];
						if(newCut<min)
							min=newCut;
					}
					min++;
					minCutMap[pair<int,int>(start,start+i-1)] = min;
				}
			}
		}
		return minCutMap[pair<int,int>(0,s.length()-1)];
	
    }
};



// word ladder II
class Solution {
public:
void rec(unordered_map<string, unordered_set<string> > &list, vector<vector<string> > &res, int wid, int step, string end, vector<string> 	vec, unordered_set<string> &dict){
	if(step==0){
     if(vec.back()==end)
          res.push_back(vec);
     return;
	}
	string curStr=vec.back();
	for(unordered_set<string>::iterator it = list[curStr].begin(); it != list[curStr].end(); it++) {
		vec.push_back(*it);
		rec(list, res, wid, step-1, end, vec, dict);
		vec.pop_back();
     }
}
int findLadderss(unordered_map<string, unordered_set<string> > &list, string start, string end, unordered_set<string> &dict) {
// Start typing your C/C++ solution below
// DO NOT write int main() function
	int step=1;
	queue<string> q;
	unordered_set<string> visited;	
	int level=1;
	int wid=start.length();
	q.push(start);
	bool found=false;
	q.push(string(""));
	while(!q.empty()){
	string s=q.front();
	q.pop();
	if(s == "") {
	if(q.empty())
	break;
	q.push(string(""));
	level++;
	continue;
	}
	visited.insert(s);
	for(unordered_set<string>::iterator it = list[s].begin(); it != list[s].end(); it++) {
     if(*it == end){
          found=true;
          return level;
     } else if(visited.find(*it) == visited.end())
          q.push(*t);
	}
	}
	if(!found)
	return -1;
}

vector<vector<string>> findLadders(string start, string end, unordered_set<string> &dict) {
//build the adjacent list
unordered_map<string, unordered_set<string> > list;
int wid=end.length();
for(unordered_set<string>::iterator it=dict.begin(); it != dict.end(); it++) {
     unordered_set<string> thisList;
     for(int i=0; i<wid; i++){
          for(int j=0; j<26; j++){
               char c='a'+j;
               string t=(*it);
               t[i]=c;
               if(dict.find(t) != dict.end() && t != (*it))
                    thisList.insert(t);
               }
          }
     list[*it] = thisList;
     }
//built the list
int level=findLadderss(list, start, end,dict); vector<vector<string> > res;
if(level==-1)
return res;
vector<string> vec;
vec.push_back(start);

rec(list, res, end.length(), level, end, vec, dict);
return res;
}
};

//////////////////////
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <queue>
using namespace std;

class Solution {
public:
bool rec(map<string, int> &bo, map<string, set<string> > &list, vector<vector<string> > &res, int wid, int step, string end, vector<string> 	vec, set<string> &dict){
	if(step==0){
     if(vec.back()==end){
          res.push_back(vec);
     	cout<<"find end"<<endl;
     	return true;
     }
     else {
     	return false;
     }
	}
	for(vector<string>::iterator vi=vec.begin(); vi != vec.end(); vi++)
		cout<<*vi<<"-";
	cout<<endl;
	string curStr=vec.back();
	if(bo[curStr] >= step)
		return false;
	bool found=false;
	for(set<string>::iterator it = list[curStr].begin(); it != list[curStr].end(); it++) {
		bool dup=false;
		for(vector<string>::iterator vi=vec.begin(); vi != vec.end(); vi++){
			if(*vi == *it)
				dup=true;
		}
		if(!dup){
		vec.push_back(*it);
		bool ret = rec(bo, list, res, wid, step-1, end, vec, dict);
		vec.pop_back();	
		found = found || ret;
		}
     }
     if(!found){
     	int t=bo[curStr];
     	if(step>t)
     		bo[curStr]=step;
     }
     return found;
}
void backSearch(vector<vector<string> > &res, vector<string> vec, int level, map<string, int> &dist, map<string, set<string> > &list, set<string> &dict) {
	//go from end to start, following the dist
	if(level==0){
		res.push_back(vec);
		return;
	}
	string curStr=vec.front();
	for(set<string>::iterator it = list[curStr].begin(); it != list[curStr].end(); it++) {
		if(dist.find(*it)!=dist.end() && dist[*it]==level) {
			vector<string> newVec=vec;
			newVec.insert(newVec.begin(), *t);
			backSearch(res, newVec, level-1, dist, list, dict);
		}
	}
}
	
	
	
int findLadderss(map<string, int> &dist, map<string, set<string> > &list, string start, string end, set<string> &dict) {
// Start typing your C/C++ solution below
// DO NOT write int main() function
	int step=1;
	queue<string> q;
	set<string> visited;	
	int level=1;
	int wid=start.length();
	q.push(start);
	bool found=false;
	q.push(string(""));
	while(!q.empty()){
	string s=q.front();
	q.pop();
	if(s == "") {
	if(q.empty())
	break;
	q.push(string(""));
	level++;
	continue;
	}
	visited.insert(s);
	dist[s]=level;
	for(set<string>::iterator it = list[s].begin(); it != list[s].end(); it++) {
     if(*it == end){
          found=true;
          return level;
     } else if(visited.find(*it) == visited.end())
          q.push(*it);
	}
	}
	if(!found)
	return -1;
}

vector<vector<string> > findLadders(string start, string end, set<string> &dict) {
//build the adjacent list
map<string, set<string> > list;
int wid=end.length();
for(set<string>::iterator it=dict.begin(); it != dict.end(); it++) {
     set<string> thisList;
     for(int i=0; i<wid; i++){
          for(int j=0; j<26; j++){
               char c='a'+j;
               string t=(*it);
               t[i]=c;
               if(dict.find(t) != dict.end() && t != (*it))
                    thisList.insert(t);
               }
          }
     list[*it] = thisList;
     }
//built the list
cout<<"built the list"<<endl;
map<string, int> dist;
int level=findLadderss(list, start, end,dict); 
vector<vector<string> > res;
cout<<"level"<<level<<endl;
if(level==-1)
return res;
vector<string> vec;
vec.push_back(end);
//map<string, int> bo;
//rec(bo, list, res, end.length(), level, end, vec, dict);
backSearch(res, vec, level, dist, list, dict);
cout<<"total ans found"<<res.size()<<endl;
return res;
}
};

int main(){
	string start="cet", end="ism";
	set<string> dict;
	char* words[]={"kid","tag","pup","ail","tun","woo","erg","luz","brr","gay","sip","kay","per","val","mes","ohs","now","boa","cet","pal","bar","die","war","hay","eco","pub","lob","rue","fry","lit","rex","jan","cot","bid","ali","pay","col","gum","ger","row","won","dan","rum","fad","tut","sag","yip","sui","ark","has","zip","fez","own","ump","dis","ads","max","jaw","out","btu","ana","gap","cry","led","abe","box","ore","pig","fie","toy","fat","cal","lie","noh","sew","ono","tam","flu","mgm","ply","awe","pry","tit","tie","yet","too","tax","jim","san","pan","map","ski","ova","wed","non","wac","nut","why","bye","lye","oct","old","fin","feb","chi","sap","owl","log","tod","dot","bow","fob","for","joe","ivy","fan","age","fax","hip","jib","mel","hus","sob","ifs","tab","ara","dab","jag","jar","arm","lot","tom","sax","tex","yum","pei","wen","wry","ire","irk","far","mew","wit","doe","gas","rte","ian","pot","ask","wag","hag","amy","nag","ron","soy","gin","don","tug","fay","vic","boo","nam","ave","buy","sop","but","orb","fen","paw","his","sub","bob","yea","oft","inn","rod","yam","pew","web","hod","hun","gyp","wei","wis","rob","gad","pie","mon","dog","bib","rub","ere","dig","era","cat","fox","bee","mod","day","apr","vie","nev","jam","pam","new","aye","ani","and","ibm","yap","can","pyx","tar","kin","fog","hum","pip","cup","dye","lyx","jog","nun","par","wan","fey","bus","oak","bad","ats","set","qom","vat","eat","pus","rev","axe","ion","six","ila","lao","mom","mas","pro","few","opt","poe","art","ash","oar","cap","lop","may","shy","rid","bat","sum","rim","fee","bmw","sky","maj","hue","thy","ava","rap","den","fla","auk","cox","ibo","hey","saw","vim","sec","ltd","you","its","tat","dew","eva","tog","ram","let","see","zit","maw","nix","ate","gig","rep","owe","ind","hog","eve","sam","zoo","any","dow","cod","bed","vet","ham","sis","hex","via","fir","nod","mao","aug","mum","hoe","bah","hal","keg","hew","zed","tow","gog","ass","dem","who","bet","gos","son","ear","spy","kit","boy","due","sen","oaf","mix","hep","fur","ada","bin","nil","mia","ewe","hit","fix","sad","rib","eye","hop","haw","wax","mid","tad","ken","wad","rye","pap","bog","gut","ito","woe","our","ado","sin","mad","ray","hon","roy","dip","hen","iva","lug","asp","hui","yak","bay","poi","yep","bun","try","lad","elm","nat","wyo","gym","dug","toe","dee","wig","sly","rip","geo","cog","pas","zen","odd","nan","lay","pod","fit","hem","joy","bum","rio","yon","dec","leg","put","sue","dim","pet","yaw","nub","bit","bur","sid","sun","oil","red","doc","moe","caw","eel","dix","cub","end","gem","off","yew","hug","pop","tub","sgt","lid","pun","ton","sol","din","yup","jab","pea","bug","gag","mil","jig","hub","low","did","tin","get","gte","sox","lei","mig","fig","lon","use","ban","flo","nov","jut","bag","mir","sty","lap","two","ins","con","ant","net","tux","ode","stu","mug","cad","nap","gun","fop","tot","sow","sal","sic","ted","wot","del","imp","cob","way","ann","tan","mci","job","wet","ism","err","him","all","pad","hah","hie","aim","ike","jed","ego","mac","baa","min","com","ill","was","cab","ago","ina","big","ilk","gal","tap","duh","ola","ran","lab","top","gob","hot","ora","tia","kip","han","met","hut","she","sac","fed","goo","tee","ell","not","act","gil","rut","ala","ape","rig","cid","god","duo","lin","aid","gel","awl","lag","elf","liz","ref","aha","fib","oho","tho","her","nor","ace","adz","fun","ned","coo","win","tao","coy","van","man","pit","guy","foe","hid","mai","sup","jay","hob","mow","jot","are","pol","arc","lax","aft","alb","len","air","pug","pox","vow","got","meg","zoe","amp","ale","bud","gee","pin","dun","pat","ten","mob"};
	int i=0;
	while(words[i]){
		cout<<words[i]<<endl;
		dict.insert(string(words[i]));
		i++;
	}
	cout<<dict.size()<<*(dict.begin())<<endl;
	Solution sol;
	sol.findLadders(start, end, dict);
}



// binary tree maximun path sum. passed both small and large.

class Solution {
public:
	map<TreeNode*, int> maxS;
	int getMaxS(TreeNode* n) {
		if(maxS.find(n) != maxS.end())
			return maxS[n];
		// compute
		int leftS, rightS;
		if(n->left)
			leftS=getMaxS(n->left);
		if(n->right)
			rightS=getMaxS(n->right);
		int ret;
		if(n->left && n->right)
			ret = n->val + max(0, max(leftS, rightS));
		else if(n->left && !n->right)
			ret = n->val + max(0,leftS);
		else if(n->right && !n->left)
			ret = n->val + max(0,rightS);
		else
			ret = n->val;
		maxS[n] = ret;
		return ret;
	}
	int maxPath;
	void combine(TreeNode* n) {
		if(!n)
			return;
		if(n->left && n->right)
			maxPath=max(maxPath, n->val + max(0,maxS[n->left])+ max(0,maxS[n->right]));
		else if(n->left && !n->right)
			maxPath=max(maxPath, n->val + max(0,maxS[n->left]));
		else if(!n->left && n->right)
			maxPath=max(maxPath, n->val + max(0,maxS[n->right]));
		else
			maxPath=max(maxPath, n->val);
		combine(n->left);
		combine(n->right);
	}
		
    int maxPathSum(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        maxPath=-999;
		maxS.clear();
		getMaxS(root);
		combine(root);
		return maxPath;
    }
};

//path sum
class Solution {
public:
	int s;
	bool found;
	void func(TreeNode* n, int sum, int curSum) {
		curSum+=n->val;
		if(!n->left && !n->right){
			if(curSum==sum)
				found=true;
			
		}
		if(n->left)
			func(n->left, sum, curSum);
		if(n->right)
			func(n->right, sum, curSum);
	}
    bool hasPathSum(TreeNode *root, int sum) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        found=false;
		s=0;
		if(root==NULL)
			return false;
		func(root, sum, 0);
		return found;
    }
};
//path sum II
class Solution {
public:
	vector<vector<int> > res;
	void func(TreeNode* n, int sum, int curSum, vector<int> curV) {
		curSum+=n->val;
		curV.push_back(n->val);
		if(!n->left && !n->right){
			if(curSum==sum)
				res.push_back(curV);
		}
		if(n->left)
			func(n->left, sum, curSum, curV);
		if(n->right)
			func(n->right, sum, curSum, curV);
	}
    vector<vector<int> > pathSum(TreeNode *root, int sum) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        res.clear();
		if(root==NULL)
			return res;
		vector<int> vec;
		func(root, sum, 0, vec);
		return res;
    }
};

//minimal depth
class Solution {
public:
    int minDepth(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(!root)
			return 0;
		if(!root->left && !root->right)
			return 1;
		if(!root->left && root->right)
			return minDepth(root->right)+1;
		if(!root->right && root->left)
			return 1+minDepth(root->left);
			
		return 1+min(minDepth(root->left), minDepth(root->right));
    }
};

// convert sorted list to binary search tree
/**
 * Definition for singly-linked list.
 * struct ListNode {
 *     int val;
 *     ListNode *next;
 *     ListNode(int x) : val(x), next(NULL) {}
 * };
 */
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    // 1 or more ... 
    ListNode *findMid(ListNode * n){
        int len=0;
        ListNode *nn=n;
        while(nn){
            nn=nn->next;
            len++;
        }
        for(int i=0;i<len/2-1;   i++)
            n=n->next;
        //ListNode *mid=n->next;
        //n->next=NULL;
        return n;
    }
    TreeNode* sortRec(ListNode *n){
        
        if(!n->next){
            TreeNode* single=new TreeNode(n->val);
            return single;
        }
        ListNode *mid=findMid(n);
        ListNode *realMid=mid->next;
        TreeNode *tn=new TreeNode(realMid->val);
        mid->next=NULL;
        ListNode *rightStart=realMid->next;
        realMid->next=NULL;
        tn->left=sortRec(n);
        if(rightStart)
            tn->right=sortRec(rightStart);
        return tn;
    }
    TreeNode *sortedListToBST(ListNode *head) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(!head)
            return NULL;
        return sortRec(head);
    }
};

//convert sorted array to binary search tree

class Solution {
public:
    TreeNode* convert(vector<int> &v, int a, int b){
        if(a>b)
            return NULL;
        if(a==b){
            TreeNode *tn=new TreeNode(v[a]);
            return tn;
        }
        int mid=(b-a)/2+a;
        TreeNode *tnn=new TreeNode(v[mid]);
        tnn->left=convert(v, a, mid-1);
        tnn->right=convert(v, mid+1, b);
        return tnn;
    }
    TreeNode *sortedArrayToBST(vector<int> &num) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        return convert(num, 0, num.size()-1);
    }
};

//binary tree level order traversal
class Solution {
public:
           vector<vector<int> > res;

    vector<vector<int> > levelOrderBottom(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
    res.clear();
    if(!root)
            return res;
        stack<vector<int> > stk;
        queue<TreeNode *> q;
        q.push(root);
        q.push(NULL);
        vector<int> v;
        //v.push_back(root->val);
        while(!q.empty()){
            TreeNode *n=q.front();
            q.pop();
            if(n==NULL){
                stk.push(v);
                v.clear();
                if(!q.empty()){
                    q.push(NULL);
                    continue;
                }else
                    break;
            }
            v.push_back(n->val);
            if(n->left)
                q.push(n->left);
            if(n->right)
                q.push(n->right);
        }
        while(stk.empty() == false){
            res.push_back(stk.top());
            stk.pop();
        }
        return res;
    }
};

//construct binary tree from inorder and post order traversal
/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    TreeNode *recBuild(vector<int> &in, int ina, int inb, vector<int> &post, int posta, int postb){
        TreeNode *ret=NULL;
        if(ina==inb){
            ret=new TreeNode(in[ina]);
            return ret;
        }
        int i=0;
        int midVal=post[postb];
        ret=new TreeNode(midVal);
        while(in[ina+i] != midVal){
            i++;
        }
        if(i>0)
            ret->left=recBuild(in, ina, ina+i-1, post, posta, posta+i-1);
        if(ina+i<inb)
            ret->right=recBuild(in, ina+i+1, inb, post, posta+i, postb-1);
        return ret;
    }
    TreeNode *buildTree(vector<int> &inorder, vector<int> &postorder) {
        if(inorder.size()==0)
            return NULL;
        return recBuild(inorder, 0, inorder.size()-1, postorder, 0, postorder.size()-1);
    }
};

// maximal depth of binary tree
class Solution {
public:
    int maxDepth(TreeNode *root) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(!root)
            return 0;
        return max(maxDepth(root->left), maxDepth(root->right)) + 1;
    }
};






// Reverse Linked List II
class Solution {
public:
    ListNode *reverseBetween(ListNode *head, int m, int n) {
        ListNode *mNode;
        ListNode *nNode;
        ListNode *b4m = NULL;
        ListNode *aftern = NULL;
        bool mIsLast=false;
        if(m==n)
            return head;
        ListNode* cur=head;
        ListNode* prev;
        for(int i=1;i<=n+1;i++){
            if(i==n+1){
                aftern=cur;
                if(!cur){
                    mIsLast=true;
                    break;
                }
            }
            ListNode* next=cur->next;
            if(i==m-1)
                b4m=cur;
            if(i==m)
                mNode=cur;
            else if(i==n){
                nNode=cur;
                nNode->next=prev;
            }
            else if(i>m && i<n){
                cur->next=prev;
            }
            prev=cur;
            cur=next;
        }
		if(!b4m)
			head=nNode;
		else
			b4m->next=nNode;
		mNode->next=aftern;
		return head;
    }
};

// merge sorted array
class Solution {
public:
    void merge(int A[], int m, int B[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int ai=m-1;
        int bi=n-1;
        int mergei=m+n-1;
        while(ai>=0 || bi>=0) {
            if(ai<0){
                A[mergei--]=B[bi--];
            } else if(bi<0) {
                A[mergei--]=A[ai--];
            }
            else if(A[ai]>=B[bi]) {
                A[mergei--]=A[ai--];
            } else {
                A[mergei--]=B[bi--];
            }
        }
    }
};


// remove dup from sorted list I
class Solution {
public:
    ListNode *deleteDuplicates(ListNode *head) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        ListNode *t=head;
        ListNode *prev=NULL;
        //bool dup=false;
        while(t) {
            if(prev && prev->val == t->val) {
                prev->next=t->next;
                t=prev->next;
            } else {
                prev=t;
                t=t->next;
            }
        }
        return head;
    }
};

// remove dup from sorted list II
class Solution {
public:
    void remove(ListNode* n, int val, ListNode* &ret) {
        if(!n) {
            while(ret && ret->val == val)
                ret=ret->next;
        } else {
            while(n->next && n->next->val == val) {
                n->next=n->next->next;
            }
        }
    }
    ListNode *deleteDuplicates(ListNode *head) {
        ListNode *cur=head;
        ListNode *prev=NULL;
        ListNode *ret=head;
        set<int> dup;
        set<int> visit;
        while(cur) {
            if(visit.find(cur->val) != visit.end())
                dup.insert(cur->val);
            else
                visit.insert(cur->val);
            cur=cur->next;
        }
        cur=head;
        while(cur) {
            if(dup.find(cur->val) != dup.end())
                remove(prev, cur->val, ret);
            else {
                prev=cur;
            }
            int value=cur->val;
            while(cur && cur->val == value)
                cur=cur->next;
        }
        return ret;
    }
};


//word search
class Solution {
public:
    bool existRec(vector<vector<char> > &visited, vector<vector<char> > &board, int x, int y, string word, int i){
        if(i==word.length())
            return true;
        if((x<0 || x>=visited.size()) || (y<0 || y>=visited[0].size()))
            //out of range
            return false;
        if(visited[x][y]=='v')
            return false;
        if(board[x][y] != word[i])
            return false;
        else{
            visited[x][y]='v';
            if(existRec(visited, board, x-1, y, word, i+1))
                return true;
            if(existRec(visited, board, x+1, y, word, i+1))
                return true;
            if(existRec(visited, board, x, y-1, word, i+1))
                return true;
            if(existRec(visited, board, x, y+1, word, i+1))
                return true;
            //unmark
            visited[x][y]='u';
            return false;
        }
        
    }
    bool exist(vector<vector<char> > &board, string word) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(board.size()<1 || board[0].size()<1)
            return false;
        vector<vector<char> > visited = board;
        for(int x=0;x<board.size();x++)
            for(int y=0; y<board[0].size(); y++)
                visited[x][y]='u';
        for(int x=0;x<board.size();x++)
            for(int y=0; y<board[0].size(); y++) {
                
                if(existRec(visited, board, x,y,word, 0))
                    return true;
            }
        return false;
    }
};

//remove dup from sorted array II
class Solution {
public:
    int removeDuplicates(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int min=A[0]-1;
        int cur=0;
        int s=0;
        int curV=A[0];
        int count=0;
        while(cur<n) {
            if(A[cur]==curV)
                count++;
            else {
                if(count>2)
                    for(int t=s+2; t<cur; t++)
                        A[t]=min;
                count=1;
                curV=A[cur];
                s=cur;
            }
            cur++;
        }
        if(count>2)
            for(int t=s+2; t<n; t++)
                    A[t]=min;
        cur=0; s=0;
        while(cur<n) {
            if(A[cur] != min){
                A[s]=A[cur];
                cur++;
                s++;
            
            }
            else
                cur++;
        }
        return s;
    }
};

//find median in sorted arrays
class Solution {
public:
  double getMed(int A[], int a, int b){
    if((b-a)%2==0)
       return A[(a+b)/2];
     else
       return (A[(a+b)/2]+A[(a+b)/2+1])/2.0;
  }
  double getMed2(int arr[], int a, int b, int val) {
	if(a==b)
		return (arr[a]+val)/2.0;
     if((b-a)%2==0) {
       int x=arr[(a+b)/2-1];
       int y=arr[(a+b)/2];
       int z=arr[(a+b)/2+1];
       if(val <= x)
          return (x+y)/(2.0);
       else if(val >= z)
          return (y+z)/(2.0);
       else
          return (val+y)/(2.0);
     } else {
       int x=arr[(a+b)/2];
       int y=arr[(a+b)/2+1];
       if(val <= x)
          return x;
       else if(val >= y)
          return y;
       else
          return val;
     }
  }
  double find(int A[], int a1,int a2, int B[], int b1,int b2) {
	if(a2<a1) {
		if((b2-b1)%2 == 0)
			return B[(b2+b1)/2];
		else
			return (B[(b2+b1)/2]+B[(b2+b1)/2+1])/2.0;
	}
	if(b2<b1) {
		if((a2-a1)%2 == 0)
			return A[(a2+a1)/2];
		else
			return (A[(a2+a1)/2]+A[(a2+a1)/2+1])/2.0;
	}
     if(a1==a2 && b1==b2)
       return (A[a1]+B[b1])/(2.0);
     if(a2==a1) {
       int aval=A[a1];
       if(aval <= B[b1])
          return getMed(B, b1, b2-1);
       else if(aval >= B[b2])
          return getMed(B, b1+1, b2);
       else {//aval within B
          return getMed2(B, b1, b2, aval);
       }
     } else if(b1==b2) {
       int bval=B[b1];
       if(bval <= A[a1])
          return getMed(A, a1, a2-1);
       else if(bval >= A[a2])
          return getMed(A, a1+1, a2);
       else {//aval within B
          return getMed2(A, a1, a2, bval);
       }
     }
     if(a1+1==a2) {
		if(A[a1]<B[b1]) {
			if(A[a2]>B[b2])
				return getMed(B, b1, b2);
			else
				return getMed2(B,b1,b2-1,A[a2]);
		} else {
			if(A[a2]>B[b2])
				return getMed2(B,b1+1,b2,A[a1]);
			else
				return find(A,a1,a2,B,b1+1,b2-1);
			}
		} else if(b1+1==b2) {
		if(B[b1]<A[a1]) {
			if(B[b2]>A[a2])
				return getMed(A, a1, a2);
			else
				return getMed2(A,a1,a2-1,B[b2]);
		} else {
			if(B[b2]>A[a2])
				return getMed2(A,a1+1,a2,B[b1]);
			else
				return find(A,a1+1,a2-1,B,b1,b2);
			}
		}
     double vam, vbm;
	 if((a2-a1)%2 == 0)
		vam = A[(a1+a2)/2];
	else
		vam = (A[(a1+a2)/2]+A[(a1+a2)/2+1])/2.0;
	 if((b2-b1)%2 == 0)
		vbm = B[(b1+b2)/2];
	else
		vbm = (B[(b1+b2)/2]+B[(b1+b2)/2+1])/2.0;
		
	
		
     int am=(a1+a2)/2;
     int bm=(b1+b2)/2;
     int cuta=a2-am;
     if((a2-a1)%2 != 0)
        cuta=cuta-1;
     int cutb=b2-bm;
     if((b2-b1)%2 != 0)
        cutb=cutb-1;
     if(vam>vbm) {
       if(cuta>cutb)
          return find(A, a1, a2-cutb, B, b1+cutb, b2);
       else
          return find(A, a1, a2-cuta, B, b1+cuta, b2);
     } else if(vam<vbm) {
       if(cuta>cutb)
          return find(A, a1+cutb, a2, B, b1, b2-cutb);
       else
          return find(A, a1+cuta, a2, B, b1, b2-cuta);
     } else { //
        if((a2-a1)%2 == 0 || (b2-b1)%2 == 0) {
            return vam;
        } 
        else {
                //both even length
                int x=A[am];
                int y=B[bm];
                int p=A[am+1];
                int q=B[bm+1];
                int z=(p>q) ? q:p;
                int zz=(x>y)?x:y;
                return (zz+z)/2.0;
        }
     }
  }
    double findMedianSortedArrays(int A[], int m, int B[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
       if(m==0)
          return getMed(B, 0, n-1);
       else if(n==0)
          return getMed(A, 0, m-1);
       return find(A, 0, m-1, B, 0, n-1);
    }
};


//longest substring without repeating characters
class Solution {
public:
    int lengthOfLongestSubstring(string s) {
        if(s.length()==0)
            return 0;
        int charIdx[26];
        for(int i=0;i<26;i++)
            charIdx[i]=-1;
        int cur=0, from=0;
        int maxL=-1;
        int nonDupStart=0;
        while(cur<s.length()){
            if(charIdx[s[cur]-'a'] == -1) {
                charIdx[s[cur]-'a'] = cur;
                maxL=max(maxL, cur-nonDupStart+1);
            }
            else {//dup
                int lastAppear=charIdx[s[cur]-'a'];
                if(lastAppear>=nonDupStart) {
                maxL=max(maxL, cur-lastAppear);
                nonDupStart=lastAppear+1;
                } else {
                    maxL=max(maxL, cur-nonDupStart+1);
                }
                charIdx[s[cur]-'a'] = cur;
                
            }
            cur++;
        }
        return maxL;
    }
};

// jump game
class Solution {
public:
    bool canJump(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(n<1)
            return false;
        int reach=0;
        for(int i=0;i<n-1;i++) {
            if(reach>=i)
            reach=max(reach, i+A[i]);
        }
        if(reach>=n-1)
            return true;
        else
            return false;
        //return reachable[n-1];
    }
};

//trap water 
class Solution {
public:
    int trap(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int leftMax[n];
        int rightMax[n];
        
        leftMax[0]=0;
        for(int i=1; i<n; i++) {
            leftMax[i]=max(leftMax[i-1], A[i-1]);
        }
        rightMax[n-1]=0;
        for(int i=n-2; i>=0; i--) {
            rightMax[i]=max(rightMax[i+1], A[i+1]);
        }
        int trapV=0;
        for(int i=1;i<n-1;i++) {
            int bound = min(leftMax[i], rightMax[i]);
            trapV+= (bound>A[i]) ? bound-A[i]:0;
        }
        return trapV;
    }
};

// first missing positive
class Solution {
public:
    int firstMissingPositive(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int cur=0;
        while(cur<n) {
            if(A[cur]==cur+1)
                cur++;
            else if(A[cur]<=0)
                cur++;
            else {
                if(A[cur]>n){
                    A[cur]=-1;
                    cur++;
                }
                else {
                    int temp=A[A[cur]-1];
                    if(temp == A[cur])
                        A[cur]=-1;
                    else{
                    A[A[cur]-1] = A[cur];
                    A[cur] = temp;
                    }
                }
            }
        }
        int i=0;
        while(A[i]==i+1)
            i++;
        return i+1;
    }
};

// search insert position
class Solution {
public:
    int searchInsert(int A[], int n, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(n==1) {
            if(A[0]>=target)
                return 0;
            else return 1;
        }
        if(n==2) {
            if(A[0]>=target)
                return 0;
            else if(A[1]>=target)
                return 1;
            else
                return 2;
        }
        int mid=n/2;
        if(A[mid]>=target)
            return searchInsert(A, mid, target);
        else
            return mid+searchInsert(A+mid, n-mid, target);
    }
};

// LRU implementation  (double linkedlist + hashmap)
class node {
	T val;
	node* prev;
	node* next;
	}
class dLL {
	node* head;
	node* tail;
	int size;
	node* remove(node* n);
	void append(node* n);
	void prepend(node* n);
	}
class LRU {
	int cacheSize;
	dLL dll;
	hashmap<val, node*> hash;
	get(T t) {
		if(node *n=hash.find(t)) {
			dll.remove(n);
			dll.prepend(n);
		} else if(dLL.size == cacheSize) {
			hash.remove(dll.tail->val);
			dll.remove(dll.tail);
			node *newNode;
			newNode->val=T;
			dll.prepend(newNode);
			hash[t]=dll.head;
		} else {
			node *newNode;
			newNode->val=T;
			dll.prepend(newNode);
			hash[t]=dll.head;
			dll.size=dll.size+1;
		}
	}
}

//字符串匹配
#include <set>
// 预处理初始化
#define NUM 29
set<long long> hash;
void initWithString(char *str) {
    int len=strlen(str);
    for(int i=1;i<=10;i++) {
        long long res[len];
        for(int start=0;start+i<=len;start++) {
            // start ... start+i-1
            long long sum=0;
            if(start == 0){
                for(int each=0;each<i;each++)
                    sum=sum*NUM+str[start + each];
                hash.insert(sum);
            } else {
                long long subtract=1;
                for(int k=0; k < i-1; k++)
                    subtract=subtract*NUM;
                sum = (res[start-1] - subtract*str[start-1])*NUM + str[start+i-1];
                hash.insert(sum);
            }
            res[start] = sum;
        }
    }
}
// 如果query是str的字串,返回true,否则返回false
bool existSubString(char *query) {
    int len=strlen(query);
    long long sum=0;
    for(int i=0; i<len; i++) {
        sum=sum*NUM+query[i];
    }
    if(hash.find(sum) != hash.end())
        return true;
    else
        return false;
    
}


// implement strstr
class Solution {
public:
    char *strStr(char *haystack, char *needle) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        char *h=haystack;
        if(!h || !needle)
            return NULL;
        if(*needle == '\0')
            return haystack;
        while(*h !='\0') {
            char *n=needle;
            char *hh=h;
            while(*n !='\0') {
                if(*n != *hh)
                    break;
                n++;
                hh++;
            }
            if(*n == '\0')
                return h;
            else
                h++;
        }
        return NULL;
    }
};

// faster version strstr
class Solution {
public:
    char *strStr(char *haystack, char *needle) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        char *h=haystack;
        if(!h || !needle)
            return NULL;
        if(*needle == '\0')
            return haystack;
        char *ahead=h;
        char *temp=needle;
        while(*temp++)
            ahead++;
        ahead--;
        while(*ahead !='\0') {
            char *n=needle;
            char *hh=h;
            while(*n !='\0') {
                if(*n != *hh)
                    break;
                n++;
                hh++;
            }
            if(*n == '\0')
                return h;
            else{
                h++;
                ahead++;
            }
        }
        return NULL;
    }
};

// remove elements
class Solution {
public:
    int removeElement(int A[], int n, int elem) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int old=0, newI=0;
        while(old<n) {
            if(A[old]==elem)
                old++;
            else
                A[newI++]=A[old++];
        }
        return newI;
    }
};

// remove duplicates from sorted array
class Solution {
public:
    int removeDuplicates(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int old=0, newI=0;
        if(n<1)
            return 0;
        int cur=A[0]-1;
        while(old<n) {
            if(A[old]==cur)
                old++;
            else {
                cur=A[old];
                A[newI++]=A[old++];
            }
        }
        return newI;
    }
};

// multiply string
class Solution {
public:
    string multiply(string num1, string num2) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int len1=num1.length();
        int len2=num2.length();
        string res;
        int r[len1+len2];
        for(int i=0; i<len1+len2; i++)
            r[i]=0;
        for(int i=len2-1; i>=0; i--)
            for(int j=len1-1; j>=0; j--) {
                r[i+j+1]=r[i+j+1]+(num1[j]-'0')*(num2[i]-'0');
            }
        for( int i=len1+len2-1; i>0; i--) {
            if(r[i]>9) {
                int c=r[i]/10;
                r[i]=r[i]%10;
                r[i-1]=r[i-1]+c;
            }
        }
        int idx=0;
        while(r[idx]==0)
            idx++;
        if(idx==len1+len2)
            res.append(1,'0');
        while(idx<len1+len2){
            res.append(1, r[idx]+'0');
            idx++;
        }
        return res;
    }
};

// compute pow
class Solution {
public:
    map<long long, double> res;
    double power(double x, long long n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        if(n==0)
            return 1;
        if(n==1)
            return x;
        if(n<0)
            return 1.0/power(x, -n);
        long long m=n/2;
        double re1=(res.find(m) != res.end()) ? res[m]:power(x, m);
        if(res.find(m)==res.end())
            res[m]=re1;
        double re2=(res.find(n-m) != res.end()) ? res[n-m]:power(x, n-m);
        if(res.find(n-m)==res.end())
            res[n-m]=re2;
        return re1 * re2;
            
    }
    double pow(double x, int n) {
        res.clear();
        return power(x, n);
    }
};

// regular expression
// (use dynamic programming
class Solution {
public:
    bool isMatch(const char *s, const char *p) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function    
        int sc=strlen(s);
        int pl=strlen(p);
        int star=0;
        for(int i=0; i<pl; i++)
            if(p[i]=='*')
                star++;
        int pc=pl-star;
        bool match[sc+1][pl+1];
        match[0][0]=true;
        for(int i=1; i<=sc; i++)
            match[i][0]=false;
        int i;
        for(i=1; i<pl; i++) {
            if(p[i]=='*'){
                match[0][i]=true;
                match[0][i+1]=true;
                i++;
            } else {
                break;
            }
        }
        while(i<=pl) {
            match[0][i]=false;
            i++;
        }
        for(int x=1; x<=sc; x++)
            for(int y=1; y<=pl; y++) {
                if((y<pl && p[y] != '*') || y==pl)
                    match[x][y]=(p[y-1]=='.' || p[y-1]==s[x-1]) ? match[x-1][y-1] : false;
                else if(p[y-1]=='.'){
                    match[x][y]=match[x][y-1] || match[x-1][y];
                    match[x][y+1]=match[x][y-1] || match[x-1][y];
                    y++;
                } else { // a* case
                    bool res=match[x][y-1] || (s[x-1]==p[y-1] && match[x-1][y]);
                    match[x][y]=res;
                    match[x][y+1]=res;
                    y++;
                }
            }
        return match[sc][pl];
    }
};

// 4 sum
class Solution {
public:
    void find2Sum(vector<int> &num, int t, int s, int e, vector<vector<int> > &res, int a, int b) {
        int prevs=b-1;
        int preve=b-1;
        while(s<e) {
            int sum=num[s] + num[e];
            if(sum>t)
                e--;
            else if(sum<t)
                s++;
            else {
                if(num[s]!=prevs && num[e]!=preve) {
                vector<int> add;
                add.push_back(a);
                add.push_back(b);
                add.push_back(num[s]);
                add.push_back(num[e]);
                res.push_back(add);
                
                prevs=num[s];
                preve=num[e];
                }
                s++;
                e--;
            }
        }
    }
    vector<vector<int> > fourSum(vector<int> &num, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<int> > res;
        int size=num.size();
       sort(num.begin(), num.end());
        int a,b,c,d;
        for(a=0; a<size-3; a++){
            if(a>0 && num[a] == num[a-1])
                continue;
            for(b=a+1; b<size-2; b++) {
                if(b>a+1 && num[b] == num [b-1])
                    continue;
                int t=target- num[a] -num[b];
                find2Sum(num, t, b+1, size-1, res, num[a], num[b]);
            }
        }
        return res;
    }
};



// simplify path
class Solution {
public:
    string simplifyPath(string path) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        stack<string> res;
    	int i=0;
		while(i<path.length()) {
			while(i<path.length() && path[i]=='/')
				i++;
			int j=i;
			while(j<path.length() && path[j]!='/')
				j++;
			j--;
			if(j<i) {
				i=j+1;
				continue;
			} else if(j==i && path[i]=='.') {
				i=i+1; continue;
			} else if(j==i+1 && path[i]=='.' && path[j]=='.'){
				if(res.empty()==false)
                    res.pop(); 
                i=j+1; continue;
			} else {
				string p; p.append(1, '/');
				p.append(path.substr(i, j-i+1));
				res.push(p);
                i=j+1;
			}
		}
		string ret;
		while(res.empty()==false) {
			string s=res.top(); res.pop();
			ret.insert(0, s);
		}
        if(ret.empty())
            ret.append(1,'/');
		return ret;
    }
};

// climb stairs
class Solution {
public:
    int climbStairs(int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int ways[n+1];
        ways[1]=1; ways[2]=2;
        for(int i=3;i<=n;i++) {
            ways[i]=ways[i-1]+ways[i-2];
        }
        return ways[n];
    }
};

// max subarray
class Solution {
public:
    int maxSubArray(int A[], int n) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int maxS;
        int prev=A[0];
        maxS=prev;
        for(int i=1; i<n; i++) {
            if(prev<0){
                maxS=max(maxS, A[i]);
                prev=A[i];
            }else{
                maxS=max(maxS, prev+A[i]);
                prev=prev+A[i];
            }
        }
        return maxS;
    }
};

// insert interval
class Solution {
public:
    vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
        int sz=intervals.size();
        vector<Interval> res=intervals;
        vector<Interval> ret;
        for(int i=0; i<sz; i++) {
            if(newInterval.start > res[i].end || newInterval.end < res[i].start){
                ret.push_back(res[i]);
                continue;
            }
            else {
                newInterval.start = min(newInterval.start, res[i].start);
                newInterval.end = max(newInterval.end, res[i].end);
            }
        }
        int retSz=ret.size();
        if(retSz==0){
            ret.push_back(newInterval);
            return ret;
        }
        if(newInterval.start<ret[0].start)
            ret.insert(ret.begin(), newInterval);
        else if(newInterval.end>ret[retSz-1].end)
            ret.insert(ret.end(), newInterval);
        else {
            vector<Interval>::iterator it=ret.begin();
            while((*(it+1)).start < newInterval.start)
                it++;
            ret.insert(it+1, newInterval);
        }
        return ret;
    }
};

// merge intervals
class Solution {
public:
    vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
        int sz=intervals.size();
        vector<Interval> res=intervals;
        vector<Interval> ret;
        for(int i=0; i<sz; i++) {
            if(newInterval.start > res[i].end || newInterval.end < res[i].start){
                ret.push_back(res[i]);
                continue;
            }
            else {
                newInterval.start = min(newInterval.start, res[i].start);
                newInterval.end = max(newInterval.end, res[i].end);
            }
        }
        int retSz=ret.size();
        if(retSz==0){
            ret.push_back(newInterval);
            return ret;
        }
        if(newInterval.start<ret[0].start)
            ret.insert(ret.begin(), newInterval);
        else if(newInterval.end>ret[retSz-1].end)
            ret.insert(ret.end(), newInterval);
        else {
            vector<Interval>::iterator it=ret.begin();
            while((*(it+1)).start < newInterval.start)
                it++;
            ret.insert(it+1, newInterval);
        }
        return ret;
    }
    vector<Interval> merge(vector<Interval> &intervals) {
        vector<Interval> ret;
        int sz=intervals.size();
        if(sz==0)
            return ret;
        ret.push_back(intervals[0]);
        
        for(int i=1; i<sz; i++) {
            ret=insert(ret, intervals[i]);
        }
        return ret;
    }
};

// construct tree from preorder and inorder traversal
class Solution {
public:
    TreeNode *build(vector<int> &pre, int prea, int preb, vector<int> &in, int ina, int inb) {
	if(prea>preb)
		return NULL;
	if(prea==preb) {
		TreeNode *t=new TreeNode(pre[prea]);
        return t;
	}
	int rootV=pre[prea];
	TreeNode *root=new TreeNode(rootV);
	int ini=ina;
	while(in[ini] != rootV)
		ini++;
	root->left=build(pre, prea+1, ini-ina+prea, in, ina, ini-1);
	root->right=build(pre, ini-ina+prea+1, preb, in, ini+1, inb);
	return root;
	}
    TreeNode *buildTree(vector<int> &preorder, vector<int> &inorder) {
		int len=preorder.size();
		return build(preorder, 0, len-1, inorder, 0, len-1);
    }
};

// search in rotated sorted array
class Solution {
public:
    int sea(int A[], int a, int b, int t) {
        if(a>b)
            return -1;
        if(a==b && A[a]==t)
            return a;
        if(a==b && A[a] != t)
            return -1;
        if(A[a]<A[b]){
            int mid=A[(a+b)/2];
            if(mid>=t)
                return sea(A,a,(a+b)/2, t);
            else
                return sea(A,(a+b)/2+1, b, t);
        } else {
            int mid=A[(a+b)/2];
            int midi=(a+b)/2;
            if(t<=A[b]) {
                if(mid>=t && mid<=A[b])
                    return sea(A, a, midi, t);
                else
                    return sea(A, midi+1, b, t);
            } else {
                if(mid>=t || mid<=A[b])
                    return sea(A,a,midi, t);
                else
                    return sea(A,midi+1,b,t);
            }
        }
    }
    int search(int A[], int n, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        return sea(A, 0, n-1, target);
    }
};

// 3 sum closest
class Solution {
public:
    int close;
    void update(int u, int target) {
        int prevGap=target-close;
        if(prevGap<0)
            prevGap=-prevGap;
        int curGap=target-u;
        if(curGap<0)
            curGap=-curGap;
        if(curGap<prevGap)
            close=u;
    }
    int threeSumClosest(vector<int> &num, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
            close=num[0]+num[1]+num[2];
            sort(num.begin(), num.end());
            for(int i=0; i<num.size()-2; i++) {
                int a=num[i];
                int j=i+1, k=num.size()-1;
                int nt=target-a;
                while(j<k) {
                    update(num[j]+num[k]+a, target);
                    if(num[j]+num[k]>nt)
                        k--;
                    else if(num[j]+num[k]<nt)
                        j++;
                    else
                        return target;
                }
            }
            return close;
    }
};


// longest common prefix
class Solution {
public:
    string longestCommonPrefix(vector<string> &strs) {
        string res;
        int sz=strs.size();
        if(sz<1)
            return res;
        if(sz==1)
            return strs[0];
        bool end=false;
        for(int x=0; x<strs[0].length(); x++) {
            char c=strs[0][x];
            for(int y=1; y<sz; y++) {
                if(strs[y].length()<=x){
                    end=true;
                    break;
                }
                if(strs[y][x] != c) {
                    end=true;
                    break;
                }
            }
            if(end)
                break;
            else
                res.push_back(c);
        }
        return res;
    }
};

// wild card matching
class Solution {
public:
    bool isMatch(const char *s, const char *p) {
        int sl=strlen(s);
        int pl=strlen(p);
        bool match[sl+1][pl+1];
        match[0][0]=true;
        int i;
        for(i=1;i<pl+1;i++)
            if(p[i-1]=='*')
                match[0][i]=true;
            else
                break;
        for(int j=i;j<pl+1;j++)
            match[0][j]=false;
        for(i=1;i<sl+1;i++)
            match[i][0]=false;
        for(int x=1; x<sl+1; x++)
            for(int y=1; y<pl+1; y++) {
                if(s[x-1] == p[y-1]){
                    match[x][y]=match[x-1][y-1];
                } else if(p[y-1]=='?')
                    match[x][y]=match[x-1][y-1];
                else if(p[y-1]=='*')
                    match[x][y]=match[x][y-1] || match[x-1][y];
                else
                    match[x][y]=false;
            }
        return match[sl][pl];
    }
};


// Permutations
class Solution {
public:
    vector<vector<int> > permute(vector<int> &num) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        vector<vector<int> > res;
    	if(num.size()==0) return res;
		if(num.size()==1) { vector<int> a; a.push_back(num[0]); res.push_back(a); return res; }
		int last=num[num.size()-1];
        vector<int> subV=num;
		subV.resize(subV.size()-1);
		vector<vector<int> > subperm=permute(subV);
		for(int i=0; i<subperm.size(); i++) {
			vector<int> v=subperm[i];
			for(int j=0; j<=v.size(); j++) {
				vector<int> vt=v;
				vt.insert(vt.begin()+j, last);
				res.push_back(vt);
			}
		}
		return res;	
    }
};


// longest increasing subsequence copied 
vector<int> find_lis(vector<int> A)
{
  vector<int> m;
  vector<int> p(A.size());
 
  int i, begin, end, mid;
  
  m.push_back(0);
 
  for (i = 1; i < A.size(); i++){
    if (A[m.back()] <= A[i]){
      p[i] = m.back();   // the order here is important
      m.push_back(i);
      continue;
    }
     
    for (begin = 0, end = m.size() -1; begin < end; ){
      mid = (begin + end) / 2;
      if (A[m[mid]] <= A[i])
        begin = mid + 1;
      else
        end = mid;
    }
     
    if (A[i] < A[m[begin]]){
      m[begin] = i;
      if (begin > 0) p[i] = m[begin - 1];
    }
  }
 
  vector<int> result;
  int pos ;
 
  pos = m.back();
  for (i = m.size(); i > 0; i--){
    result.push_back(A[pos]);
    pos =  p[pos];
  }
   
  reverse(result.begin(), result.end());
  return result;
}

		
//  LRU
#ifndef myLRU
#define myLRU
#include <iostream>
#include <map>
using namespace std;

template<class K, class V>
class LRU {
	struct Node {
		K key;
		V value;
		Node* left;
		Node* right;
		Node() {left=right=NULL;}
	};
	public:
	LRU(unsigned c) {cap = c; curSize=0; 
		head=new Node; tail =new Node;
		head->right = tail; tail->left=head;}
	~LRU() {
		Node *t=head;
		while(t) {
			Node* nextT=t->right;
			delete t;
			t=nextT;
		}
	}
	V find(K key, bool &found) {
		if(index.find(key) != index.end()) {
			Node* n = index[key];
			remove(n);
			prepend(n);
			found=true;
			return n->value;
		} else {
			found = false;
			return *(new V);
		}
	}
	void insert(K key, V value) {
		if(index.find(key) != index.end()) { // update the value
			Node *n = index[key];
			remove(n);
			prepend(n);
			n->value = value;
			index[key] = n;
		} else if(curSize < cap) {
			Node *n = new Node;
			n->key = key; n->value = value;
			prepend(n);
			index[key] = n;
			curSize++;
		} else {
			Node *change = tail->left;
			remove(change);
			index.erase(change->key);
			change->key = key;
			change->value = value;
			prepend(change);
			index[key]=change;
		}
	}
	void reportStat() {
		cout<<"current size="<<curSize<<" of capacity "<<cap<<endl;
	}
	private:
	void prepend(Node* n) {
		n->left=head;
		n->right=head->right;
		head->right=n;
		n->right->left=n;
	}
	void remove(Node* n) {
		n->left->right=n->right;
		n->right->left=n->left;
		n->left=NULL; n->right=NULL;
	}
	unsigned curSize;
	unsigned cap;
	Node* head;
	Node* tail;
	map<K, Node*> index;
};
#endif


#ifndef MYHASH
#define MYHASH
#include <iostream>
template <class K, class O>
class Hash {
private:
	struct Node {
		K key;
		O obj;
	};
	int sz;
	unsigned hashFunc(K key);
	static const int PRIME=7919;
	O nullObj;
	Node** nodes;
public:
	Hash() {nodes = new Node*[PRIME]; memset(nodes, 0, PRIME); sz=PRIME;} // choose a large prime
	unsigned size();
	O find(K key);
	bool put(K key, O obj);
	~Hash() {
		for(int i=0; i<PRIME; i++)
			if(nodes[i])
				delete nodes[i];
		delete nodes;
	}
	O end() {return nullObj;}
};
template<class K, class O>
unsigned Hash<K,O>::hashFunc(K key) {
	unsigned kSize=sizeof(K);
	char* ptr=(char*)&key;
	unsigned h=0;
	for(int i=0; i<kSize; i++) {
		h=h*7 + *(ptr+i);
	}
	h=h/sz;
	return h;
}
template<class K, class O>
unsigned Hash<K,O>::size() {
	return sz;
}
template<class K, class O>
O Hash<K,O>::find(K key) {
	unsigned pos=hashFunc(key);
	if(nodes[pos] == NULL)
		return nullObj;
	return nodes[pos]->obj;
}
template<class K, class O>
bool Hash<K,O>::put(K k, O o) {
	unsigned pos=hashFunc(k);
	Node* n = new Node;
	n->key=k; n->obj=o;
	if(nodes[pos] == NULL) {
		nodes[pos] = n;
		std::cout<<"inserted at position "<<pos<<endl;
		return true;
	} else { // linear probing
		return false;
	}
}

#endif			

