// PreOrder Option #2 Non-recursive static method
public static IEnumerable<T> Preorder<T>(Node<T> root)
{
    var stack = new Stack<Node<T>>();
    stack.Push(root);

    while (stack.Count > 0)
    {
        var node = stack.Pop();
        yield return node.Data;
        if (node.RightChild != null)
            stack.Push(node.RightChild);
        if (node.LeftChild != null)
            stack.Push(node.LeftChild);
    }
}

// non recursive post order using visited information
nonRecursivePostorder(rootNode)
  nodeStack.push(rootNode)
  while (! nodeStack.empty())
    currNode = nodeStack.peek()
    if ((currNode.left != null) and (currNode.left.visited == false))
      nodeStack.push(currNode.left)
    else 
      if ((currNode.right != null) and (currNode.right.visited == false))
        nodeStack.push(currNode.right)
      else
        print currNode.value
        currNode.visited := true
        nodeStack.pop()

// non recursive post order without visited information
void postOrderTraversalIterative(BinaryTree *root) {
  if (!root) return;
  stack<BinaryTree*> s;
  s.push(root);
  BinaryTree *prev = NULL;
  while (!s.empty()) {
    BinaryTree *curr = s.top();
    if (!prev || prev->left == curr || prev->right == curr) {
      if (curr->left)
        s.push(curr->left);
      else if (curr->right)
        s.push(curr->right);
    } else if (curr->left == prev) {
      if (curr->right)
        s.push(curr->right);
    } else {
      cout << curr->data << " ";
      s.pop();
    }
    prev = curr;
  }
}		
		
// non rec inorder using visited
void foo(node *root){
	stack s;
	s.push(root);
	while(s.empty==false){
		node *n=s.peek();
		if(n.visited){
			printf(n.data);
			s.pop();
			if(n.right)
				s.push(n.right);
		}
		else{
			n.visited=true;
			if(n.left)
				s.push(n.left);
		}
	}
}
// non rec inorder without using visited
// 1
void in_order_traversal_iterative(BinaryTree *root) {
  stack<BinaryTree*> s;
  BinaryTree *current = root;
  while (!s.empty() || current) {
    if (current) {
      s.push(current);
      current = current->left;
    } else {
      current = s.top();
      s.pop();
      cout << current->data << " ";
      current = current->right;
    }
  }
}
// 2
void iterative_inorder(Node * node) {
Stack s; // Initially empty
while(1) {
	while(node) {
		s.push(node);
		node = node->left;
		}
	if(s.isempty()) break;
	else {
		Node * temp = s.pop();
		printf("%d ", temp->value);
		node = temp->right;
		}
	}
}


//DFS
So the basic structure will look something like this:
dfs(node start) {
 stack s;
 s.push(start);
 while (s.empty() == false) {
  top = s.top();
  s.pop();
  mark top as visited;
  check for termination condition
  add all of top's unvisited neighbors to the stack.
  //mark top as not visited;
 }
}
Alternatively we can define the function recursively as follows:
dfs(node current) {
 mark current as visited;
 visit all of current's unvisited neighbors by calling dfs(neighbor)
 mark current as not visited;
}

//BFS
The basic structure of a breadth first search will look this:
void bfs(node start) {
 queue s;
 s.push(start);
 while (s.empty() == false) {
  top = s.front();
  s.pop();
  mark top as visited;
check for termination condition (have we reached the node we want to?) add all of top's unvisited neighbors to the stack.
 }
}

//tree succeessor
TREE-SUCCESSOR.x/
1 if x:right ก่ NIL
2 return TREE-MINIMUM.x:right/
3 y D x:p
4 while y ก่ NIL and x == y:right
5 x D y
6 y D y:p
7 return y

//tree predecessor symmetry of above




//regular expression . *
class Solution {
public:
    bool isMatch(const char *s, const char *p) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function    
        if(!s || *s=='\0'){
			if(!p || *p=='\0')
				return true;
			else if((p+1) && *(p+1)=='*')
				return isMatch(s, p+2);
			else
				return false;
		}
		else if(!p || *p=='\0')
			return false;
		//now both have content
		else if(*s != *p && *p != '.'){
			if(*(p+1) !='*')
				return false;
			else
				return isMatch(s, p+2);
			}
		else if(*s != *p && *p == '.'){
			if(*(p+1) != '*')
				return isMatch(s+1, p+1);
			else
				return isMatch(s, p+2) || isMatch(s+1, p);
			}
		else if(*s == *p){
			if(*(p+1) != '*')
				return isMatch(s+1, p+1);
			else
				return isMatch(s+1, p) || isMatch(s, p+2);
			}
    }
};



// Maximum_subarray_problem
// dp using O(n) space
int foo(int a[], int len)
{
	// max subarray ending at i = max[i]
	int max[len]={0};
	max[0]=a[0];
	for(int i=1;i<len;i++){
		if(max[i-1]<=0)
			max[i]=a[i];
		else
			max[i]=a[i]+max[i-1];
		}
	//search in max and find the max value
	int res=-10000;
	for(int j=0;j<len;j++)
		if(max[i]>res)
			res=max[i];
	
	return res;
}
// modified using O(1) space
int foo(int a[], int len)
{
	// max subarray ending at i = max[i]
	int max, currentMax;
	max=a[0];
	currentMax=a[0];
	for(int i=1;i<len;i++){
		currentMax=currentMax<0?a[i]:a[i]+currentMax;
		max=currentMax>max?currentMax:max;
		}
	return max;
}

// smart pointer
1 template <class T> class SmartPointer {
2 public:
3 SmartPointer(T * ptr) {
4 ref = ptr;
5 ref_count = (unsigned*)malloc(sizeof(unsigned));
6 *ref_count = 1;
7 }
8 SmartPointer(SmartPointer<T> & sptr) {
9 ref = sptr.ref;
10 ref_count = sptr.ref_count;
11 ++*ref_count;
12 }
13 SmartPointer<T> & operator=(SmartPointer<T> & sptr) {
14 if (this != &sptr) {
15 ref = sptr.ref;
16 ref_count = sptr.ref_count;
17 ++*ref_count;
18 }
19 return *this;
20 }
21 ~SmartPointer() {
22 --*ref_count;
23 if (*ref_count == 0) {
24 delete ref;
25 free(ref_count);
26 ref = NULL;
27 ref_count = NULL;
28 }
29 }
30 T getValue() { return *ref; }
31 protected:
32 T * ref;
33 unsigned * ref_count;
34 };



// BINARY SEARCH
void searchRange(int A[], int n, int target) {
        // Start typing your C/C++ solution below
        // DO NOT write int main() function
        int s=-1,e=-1;
        int bs=0,be=n-1;
        while(1){
            int mid=(bs+be)/2;
            if(bs==be && A[bs]!=target)
            	break;
           	else        if(A[mid]<target)
                bs=mid+1;
            else if(A[mid]>target)
                be=mid-1;
            else if(mid>0 && A[mid-1]==target)
                be=mid-1;
            else{
                s=mid;
                break;
            }
        
        }
        bs=0; be=n-1;
        while(1){
            int mid=(bs+be)/2;
            if(bs==be && A[bs]!=target)
                break;
               else if(A[mid]<target)
                bs=mid+1;
            else if(A[mid]>target)
                be=mid-1;
            else if(mid<n-1 && A[mid+1]==target)
                bs=mid+1;
         
            else{
                e=mid;
                break;
            }
        }
        cout<<s<<" "<<e;
        /*
        vector<int> res;
        res.insert(res.end(),s);
        res.insert(res.end(),e);
        return res;
        */
}