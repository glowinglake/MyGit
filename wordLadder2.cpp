class Solution {
public:
void rec(unordered_map<string, int> &visited, unordered_map<string, unordered_set<string> > &list, 
vector<vector<string> > &res, int wid, int curS, int step, string end, 
vector<string> vec, unordered_set<string> &dict){
if(step<=0){
    if(vec.back()==end)
		res.push_back(vec);
	return;
}
string curStr=vec.back();
for(unordered_set<string>::iterator it = list[curStr].begin(); it != list[curStr].end(); it++) {
	if(visited.find(*it)!=visited.end() && visited[*it]==curS+1){
		vec.push_back(*it);
		rec(visited, list, res, wid, curS+1, step-1, end, vec, dict);
		vec.pop_back();
		}
	}
}
int findLadderss(unordered_map<string, int> &visited, unordered_map<string, unordered_set<string> > &list, string start, string end, unordered_set<string> &dict) {
int step=1;
queue<string> q;
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
if(found)
break;
q.push(string(""));
level++;
continue;
}
//visited[s]=level;
visited.insert(pair<string, int>(s, level));
for(unordered_set<string>::iterator it = list[s].begin(); it != list[s].end(); it++) {
	if(*it == end){
		found=true;
		//return level;
	} else if(visited.find(*it) == visited.end())
		q.push(*it);
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
unordered_map<string, int> visited;
int level=findLadderss(visited, list, start, end,dict); vector<vector<string> > res;
if(level==-1)
return res;
vector<string> vec;
vec.push_back(start);

rec(visited, list, res, end.length(), 1, level, end, vec, dict);
return res;
}
};

