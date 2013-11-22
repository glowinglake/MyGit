#include <iostream>
using namespace std;
bool compute24(int nums[],int length, int target){
	int i,j;
	for(i=0;i<length;i++)
		for(j=0;j<length;j++)
			if(i!=j){
				int val;
				if(length==2){
					//see the result to determine if there is a solution
					if(nums[i]+nums[j]==target)
						return true;
					if(nums[i]-nums[j]==target)
						return true;
					if(nums[i]*nums[j]==target)
						return true;
					if(nums[i]/nums[j]==target)
						return true;
						
					return false;
					}
				else if(length>2){
					//recursion
					int newnums[length-1];
					int newi=0;
					//put the other numbers in the new array
					for(int t=0;t<length;t++){
						if(t!=i && t!=j){
							newnums[newi++]=nums[t];
						}
					newnums[length-2]=nums[i]+nums[j];
					if(compute24(newnums,length-1,target))
						return true;
					
					newnums[length-2]=nums[i]-nums[j];
					if(compute24(newnums,length-1,target))
						return true;
					
					newnums[length-2]=nums[i]*nums[j];
					if(compute24(newnums,length-1,target))
						return true;
					
					newnums[length-2]=nums[i]/nums[j];
					if(compute24(newnums,length-1,target))
						return true;
					
					return false;
					}
				}
			}
}

int main(){
	int nums[4]={1,2,3,5};
	
	if(compute24(nums,4,24))
		cout<<"found\n";
}