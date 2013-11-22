//implementation of LRU

template <T>
class LRU{
	DoubleLinkList *dlink;
	T cache[size];
	HashMap *cacheMap;
	
	T Get(T* item){
		if(cacheMap->find(item)){
			//move the item to the first in list
			dlink->MoveFirst(item);	
			return T;
			}
		else{
			return NULL;
			}
	}
	bool Put(T* item){
		if((cahceMap->find(item)){
			//already in the cache
			return true;
			}
		for(int i=0;i<size;i++)
	
	}