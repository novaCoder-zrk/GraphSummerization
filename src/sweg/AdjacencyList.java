package sweg;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import sweg.util.Int2IntHashMap;
// 邻接链表 用一个数组 adjList 每个数组的元素是一个hashmap
// adjList.get(src).get(dst) 代表了从src到dst的连接边的数量 用来存储超点
public class AdjacencyList {//邻接链表
    private int n = 0;
    //一个list 每个元素是一个hashmap
    private ObjectArrayList<Int2IntHashMap> adjList = new ObjectArrayList<Int2IntHashMap>(0);

    public AdjacencyList(){}//初始化

    public ObjectArrayList<Int2IntHashMap> getAdjList(){
        return adjList;
    }//返回邻接链表

    public void expand(){
        //扩展，多创建一个hashmap
        n += 1;
        Int2IntHashMap newMap = new Int2IntHashMap(0);
        newMap.defaultReturnValue(0);
        adjList.add(newMap);
    }

    public void updateEdge(final int src, final int dst, final int incr){
        assert(src < n && dst < n && src > -1 && dst > -1);//src 和 dst 都 <n >-1
        //adjList [src]  [dst] +1
        int prevValue = adjList.get(src).addTo(dst, incr);
        if(prevValue == -incr) adjList.get(src).remove(dst);
    }

    public void addEdge(final int src, final int dst) {
        updateEdge(src, dst, 1);
    }

    public void deleteEdge(final int src, final int dst){
        updateEdge(src, dst, -1);
    }

    public int getEdgeCount(final int src, final int dst){
        return adjList.get(src).getOrDefault(dst, 0);
    }

    public IntSet getNeighbors(final int src){
        return adjList.get(src).keySet();
    }

    public Int2IntMap.FastEntrySet getNeighborsAndWeights(final int src){
        return adjList.get(src).int2IntEntrySet();
    }

    public int getDegree(final int src) { return adjList.get(src).size(); }
}
