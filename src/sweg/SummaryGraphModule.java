package sweg;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import jdk.jshell.spi.ExecutionControl;

import java.util.concurrent.ThreadLocalRandom;

public class SummaryGraphModule implements IncrementalInterface{
    //1、原图是怎么存的？ 是通过子类的重构，依旧是邻接链表
    // Core data structure
    protected IntArrayList V = new IntArrayList(0);//存放的是对应点的超点号
    protected AdjacencyList P = new AdjacencyList();
    protected AdjacencyList Cp = new AdjacencyList();
    //protected AdjacencyList Cm = new AdjacencyList();
    protected int n = 0;


    // helper
    private final ThreadLocalRandom random;
    private IntArrayList idxs = new IntArrayList(0);//存放实际顶点号
    private Int2IntOpenHashMap vmap = new Int2IntOpenHashMap(0); // 从真实编号向虚的编号的映射

    // argument
    protected boolean directed;

    public SummaryGraphModule(boolean directed){
        this.random = ThreadLocalRandom.current();
        this.directed = directed;
    }

    public void addVertex(int idx) {
        vmap.addTo(idx, n);//给顶点重新编号
        idxs.add(idx);// 保存实际顶点号
        V.add(n);//V 保存虚标号
        n += 1;
        P.expand(); Cp.expand();// Cm.expand();
    }

    @Override
    public void processAddition(final int src, final int dst){
        //加边操作
        //首先确定，是否有顶点没有加入顶点集，如果有新的顶点加入，为其分配新的编号，通过vmap储存映射关系
        //processEdge 将原图中的边插入
        if(!vmap.containsKey(src)) addVertex(src);
        if(!vmap.containsKey(dst)) addVertex(dst);
        int[] vIdx = {vmap.getOrDefault(src, -1), vmap.getOrDefault(dst, -1)};
        processEdge(vIdx[0], vIdx[1], true);
    }
    //weighted 版本
    public void processWeightedAddition(final int src,final int dst,final double weight){
        if(!vmap.containsKey(src)) addVertex(src);
        if(!vmap.containsKey(dst)) addVertex(dst);
        int[] vIdx = {vmap.getOrDefault(src, -1), vmap.getOrDefault(dst, -1)};
        processWeightedEdge(vIdx[0],vIdx[1],weight,true);
    }

    public void processWeightedEdge(final int src, final int dst, final double weight,final boolean add) {
        //这个函数要在下一层被重构
        try {
            throw new ExecutionControl.NotImplementedException("processWeightedEdge: NOT_IMPLEMENTED");
        } catch (ExecutionControl.NotImplementedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }


    @Override
    public void processDeletion(final int src, final int dst){
        //删除边
        if(!vmap.containsKey(src)) addVertex(src);
        if(!vmap.containsKey(dst)) addVertex(dst);
        int[] vIdx = {vmap.getOrDefault(src, -1), vmap.getOrDefault(dst, -1)};
        processEdge(vIdx[0], vIdx[1], false);
    }

    public void processBatch(){
        //这个在下一层也要被重构
        try {
            throw new ExecutionControl.NotImplementedException("processBatch: NOT_IMPLEMENTED");
        } catch (ExecutionControl.NotImplementedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }


    public ObjectArrayList<Int2DoubleOpenHashMap> returnWeight(){
        try {
            throw new ExecutionControl.NotImplementedException("returnWeight: NOT_IMPLEMENTED");
        } catch (ExecutionControl.NotImplementedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
    }

    protected double randDouble(){
        return random.nextDouble();
    }

    protected int randInt(final int from, final int to){
        // return generated random number in [from, to] (close interval)//闭区间
        return from + random.nextInt(to - from + 1);
    }

    protected long randLong(final long from, final long to){
        // return generated random number in [from, to] (close interval)
        return from + random.nextLong(to - from + 1);
    }

    public void processEdge(final int src, final int dst, final boolean add) {
        //这个函数要在下一层被重构
        try {
            throw new ExecutionControl.NotImplementedException("processEdge: NOT_IMPLEMENTED");
        } catch (ExecutionControl.NotImplementedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
    public ObjectArrayList<Double> returnWeightSet(int su,int sv){
        try {
            throw new ExecutionControl.NotImplementedException("returnWeightSet: NOT_IMPLEMENTED");
        } catch (ExecutionControl.NotImplementedException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
    }

    public Tuple<IntArrayList, AdjacencyList[], IntArrayList> getCompressedResult(){
        //return new Tuple<IntArrayList, AdjacencyList[], IntArrayList>(V, new AdjacencyList[]{P, Cp, Cm}, idxs);
        return new Tuple<IntArrayList, AdjacencyList[], IntArrayList>(V, new AdjacencyList[]{P, Cp}, idxs);
    }
}
