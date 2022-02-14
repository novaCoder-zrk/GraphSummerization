package sweg.algorithm;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import jdk.jshell.spi.ExecutionControl;
import sweg.SummaryGraphModule;

import java.time.Instant;
import java.util.ArrayList;

import static java.lang.Integer.min;
import static java.lang.Long.min;

public class SWeG extends SummaryGraphModule {

    private static int CMDFAULT = -1;
    private int T = 20; // number of iterations
    private double eps = 0; // error bound

    //input 用来存放原图的邻接链表  每个点的临接点用set存储 （之后为了效率，将其转化为数组）
    private ObjectArrayList<IntOpenHashSet> input = new ObjectArrayList<>(0);
    // 用来存原图的邻接链表，由input转化而来
    private IntArrayList[] adjList;


    //用来存放超点
    private IntArrayList[] S;
    private int[] deg;
    //Super2Node 保存了每个超点的所有邻接点
    //Node2Super 记录了每个点所邻接的超点
    private Int2IntOpenHashMap[] Super2Node, Node2Super;

    //weightMap
    private ObjectArrayList<ObjectArrayList<Double>> weightSet = new ObjectArrayList<>();
    private Int2IntOpenHashMap[] SuperEdge2Weight;
    int wn = 0;

    private ObjectArrayList<Int2DoubleOpenHashMap> weightMap = new ObjectArrayList<>();
    private int[] h;//用来储存随机排序pi
    private int[] shingles;
    private ObjectArrayList<IntArrayList> groups;//存放divide分组结果
    private int threshold = 500, max_step = 10;//在divide步骤中使用到的数据
    private int iter_max_step = 0;
    private Instant _start = Instant.now();

    private long uniq = 0;
    private long[] inGroup;
    private int[] val = new int[threshold];

    private int[] candidates; private int sn;//度大于1的点会被归为candidates
    private long[][] sep; private int sepn;//不明确

    public SWeG(boolean directed){
        super(directed);
        if(directed){
            try {
                throw new ExecutionControl.NotImplementedException("Directed version is NOT_IMPLEMENTED");
            } catch (ExecutionControl.NotImplementedException e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }
    }

    @Override
    public void addVertex(int idx) {
        super.addVertex(idx);

        //modified
        //input.add(new IntOpenHashSet(1));
        weightMap.add(new Int2DoubleOpenHashMap());
    }

    @Override
    public void processEdge(final int src, final int dst, final boolean add) {
        if(add){
            input.get(src).add(dst);
            input.get(dst).add(src);
        }else{
            input.get(src).remove(dst);
            input.get(dst).remove(src);
        }
    }

    @Override
    public void processWeightedEdge(final int src, final int dst,final double weight,final boolean add){
        if(add){
            weightMap.get(src).put(dst,weight);
            weightMap.get(dst).put(src,weight);
        }else{
            weightMap.get(src).remove(dst);
            weightMap.get(dst).remove(src);
        }
    }

    @Override
    public ObjectArrayList<Int2DoubleOpenHashMap> returnWeight(){
        return weightMap;
    }

    private void addSuperEdge(int u, int v){
        // add directed edge u -> v
        if(u == v && getSize(u) == 1) return;
        P.addEdge(u, v);

        //不再需要Cm了
        /*for(int _u: S[u]){
            for(int _v: S[v]){
                if(_u == _v) continue;
                if(!adjList[_u].contains(_v)){
                    Cm.addEdge(_u, _v);
                }
            }
        }*/
    }

    private int getSize(final int Sv){
        return S[Sv].size();
    }

    private long getPi(final int Su, final int Sv){
        long pi = getSize(Su); pi *= getSize(Sv);
        if(Su == Sv){
            pi -= getSize(Sv);
            pi /= 2;
        }
        return pi;
    }

    private Int2IntOpenHashMap getSuperDegree(int Sv){
        Int2IntOpenHashMap SuperDeg = new Int2IntOpenHashMap();
        for(Int2IntMap.Entry v: Super2Node[Sv].int2IntEntrySet()){
            SuperDeg.addTo(V.getInt(v.getIntKey()), v.getIntValue());
        }
        return SuperDeg;
    }

    private long getCostInner(int v, int vSize, Int2IntOpenHashMap Nv, boolean add){
        long cost = 0;
        //遍历每个neighbor  相邻的supernode
        for(Int2IntMap.Entry nbr: Nv.int2IntEntrySet()){
            long pi, edgeCount;
            if(v == nbr.getIntKey()){
                pi = vSize; pi *= (vSize - 1);
                edgeCount = nbr.getIntValue();
                pi /= 2; edgeCount /= 2;
            }else{
                pi = vSize; pi *= getSize(nbr.getIntKey());
                edgeCount = nbr.getIntValue();
            }
            cost += min(pi - edgeCount + 1, edgeCount);

            if(add && (pi - edgeCount + 1) < edgeCount) {
                addSuperEdge(v, nbr.getIntKey());
                //添加超边
                if(v <= nbr.getIntKey()) //防止重复加边
                    processWeightSet(v,nbr.getIntKey());
            }

        }
        return cost;
    }


    private long getCost(int v){

        return getCostInner(v, getSize(v), getSuperDegree(v), false);
    }

    private long getMergeCost(int u, int v){
        Int2IntOpenHashMap Nw = getSuperDegree(u);
        for(Int2IntMap.Entry nbr: getSuperDegree(v).int2IntEntrySet()){
            Nw.addTo(nbr.getIntKey(), nbr.getIntValue());
        }
        Nw.addTo(v, Nw.getOrDefault(u, 0));
        Nw.remove(u);
        return getCostInner(v, getSize(u) + getSize(v), Nw, false);
    }

    private double saving(int A, int B){
        long before = getCost(A) + getCost(B);
        long after = getMergeCost(A, B);
        long pi = getSize(A); pi *= getSize(B);
        long edgeCount = 0;
        for(int v: S[A]){
            edgeCount += Node2Super[v].getOrDefault(B, 0);
        }
        before -= min(edgeCount, pi - edgeCount + 1);
        //System.out.println(A + " vs " + B + " : " + before + " -> " + after);
        return 1 - (after / (double)before);
    }

    private void divideInner(int step){
        //一个递归程序，如果分到的某一组中顶点的数量大于threshold，会进行第二次分组，
        // 直到max_step次递归时，如果依旧是大于threshold，直接threshold个顶点一组
        int me = step % 2;
        if(step == max_step){
            //到最后，如果每个组的大小还大于threshold，直接每threshold一组
            for(int ii=0;ii<sepn;ii++){
                //sep 的一个元素是一个long类型 long是64位 int 32位 一个long可以记录一段区间的开头结尾
                int s = (int)(sep[me][ii] >> 32), e = (int)(sep[me][ii] & 0x7FFFFFFFL);
                for(int i=s;i<=e;i+=threshold){
                    //以每threshold个顶点为一组，最后一组不足threshold个点
                    if(i + (threshold - 1) <= e){
                        groups.add(new IntArrayList(candidates, i, threshold));
                    }else{
                        groups.add(new IntArrayList(candidates, i, (e-i+1)));
                    }
                }
            }
            return;
        }
        //生成排序 pi h用来储存排序
        h[0] = 1;
        for (int i = 1; i < n; i++) {
            h[i] = i+1;
            int randIdx = randInt(0, i);
            h[i] = h[randIdx];
            h[randIdx] = i+1;
        }

        //sepn说明了有几个要分组的区间
        // 可能有多组在经过前一轮分组后，每个组内的元素数大于threshold
        int nxt_sepn = 0;
        for(int ii=0;ii<sepn;ii++){
            //确定每一组的起点终点
            int s = (int)(sep[me][ii] >> 32), e = (int)(sep[me][ii] & 0x7FFFFFFFL);

            for(int i=s;i<=e;i++){
                //遍历每一个组的candidate
                final int A = candidates[i];//将candidates 的single依次算出
                int minHash = 0x7FFFFFFF;
                //计算超点的shingle
                for(int v: S[A]){
                    minHash = min(minHash, h[v]);
                    for(int u: adjList[v]){
                        minHash = min(minHash, h[u]);
                    }
                }
                shingles[A] = minHash;
            }
            //快排 将candidates 按照shingles进行排序
            IntArrays.parallelQuickSort(candidates, s, e + 1, new IntComparator() {
                @Override
                public int compare(int i, int i1) {
                    return Integer.compare(shingles[i], shingles[i1]);
                }
            });
            //将candidate按照 shingles排过了
            int prv = s;
            for(int i=s+1;i<=e;i++){
                if(shingles[candidates[i]] != shingles[candidates[i-1]]){
                    if(i - prv <= threshold){
                        //如果这一个组的规模小于 threshold 直接归为一组
                        groups.add(new IntArrayList(candidates, prv, i-prv));
                    }else{
                        //进行第二步 再分 记录下起始位置，记录下终止位置 ，放到一个long数组里
                        sep[1-me][nxt_sepn++] = (((long)prv) << 32) + (i-1);
                    }
                    prv = i;//记录上一个分解
                }
            }
            if((e + 1) - prv <= threshold){
                groups.add(new IntArrayList(candidates, prv, (e+1) - prv));
            }else{
                sep[1-me][nxt_sepn++] = (((long)prv) << 32) + e;
            }
        }
        if(nxt_sepn > 0){
            sepn = nxt_sepn;
            divideInner(step + 1);
        }
    }

    private void divide(){
        // Input: input graph G = (V, E), current supernodes S
        // Output: disjoint groups of supernodes (shingle)

        //度大于0的点被归为candidates  deg是超点的度？？？

        sn = 0;
        for(int i=0;i<n;i++){
            if(deg[i] > 0) candidates[sn++] = i;//通过判断deg是否大于零，来筛选超点
        }
        //sepn 指明分段区间的个数
        sepn = 1;
        sep = new long[2][(sn + threshold + threshold) / threshold];//一个滚动数组，用于记录进行分组的区间
        sep[0][0] = sn-1;
        //以threshold为界，每一组最多有threshold个点
        groups = new ObjectArrayList<>(n / threshold);
        //进一步进行
        divideInner(0);
        /*
        int _sn = 0;
        IntArrayList checky = new IntArrayList(n);
        for(int i=0;i<n;i++){
            checky.add(0);
        }
        for(IntArrayList _Q: groups){
            _sn += _Q.size();
            if(_Q.size() > threshold) System.out.println("@");
            for(int v: _Q){
                if(checky.getInt(v) > 0){ System.out.println("?"); }
                checky.set(v, 1);
            }
        }
        System.out.println(sn + " " + _sn);
         */
    }

    private void merge(int iter){
        //将分组后的超点取出
        for(IntArrayList _Q: groups){
            sn = 0;//暂时不知到
            IntArrayList Q = new IntArrayList(_Q);
            for(int q: Q) {
                //Super2Node 为每一个超点建立一个map
                Super2Node[q] = new Int2IntOpenHashMap(0);
            }
            for(int q: Q){
                //遍历一组中的每一个超点
                //但遍历到一个普通点，就检查是否有Node2Super的HashMap
                for(int v: S[q]){
                    //遍历SuperNode中的每一个点 v
                    if(Node2Super[v] == null){
                        //如果Node2Super 为空
                        candidates[sn++] = v; //将普通点加入candidate
                        Node2Super[v] = new Int2IntOpenHashMap(0);
                    }
                    for(int u: adjList[v]){
                        //遍历点v的所有邻接点
                        if(Node2Super[u] == null){
                            candidates[sn++] = u;
                            Node2Super[u] = new Int2IntOpenHashMap(0);
                        }
                        // u作为超点q的邻接点加入Super2Node
                        // 超点q是点u所邻接的超点，加入Node2Super
                        Super2Node[q].addTo(u, 1);
                        Node2Super[u].addTo(q, 1);
                        //V储存的是对应的超点号
                        if(Super2Node[V.getInt(u)] == null) Node2Super[v].addTo(V.getInt(u), 1);

                    }
                }
            }
            int sz = Q.size();
            while(sz > 1){
                uniq += 1;
                // pick and remove random supernode A from Q
                int randIdx = randInt(0, sz-1);
                int A = Q.getInt(randIdx);
                Q.set(randIdx, Q.getInt(sz-1));
                Q.popInt(); sz -= 1;

                for(int i=0;i<sz;i++){
                    //inGroup是干嘛的？
                    inGroup[Q.getInt(i)] = uniq * (threshold + 1) + (i + 1);
                    val[i] = 0;
                }
                for(Int2IntMap.Entry v: Super2Node[A].int2IntEntrySet()){
                    for(Int2IntMap.Entry U: Node2Super[v.getIntKey()].int2IntEntrySet()){
                        //v.getIntKey() 是超点A的neighbor
                        if(inGroup[U.getIntKey()] / (threshold + 1) == uniq){
                            //为了保证找到的点都在同一个 group 里， A B 属于同一个group
                            // Compute sum of min(w(A, v), w(B, v))
                            // val储存的是 sum of min
                            val[(int)(inGroup[U.getIntKey()] % (threshold + 1)) - 1] += min(v.getIntValue(), U.getIntValue());
                        }
                    }
                }
                double maxSuperJaccard = -1.0;
                int argMax = -1;
                for(int i=0;i<sz;i++){
                    double superJaccard = val[i] / (double) (deg[A] + deg[Q.getInt(i)] - val[i]);
                    if(maxSuperJaccard < superJaccard){
                        maxSuperJaccard = superJaccard;
                        argMax = i;
                    }
                }
                int B = Q.getInt(argMax);
                Q.set(argMax, Q.getInt(sz-1));
                Q.set(sz-1, B);
                //找到jaccard相似系数最大的两个超点 A,B
                double threshold = (iter < T) ? (1.0 / (double) (1 + iter)) : 0;
                double result = saving(A, B);
                if(result >= threshold){
                    if(getSize(A) + Super2Node[A].size() < getSize(B) + Super2Node[B].size()){
                        int tmp = A; A = B; B = tmp;
                    }
                    // merge A and B: A <- B
                    for(int v: S[B]){
                        S[A].add(v);
                        V.set(v, A);  //破案了 在虚点号为v的位置，V[v]储存所在超点的编号
                    }
                    deg[A] += deg[B];//deg只出储存超点的度
                    deg[B] = 0;      //非超点的度被清零
                    for(Int2IntMap.Entry v: Super2Node[B].int2IntEntrySet()){
                        Node2Super[v.getIntKey()].addTo(A, v.getIntValue());
                        Node2Super[v.getIntKey()].remove(B);
                        Super2Node[A].addTo(v.getIntKey(), v.getIntValue());
                    }

                    S[B] = null;
                    Super2Node[B] = null;
                    Q.set(sz-1, A);
                }
            }
            // Clear Super2Node and Node2Super
            for(int q: _Q){
                Super2Node[q] = null;
            }
            for(int i=0;i<sn;i++){
                Node2Super[candidates[i]] = null;
            }
        }
        groups.clear();
    }
    private void processWeightSet(int su,int sv){
        //将 超边su sv 上的边的属性记录在集合中
        if(SuperEdge2Weight[su] == null) SuperEdge2Weight[su] = new Int2IntOpenHashMap();
        if(SuperEdge2Weight[sv] == null) SuperEdge2Weight[sv] = new Int2IntOpenHashMap();
        SuperEdge2Weight[su].put(sv,wn);
        SuperEdge2Weight[sv].put(su,wn);
        wn++;
        weightSet.add(new ObjectArrayList<Double>());

        if(su == sv){
            int size  = S[su].size();
            for(int i = 0;i < size;i++)
                for(int j = i+1;j < size;j++){
                    int a = S[su].getInt(i);
                    int b = S[su].getInt(j);
                    //在这里设置default
                    double w = weightMap.get(a).getOrDefault(b,-1);
                    weightSet.get(wn-1).add(w);
                }
        }


        if(su < sv){
            int size_1 = S[su].size();
            int size_2 = S[sv].size();
            for(int i = 0;i < size_1;i++)
                for(int j = 0;j < size_2;j++){
                    int a = S[su].getInt(i);
                    int b = S[sv].getInt(j);
                    //在这里设置default
                    double w = weightMap.get(a).getOrDefault(b,CMDFAULT);
                    weightSet.get(wn-1).add(w);
                }
        }
    }
    @Override
    public ObjectArrayList<Double> returnWeightSet(int su,int sv){
        int a = SuperEdge2Weight[su].get(sv);
        return weightSet.get(a);
    }

    private void encode(){
        for(int i=0;i<n;i++){
            for(int j: adjList[i]){
                //此时Super2Node应该已经都被清空了
                if(Super2Node[V.getInt(i)] == null) Super2Node[V.getInt(i)] = new Int2IntOpenHashMap(0);
                Super2Node[V.getInt(i)].addTo(V.getInt(j), 1);
            }
        }
        // add superedges (P and Cm) 在需要Cm了
        for(int i=0;i<n;i++){
            //getsize 返回超点的大小
            if(deg[i] > 0) getCostInner(i, getSize(i), Super2Node[i], true);
        }
        // add correction edges (Cp)
        for(int i=0;i<n;i++){
            int Si = V.getInt(i);
            for(int v: adjList[i]){
                int Sv = V.getInt(v);
                if(!P.getNeighbors(Si).contains(Sv)){
                    Cp.addEdge(i, v);
                }
            }
        }
    }

    @Override
    public void processBatch(){
        System.out.println("|V|: " + n);
        //初始化邻接链表 此时一共有n个节点 从set换成了IntArrayList 是为了效率？
        adjList = new IntArrayList[n];

        /*for(int i=0;i<n;i++){
            adjList[i] = new IntArrayList(input.get(i));
            input.set(i, null);

        }*/
        //input = null;
        for(int i = 0;i < n;i++){
            adjList[i] = new IntArrayList(weightMap.get(i).keySet());
        }

        // Overview of SWeG
        // Input: input graph G = (V, E), iterations T, error bound eps
        // Output: Summary graph G' = (S, P), corrections C+, C-

        S = new IntArrayList[n];

        // for divide
        h = new int[n];//不清楚
        shingles = new int[n];//记录shingles
        candidates = new int[n+1];//用于筛选出超点
        //Super2Node 和 Node2Super 分别记录 超点所邻接的点 点所邻接的超点
        Super2Node = new Int2IntOpenHashMap[n];
        Node2Super = new Int2IntOpenHashMap[n];
        inGroup = new long[n];
        deg = new int[n];

        //for weight
        SuperEdge2Weight = new Int2IntOpenHashMap[n];


        // initialize supernodes S to {{v}: v \in V}
        for(int i=0;i<n;i++){
            S[i] = new IntArrayList(1);
            S[i].add(i);
            Super2Node[i] = null;
            Node2Super[i] = null;
            inGroup[i] = 0;//inGroup暂不清楚
            deg[i] = adjList[i].size();//记录每个SuperNode的度
            SuperEdge2Weight[i] = null;
        }

        for(int iter=1;iter<=T;iter++){
            // divide S into disjoint groups
            divide();
            // merge some supernodes within each group
            merge(iter);
        }
        // encode edges E into superedges P and corrections C
        encode();
        // We only implemented lossless version of SWeG,
        // since our proposed algorithm MoSSo only supports lossless summarization.
        // System.out.println(uniq * (threshold + 1) + (threshold + 1));
    }
}
