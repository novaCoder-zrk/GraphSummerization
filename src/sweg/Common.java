package sweg;

import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;

public class Common {

    public Common(){}

    private static int[] parseEdge(String line, String delim){
        String[] tokens = line.split(delim);
        try {
            int src = Integer.valueOf(tokens[0]);
            int dst = Integer.valueOf(tokens[1]);
            int add = Integer.valueOf(tokens[2]);
            return new int[]{src, dst, add};
        }catch(Exception e){
            return null;
        }
    }
    private static double[] parseDoubleEdge(String line,String delim){
        String[] tokens = line.split(delim);
        try{
            double src = Double.valueOf(tokens[0]);
            double dst = Double.valueOf(tokens[1]);
            double weight = Double.valueOf(tokens[2]);
            return new double[]{src,dst,weight};
        }catch (Exception e){
            return null;
        }
    }

    public static int execute(final SummaryGraphModule module, final String inputPath, final String delim) throws IOException{
        //开始执行 设定时间 打开文件
        int count = 0; //
        long start = System.currentTimeMillis();
        BufferedReader br = new BufferedReader(new FileReader(inputPath));

        while(true){
            //读入每一行进行处理
            final String line = br.readLine();
            if(line == null) break;
            final int[] edge = parseEdge(line, delim);//对一行的输入进行解析
            if(edge == null) break;
            count += 1; //对总的边数进行计算
            //if(edge[2] > 0)
            module.processWeightedAddition(edge[0], edge[1],edge[2]); //加边操作
            ///else module.processDeletion(edge[0], edge[1]); //减边操作
        }
        //问题到目前还没看见在哪里确定边的关系
        module.processBatch(); //批处理

        br.close();
        long end = System.currentTimeMillis();
        System.out.println("Execution time: " + (end - start) / 1000.0 + "s.");//输出时间
        return count;
    }

    public static void writeOutputs(final SummaryGraphModule module, final String outputPath, final int edgeCount) throws IOException{
        /*
        Output:
        line 1: |V|, |E|
        line 2 ~ |V|+1: (real_index, membership(block) info)
        line |V|+2 ~ |V|+|P|+1: superedges in P
        line |V|+|P|+2: -1 -1
        line |V|+|P|+3 ~ |V|+|P|+|Cp|+2: edges in Cp
        line |V|+|P|+|Cp|+3: -1 -1
        line |V|+|P|+|Cp|+4 ~ |V|+|P|+|Cp|+|Cm|+3: edges in Cm
        line |V|+|P|+|Cp|+|Cm|+4: -1 -1
        */
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
        Tuple<IntArrayList, AdjacencyList[], IntArrayList> compressedResult = module.getCompressedResult();
        int n = compressedResult.third.size();
        int[] summaryCount = {0, 0, 0};
        //List<String> listNames = Arrays.asList("P", "C_plus", "C_minus");
        List<String> listNames = Arrays.asList("P", "C_plus");
        ListIterator<String> it = listNames.listIterator();
        IntArrayList realIdxs = compressedResult.third;

        // line 1
        bw.write(n + "\t" + edgeCount + "\n");

        // line 2 ~ |V|+1
        //first中记录了每个点的所属的组号
        for(int i=0;i<n;i++){
            bw.write(realIdxs.getInt(i) + "\t" + compressedResult.first.getInt(i) + "\n");
        }
        //second中记录了邻接链表   有p cp cm 的邻接链表

        ObjectArrayList<Int2DoubleOpenHashMap>weightMap = module.returnWeight();

        /*while(it.hasNext()){
            int idx = it.nextIndex();
            AdjacencyList target = compressedResult.second[idx];
            it.next();
            for(int i=0;i<n;i++){
                for(int v: target.getNeighbors(i)){
                    if(i <= v){//防止输出两个 u，v   v,u 只输出一次
                        bw.write(i + "\t" + v +"\n");
                        summaryCount[idx]++;
                    }
                }
            }
            bw.write("-1\t-1\n");
        }*/
        // print out P and wright set
        AdjacencyList targetP = compressedResult.second[0];
        for(int i = 0;i < n;i++){
            for(int v:targetP.getNeighbors(i)){
                if(i <= v){
                    bw.write(i + "\t" + v + "\t" +":"+"\t");
                    ObjectArrayList<Double> weightSet = module.returnWeightSet(i,v);
                    for(double w:weightSet)
                        bw.write(w + "\t");
                    bw.write("\n");
                    summaryCount[0]++;
                }
            }
        }

        //print out Cp and weight
        bw.write("-1\t-1\n");

        AdjacencyList targetCp = compressedResult.second[1];
        for(int i=0;i<n;i++){
            for(int v: targetCp.getNeighbors(i)){
                if(i <= v){//防止输出两个 u，v   v,u 只输出一次
                    double w = weightMap.get(i).getOrDefault(v,-1);
                    bw.write(i + "\t" + v +"\t"+w+"\n");
                    summaryCount[1]++;
                }
            }
        }

        System.out.println("SUMMARY:");
        for(int i=0;i<2;i++){
            System.out.println("|" + listNames.get(i) + "| (Including Self-loop): " + summaryCount[i]);
        }
        /*int totalCount = (summaryCount[0] + summaryCount[1] + summaryCount[2]);
        System.out.println("total: " + totalCount);
        System.out.println("compression_rate **(Including Self-loop)**: " + totalCount / (double) edgeCount);
        System.out.println();*/
        bw.close();
    }
}
