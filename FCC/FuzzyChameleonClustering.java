package FC;

import java.io.*;
import java.util.*;

public class FuzzyChameleonClustering {
    private double[][] fuzzyUtilityRatioMatrix;
    private HashMap<Integer, List<FeatureFUF>> k_Nearest_Neighbor;
    private HashMap<Integer, List<FeatureFUF>> initialization_Subcluster;
    private HashMap<String, Map<String, Double>> featureInstance; //特征及其实例与效用集
    private HashMap<String, Integer> featureCount; //特征的实例数量
    private List<Integer> list;
    private List<List<Integer>> res;

    public FuzzyChameleonClustering(){
        featureInstance = new HashMap<>();
        featureCount = new HashMap<>();
    }

    //初始化，先生成fuzzyUtilityRatioMatrix矩阵
    public void Initialization(String path1, String path2){
        readInstance(path1);
        FuzzyUtilityRatio fur = new FuzzyUtilityRatio();
        fur.Initialization(path1, path2);
        fuzzyUtilityRatioMatrix = fur.fuzzyUtilityRatioMatrix();
    }
    //读取文件中的特征的实例数量数量和所有特征及其实例与效用的集合
    //path中的文件格式为“id f uti(f.id) x坐标 y坐标 z坐标”
    private void readInstance(String path){
        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(path)));
            String line = null;

            while ((line = reader.readLine()) != null){
                if (line.equals("")){
                    continue;
                }

                String[] tmp = line.split("\t");
                String feature = tmp[1];
                String id = tmp[0];
                String instance = tmp[1] + "." + tmp[0];
                double utility = Double.parseDouble(tmp[2]);

                if (!featureInstance.containsKey(feature)){
                    featureInstance.put(feature, new HashMap<>());
                }
                featureInstance.get(feature).put(instance, utility);

                if (!featureCount.containsKey(feature)){
                    featureCount.put(feature, 0);
                }
                featureCount.put(feature, featureCount.get(feature) + 1); //HashMap当key值相同时，前一个value值会被后一个value值覆盖
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    //构建基于特征间的模糊效用率的k邻近图
    public void gen_fuzzy_k_nearest_neighbor(int k){
        k_Nearest_Neighbor = new HashMap<>();
        for (int i = 0; i < fuzzyUtilityRatioMatrix.length; i++){
            List<FeatureFUF> furList = new ArrayList<>();
            for (int j = 0; j < fuzzyUtilityRatioMatrix.length; j++){
                furList.add(new FeatureFUF(j, fuzzyUtilityRatioMatrix[i][j]));
            }

            //根据指定比较器诱导的顺序对指定列表进行排序
            //以从大到小的顺序排列
            Collections.sort(furList, new Comparator<FeatureFUF>() {
                @Override
                public int compare(FeatureFUF o1, FeatureFUF o2) {
                    if (o1.getFuf() < o2.getFuf()){
                        return 1; //交换位置
                    }
                    if (o1.getFuf() == o2.getFuf()){
                        return 0; //不交换位置
                    }
                    return -1; //不交换位置
                }
            });
            List<FeatureFUF> knn = new ArrayList<>();
            for (int j = 0; j < k; j++){
                knn.add(furList.get(j)); //取出前k个特征及其模糊效用率
            }
            k_Nearest_Neighbor.put(i, knn);
        }
    }

    //初始化子簇
    public void gen_initialization_subcluster(double min_u, int m){
        initialization_Subcluster = new HashMap<>();
        for (int i = 0; i < k_Nearest_Neighbor.size(); i++){
            int c = i;
            List<FeatureFUF> fList = new ArrayList<>();
            fList.addAll(k_Nearest_Neighbor.get(c)); //记录该聚类中心对应的全部k近邻
            for (int j = 0; j < fList.size(); j++){

                int f = fList.get(j).getFeature(); //fuf第j大的特征
                double fuf = fList.get(j).getFuf();

                double u = computeMembershipofCluster(f, fuf, m); //计算特征对簇的隶属度

                if (u < min_u){
                    fList.remove(new FeatureFUF(f, fuf));
                }
            }
            initialization_Subcluster.put(i, fList);
        }
    }

    //计算特征对簇的隶属度
    private double computeMembershipofCluster(int f, double fuf, int m) {
        double u = 0.0;
        for (int k = 0; k < k_Nearest_Neighbor.size(); k++){
            List<FeatureFUF> fList = new ArrayList<>();
            fList.addAll(k_Nearest_Neighbor.get(k));
            for (int j = 0; j < fList.size(); j++){
                if (fList.get(j).getFeature() == f){
                    u += Math.pow(fuf / fList.get(j).getFuf(), 2 / (m - 1));
                    break;
                }
            }
        }
        return u;
    }


    //合并子簇
    public HashMap<Set<Integer>, Double> gen_merging_subcluster(double min_fuf, double min_sim, int a){
        boolean flag = true;
        List<List<FeatureFUF>> CL = new ArrayList<>();
        HashMap<Set<Integer>, Double> FUP = new HashMap<>(); //存储Fuzzy Utility Patterns
        HashMap<Set<Integer>, Double> FHUP = new HashMap<>(); //存储模糊高效用模式
        double FUoC = 0.0; //模式的模糊模式效用度
        double AUoC = 0.0; //模式特征簇的平均效用度
        for (int i = 0; i < initialization_Subcluster.size(); i++){
            CL.add(initialization_Subcluster.get(i));
        }

        //初始化边的集合，EdgeFUF保存边和边上的权重，即特征间的模糊效用率
        List<List<EdgeFUF>> EdgeSet = new ArrayList<>();
        initializationEdgeSet(EdgeSet);
//        System.out.println(EdgeSet);

        while (flag){
            flag = false;
            double maxSim = 0.0;
            String simid = null;
            List<EdgeFUF> newEdge = new ArrayList<>();
            for (int i = 0; i < CL.size(); i++){
                List<FeatureFUF> CL1 = new ArrayList<>();
                CL1.addAll(CL.get(i)); //第i个簇
                for (int j = 0; j < CL.size(); j++){
                    double sim = 0.0;
                    if (i == j){
                        continue;
                    }
                    List<FeatureFUF> CL2 = new ArrayList<>();
                    CL2.addAll(CL.get(j)); //第j个簇
                    Map<String, Double> NFS = gen_neighbor_feature_pairs(CL1, CL2, min_fuf);
//                    System.out.println(NFS);
                    if (NFS != null){
                        sim = computeSimilarity(CL1, CL2, a, NFS);
//                        System.out.println(sim);
                    }
                    if (sim > min_sim && sim > maxSim){
                        System.out.println(sim);
                        flag = true;
                        maxSim = sim;
                        simid = i + "," + j;
//                        System.out.println(simid);

                        newEdge.clear();
                        for(Map.Entry<String, Double> entry : NFS.entrySet()) {
                            EdgeFUF edgefuf = new EdgeFUF(entry.getKey(), entry.getValue());
                            newEdge.add(edgefuf);
                        }
                    }
                }
            }
            if (simid != null){
                int index = simid.lastIndexOf(",");
                int id1 = Integer.parseInt(simid.substring(0, index));
                int id2 = Integer.parseInt(simid.substring(index + 1));
//                System.out.println(id1 + "," + id2);
                List<FeatureFUF> CL1 = CL.get(id1); //取出id1代表的簇
                List<FeatureFUF> CL2 = CL.get(id2); //取出id2代表的簇
                List<FeatureFUF> newCL = mergingCL(CL1, CL2);
                CL.remove(CL1);
                CL.remove(CL2);
                CL.add(newCL);
                List<EdgeFUF> EdgeList1 = EdgeSet.get(id1); //取出id1代表的边的集合
                List<EdgeFUF> EdgeList2 = EdgeSet.get(id2); //取出id2代表的边的集合
                List<EdgeFUF> newEdgeList = mergingEdge(EdgeList1, EdgeList2, newEdge);
                EdgeSet.remove(EdgeList1);
                EdgeSet.remove(EdgeList2);
                EdgeSet.add(newEdgeList);
            }
        }

        //生成模式特征簇及计算每个模式的模糊模式效用度
        for(int i = 0; i < EdgeSet.size(); i++) {
            List<EdgeFUF> list = EdgeSet.get(i);
            FUoC = computeFUoC(list);
            Set<Integer> FUC = gen_fuzzy_utility_patterns(list);
            FUP.put(FUC, FUoC);
        }

        //计算模式特征簇的平均效用度，并生成模糊高效用模式
        double sumAUoC = 0.0;
        for(Map.Entry<Set<Integer>, Double> entry : FUP.entrySet()) {
            sumAUoC += entry.getValue();
        }
        AUoC = sumAUoC / FUP.size();
        for(Map.Entry<Set<Integer>, Double> entry : FUP.entrySet()) {
            if(entry.getValue() >= AUoC) {
                FHUP.put(entry.getKey(), entry.getValue());
            }
        }

        return FHUP;
    }

    private Set<Integer> gen_fuzzy_utility_patterns(List<EdgeFUF> list) {
        Set<Integer> FUC = new HashSet<>();
        for(EdgeFUF edgefuf : list) {
            String edge = edgefuf.getEdge();
            String[] tmp = edge.split(",");
            int feature1 = Integer.parseInt(tmp[0]);
            FUC.add(feature1);
            int feature2 = Integer.parseInt(tmp[1]);
            FUC.add(feature2);
        }
        return FUC;
    }

    private double computeFUoC(List<EdgeFUF> list) {
        double sumFUoC = 0.0;
        for(EdgeFUF edgefuf : list) {
            sumFUoC += edgefuf.getFuf();
        }
        double FUoC = sumFUoC / list.size();
        return FUoC;
    }

    private List<EdgeFUF> mergingEdge(List<EdgeFUF> EdgeList1, List<EdgeFUF> EdgeList2, List<EdgeFUF> newEdge) {
        List<EdgeFUF> newEdgeList = new ArrayList<>();
        newEdgeList.addAll(EdgeList1);

        for (EdgeFUF edgefuf : EdgeList2) {
            boolean flag = true;
            String edge1 = edgefuf.getEdge();
            String[] tmp = edge1.split(",");
            String edge2 = tmp[1] + "," + tmp[0];
            for (int i = 0; i < EdgeList1.size(); i++) {
                String edge3 = EdgeList1.get(i).getEdge();
                if(edge1.equals(edge3) || edge2.equals(edge3)) {
                    flag = false; //EdgeList1中已经存在该边
                    break;
                }
            }
            if(flag) {
                newEdgeList.add(edgefuf);
            }
        }

        List<EdgeFUF> EdgeList3 = new ArrayList<>();
        EdgeList3.addAll(newEdgeList);
        for(EdgeFUF edgefuf : newEdge) {
            boolean flag = true;
            String edge1 = edgefuf.getEdge();
            String[] tmp = edge1.split(",");
            String edge2 = tmp[1] + "," + tmp[0];
            for (int i = 0; i < EdgeList3.size(); i++) {
                String edge3 = EdgeList3.get(i).getEdge();
                if(edge1.equals(edge3) || edge2.equals(edge3)) {
                    flag = false; //EdgeList1中已经存在该边
                    break;
                }
            }
            if(flag) {
                newEdgeList.add(edgefuf);
            }
        }

        return newEdgeList;
    }

    private void initializationEdgeSet(List<List<EdgeFUF>> EdgeSet) {
        for(int i = 0; i < initialization_Subcluster.size(); i++) {
            int f1 = i; //初始化子簇的中心特征
            List<FeatureFUF> list1 = initialization_Subcluster.get(i);
            List<EdgeFUF> list2 = new ArrayList<>();
            for (int j = 0; j < list1.size(); j++) {
                int f2 = list1.get(j).getFeature();
                double FUF = list1.get(j).getFuf();
                String edge = f1 + "," + f2;
                EdgeFUF edgefuf = new EdgeFUF(edge, FUF);
                list2.add(edgefuf);
            }
            EdgeSet.add(list2);
        }
    }

    private void initializationCRtemp(List<List<Integer>> CRtemp) {
        for (int i = 0; i < initialization_Subcluster.size(); i++){
            List<Integer> s = new ArrayList<>();
            int cf = i; //中心特征
            s.add(cf);
            List<FeatureFUF> list = initialization_Subcluster.get(i);
            for (int j = 0; j < list.size(); j++){
                int f = list.get(j).getFeature();
                s.add(f); //簇中除了中心特征的其他特征
            }
            CRtemp.add(s); //s中是包含中心特征及其初始化后knn中的所有特征
        }

    }

    private List<FeatureFUF> mergingCL(List<FeatureFUF> CL1, List<FeatureFUF> CL2) {
        List<FeatureFUF> newCL = new ArrayList<>();
        newCL.addAll(CL1);
        newCL.addAll(CL2);
        return newCL;
    }

    private double computeSimilarity(List<FeatureFUF> CL1, List<FeatureFUF> CL2, int a, Map<String, Double> NFS) {
        double ri = computeRelativeInterconnectivity(CL1, CL2, NFS);
//        System.out.println(ri);
        double rc = computeRelativeClosenss(CL1, CL2, NFS);
//        System.out.println(rc);
        double sim = Math.pow(ri * rc, a);
        return sim;
    }

    private double computeRelativeClosenss(List<FeatureFUF> CL1, List<FeatureFUF> CL2,  Map<String, Double> NFS) {
        double rc = 0.0;
        double sumsfe = 0.0, sfe = 0.0, sumse1 = 0.0, sumse2 = 0.0;
        List<String> fmember = new ArrayList<>();
        fmember.addAll(NFS.keySet());
        for (int i = 0; i < NFS.size(); i++){
            sumsfe += NFS.get(fmember.get(i));
        }
//        System.out.println(sumsfe);
        sfe = sumsfe / NFS.size();
        for (int i = 0; i < CL1.size(); i++){
            sumse1 += CL1.get(i).getFuf();
        }
//        System.out.println(sumse1);
        double se1 = sumse1 / CL1.size();
        for (int i = 0; i < CL2.size(); i++){
            sumse2 += CL2.get(i).getFuf();
        }
//        System.out.println(sumse2);
        double se2 = sumse2 / CL2.size();
        double part1 = CL1.size() * se1 / (CL1.size() + CL2.size());
//        System.out.println(part1);
        double part2 = CL2.size() * se2 / (CL1.size() + CL2.size()) ;
//        System.out.println(part2);
        rc = sfe / (part1 + part2);
        return rc;
    }

    private double computeRelativeInterconnectivity(List<FeatureFUF> CL1, List<FeatureFUF> CL2,  Map<String, Double> NFS) {
        double ri = 0.0;
        double ef = 0.0, sumse1 = 0.0, sumse2 = 0.0;
        List<String> fmember = new ArrayList<>();
        fmember.addAll(NFS.keySet());
        for (int i = 0; i < NFS.size(); i++){
            ef += NFS.get(fmember.get(i));
        }
//        System.out.println(ef);
        for (int i = 0; i < CL1.size(); i++){
            sumse1 += CL1.get(i).getFuf();
        }
//        System.out.println(sumse1);
        double se1 = sumse1 / CL1.size();
        for (int i = 0; i < CL2.size(); i++){
            sumse2 += CL2.get(i).getFuf();
        }
//        System.out.println(sumse2);
        double se2 = sumse2 / CL2.size();
        ri = ef / ((se1 + se2) * NFS.size());
        return ri;
    }

    //生成邻近特征对集合
    private Map<String, Double> gen_neighbor_feature_pairs(List<FeatureFUF> CL1, List<FeatureFUF> CL2, double min_fuf) {
        List<String> fList = new ArrayList<>();
        fList.addAll(featureCount.keySet());
        Map<String, Double> eSet = new HashMap<>();
        for (int i = 0; i < CL1.size(); i++){
            int f1 = CL1.get(i).getFeature(); //特征的编号为f1
            for (int j = 0; j < CL2.size(); j++){
                int f2 = CL2.get(j).getFeature(); //特征的编号为f2
                if(f1 != f2) {
                    double fuf = fuzzyUtilityRatioMatrix[f1][f2];
                    String edge = f1 + "," + f2;
                    if (! eSet.containsKey(f1 + "," + f2) && !eSet.containsKey(f2 + "," + f1) && fuf > min_fuf){
                        eSet.put(edge, fuf);
                    }
                }
            }
        }
        return eSet;
    }

    //将聚类结果由Integer类型转化为String类型，并将结果写入到文件中
    public void integertoString(HashMap<Set<Integer>, Double> FHUP){
        Map<Set<String>, Double> clusterResult = new HashMap<>();
        List<String> fList = new ArrayList<>();
        fList.addAll(featureCount.keySet()); //所有特征的集合
        System.out.println("特征集为：" + fList);
        for(Map.Entry<Set<Integer>, Double> clusterResultofInt : FHUP.entrySet()) {
            Set<String> resSet = new HashSet<>();
            Set<Integer> cluster = clusterResultofInt.getKey();
            for(int i : cluster) {
                resSet.add(fList.get(i));
            }
            clusterResult.put(resSet, clusterResultofInt.getValue());
        }
        System.out.println("聚类结果为（字母表示）：");
        Set<Set<String>> clusters = clusterResult.keySet();
        for(Set<String> cluster : clusters) {
            System.out.println(cluster);
        }
        writeResult(clusterResult, "result/out_FCH_data-20-8000(utility).txt"); //将结果写入文件中
    }

    private void writeResult(Map<Set<String>, Double> clusterResult, String path) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(path)));
            writer.write("特征簇\t模糊模式效用度");
            writer.newLine();
            for (Map.Entry<Set<String>, Double> entry : clusterResult.entrySet()){
                writer.write(entry.getKey() + "\t" + entry.getValue());
                writer.newLine();
            }
            writer.flush();
            writer.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    //生成模糊高效用模式  //FHUC为模糊高效用簇,min_fui为模糊模式效用度阈值
    public Map<Set<Integer>, Double> gen_fuzzy_high_utility_patterns(HashMap<Set<Integer>, Double> FHUC, double min_fui){
        Map<Set<Integer>, Double> FHUP = new HashMap<>(); //模糊高效用模式
        for (Map.Entry<Set<Integer>, Double> entry1 : FHUC.entrySet()){
            Set<Set<Integer>> candidate = gen_candidate_pattern(entry1.getKey());
            Map<Set<Integer>, Double> candidateUtility = computeCandidateFuzzyUtility(candidate);
            for (Map.Entry<Set<Integer>, Double> entry2 : candidateUtility.entrySet()){
                if (entry2.getValue() >= min_fui){
                    FHUP.put(entry2.getKey(), entry2.getValue());
                }
            }
        }

        return FHUP;
    }

    //生成候选模式集
    private Set<Set<Integer>> gen_candidate_pattern(Set<Integer> item) {
        List<Integer> listtmp  = new ArrayList<>();
        listtmp.addAll(item);
        list = new ArrayList<>();
        res = new ArrayList<>();
        backtracking(listtmp, 0);

        Set<Set<Integer>> set = new HashSet<>();

        for (List<Integer> re : res) {
            set.add(new HashSet<>(re));
        }
        return set;
    }

    private void backtracking(List<Integer> listtmp, int startIndex) {
        if(list.size()>=2){
            res.add(new ArrayList<>(list));
        }
        for(int i = startIndex ; i < listtmp.size(); i++){
            list.add(listtmp.get(i));
            backtracking(listtmp, i+1);
            list.remove(list.size()-1);
        }
    }

    //计算每个候选模式集的效用
    private Map<Set<Integer>, Double> computeCandidateFuzzyUtility(Set<Set<Integer>> candidate) {
        Map<Set<Integer>, Double> map = new HashMap<>();
        double sumfuf = 0.0;
        for (Set<Integer> set : candidate){
            for (int i : set){
                for (int j : set){
                    if (i != j){
                        sumfuf += fuzzyUtilityRatioMatrix[i][j];
                    }
                }
            }
            double fuzzyUtilityIndex = sumfuf / Math.pow(set.size(), 2);
            map.put(set, fuzzyUtilityIndex);
        }
        return map;
    }

    private double computeAvgFUI(HashMap<Set<Integer>, Double> FHUC){
        double sumfui = 0.0;
        int countsize = 0;
        for (Map.Entry<Set<Integer>, Double> entry1 : FHUC.entrySet()){
            Set<Set<Integer>> candidate = gen_candidate_pattern(entry1.getKey());
            Map<Set<Integer>, Double> candidateUtility = computeCandidateFuzzyUtility(candidate);
            countsize += candidate.size();
            for (Map.Entry<Set<Integer>, Double> entry2 : candidateUtility.entrySet()){
                sumfui += entry2.getValue();
            }
        }
        double avg_fui = sumfui /countsize;
        return avg_fui;
    }

    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        FuzzyChameleonClustering fcc = new FuzzyChameleonClustering();
        fcc.Initialization("gen/data-20-8000(utility).txt", "gen/fuzzyneighbor(50)_data-20-8000(utility).txt");
        fcc.gen_fuzzy_k_nearest_neighbor(6); //生成k近邻
        fcc.gen_initialization_subcluster(0.3, 2); //初始化子簇
        HashMap<Set<Integer>, Double> FHUC =  fcc.gen_merging_subcluster(0.7, 0.5, 1); //合并子簇
        fcc.integertoString(FHUC);
        long end = System.currentTimeMillis();
        System.out.println("runtime:"+(end-start)/1000.0+"s");
    }
}