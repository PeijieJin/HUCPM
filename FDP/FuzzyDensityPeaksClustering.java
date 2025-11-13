package FDP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FuzzyDensityPeaksClustering {
    private HashMap<String, Map<String, Double>> featureInstance; //特征及其实例与效用集<feature,<instance,utility>>
    private HashMap<String, Integer> featureCount; //特征的实例数量
    private double[][] fuzzyProximityMatrix;
    private double[][] fuzzyUtilityRatioMatrix;
    private HashMap<Integer, List<Double>> K_Nearest_Neighbor; //K近邻集
    private double pmax;
    private double gmax;

    public FuzzyDensityPeaksClustering(){
        featureInstance = new HashMap<>();
        featureCount = new HashMap<>();
        pmax = 0;
        gmax = Double.MAX_VALUE;
    }

    //初始化
    public void Initialization(String path1, String path2){
        readInstance(path1);
        System.out.println("每个特征及其实例总数为：" + featureCount);
        FeatureProximity f = new FeatureProximity();
        f.Initialization(path1, path2);
        fuzzyProximityMatrix = f.fuzzy_proximity_matrix(); //生成模糊邻近度矩阵
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

    //生成k近邻
    public void gen_K_Nearest_Neighbor(int k){
        K_Nearest_Neighbor = new HashMap<>();
        for (int i = 0; i < fuzzyUtilityRatioMatrix.length; i++){
            List<Double> fList = new ArrayList<>();
            for (int j = 0; j < fuzzyUtilityRatioMatrix.length; j++){
                fList.add(fuzzyUtilityRatioMatrix[i][j]);
            }
            //将第i个特征与其他特征的模糊邻近度保存到fList,将fList按升序排列
            Collections.sort(fList); //根据元素的自然顺序对指定fList集合的元素按升序进行排序
            List<Double> knn = new ArrayList<>();
            //将fList的k个最大值存入knn
            for (int j = fList.size() - 1; j > fList.size() - k - 1; j--){
                knn.add(fList.get(j));
            }
            K_Nearest_Neighbor.put(i, knn); //将每一个特征的k近邻保存到K_Nearest_Neighbor集中
        }
    }

    //计算每个特征的局部密度
    public HashMap<Integer, Double> computeLocalDensity(int k){
        HashMap<Integer, Double> localDensity = new HashMap<>(); //存储每个特征的局部密度
        gen_K_Nearest_Neighbor(k);
        for (int i = 0; i < fuzzyUtilityRatioMatrix.length; i++){
            double sump = 0.0, p=0.0;
            for (int j = 0; j < K_Nearest_Neighbor.get(i).size(); j++){
                sump += K_Nearest_Neighbor.get(i).get(j);
            }
            p = sump / k; //第i个特征的局部密度
            localDensity.put(i, p);
        }
        return localDensity;
    }

    //对输入的数据按降序排列
    public List<Integer> sort(HashMap<Integer, Double> listIn){
        List<Integer> listOut = new ArrayList<>(); //存储按照局部密度降序排列好的特征集合
        List<Map.Entry<Integer, Double>> list = new ArrayList<Map.Entry<Integer, Double>>(listIn.entrySet());
        Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {
            @Override
            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
                return o2.getValue().compareTo(o1.getValue()); //将list中的数据按照value值进行降序排列
            }
        });
        for (int i = 0; i < list.size(); i++){
            listOut.add(list.get(i).getKey());
        }
        System.out.println("按照局部密度降序排列好的特征集合:" + listOut);
        return listOut;
    }

    //计算delta值，生成特征及其delta值的集合gList.pfList是按局部密度降序排好的特征集合
    public HashMap<Integer, Double> computeRelativeDistance(List<Integer> pfList, HashMap<Integer, Double> localDensity){
        HashMap<Integer, Double> deltaList = new HashMap<>();
        pmax = localDensity.get(pfList.get(0)); //pList的第一个特征的局部密度就是最大局部密度
        double gmin = 0;
        //确定特征pList.get(0)的delta值
        for (int i = 0; i < pfList.size(); i++){
            //计算特征pList.get(0)与集合pList的其他特征之间最小特征邻近度
            if (fuzzyUtilityRatioMatrix[pfList.get(0)][pfList.get(i)] < gmin){
                gmin = fuzzyUtilityRatioMatrix[pfList.get(0)][pfList.get(i)];
            }
        }
        deltaList.put(pfList.get(0), 1 - gmin);

        //确定特征pList.get(1)到特征pList.get(pList.size()-1)的delta值
        for (int i = pfList.size() - 1; i > 0; i--){
            double max = 0.0;
            double min = Double.MAX_VALUE;
            //如果当前特征的局部密度等于pmax证明pList中所有特征的局部密度小于等于当前特征局部密度，满足delta求值公式的1-min
            if (localDensity.get(pfList.get(i)) == pmax){
                for (int j = 0; j < pfList.size(); j++){
                    if (fuzzyUtilityRatioMatrix[pfList.get(i)][pfList.get(j)] < min){
                        min = fuzzyUtilityRatioMatrix[pfList.get(i)][pfList.get(j)];
                    }
                }
                deltaList.put(pfList.get(i), 1 - min);
            }
            //如果当前特征的局部密度不等于pmax证明pList中必然存在其他特征的局部密度大于当前特征的局部密度，满足delta求值公式的1-max
            else {
                for (int j = i - 1; j > 0; j--){
                    if (fuzzyUtilityRatioMatrix[pfList.get(i)][pfList.get(j)] > max){
                        max = fuzzyUtilityRatioMatrix[pfList.get(i)][pfList.get(j)];
                    }
                }
                deltaList.put(pfList.get(i), 1 - max);
            }
        }
        return deltaList;
    }

    //计算Score值，选取聚类中心
    public HashMap<Integer, Double> computeScore(List<Integer> pfList, HashMap<Integer, Double> localDensity, HashMap<Integer, Double> deltaList, int n){
        double score = 0.0;
        HashMap<Integer, Double> scoreList = new HashMap<>();
        HashMap<Integer, Double> cList = new HashMap<>();
        List<Integer> gfList = sort(deltaList); //按delta值降序保存特征
        gmax = deltaList.get(gfList.get(0)); //取出最大delta值
        //计算第i个特征的Score
        for (int i = 0; i < pfList.size(); i++){
            score = (localDensity.get(pfList.get(i)) / pmax) * (deltaList.get(pfList.get(i)) / gmax);
            scoreList.put(pfList.get(i), score);
        }

        //将scoList排好序，选出其中score值最大的作为聚类中心
        List<Integer> scoreList1 = sort(scoreList);
        //取出前n个聚类中心及Score值
        for (int i = 0; i < n; i++){
            cList.put(scoreList1.get(i), scoreList.get(scoreList1.get(i)));
        }
        System.out.println("聚类中心及Sorce值为：" + cList);
        return cList;
    }

    //计算每个特征对聚类中心的隶属度
    public double[][] computeMembershipofFeature(HashMap<Integer, Double> cList){
        int length = fuzzyUtilityRatioMatrix.length;
        double[][] membershipMatrix = new double[length][length];
        Set<Integer> cfList = cList.keySet(); //cfList存储作为聚类中心的特征
        for (int i : cfList){
            for (int j = 0; j < length; j++){
                if (fuzzyUtilityRatioMatrix[j][i] == 1){
                    membershipMatrix[j][i] = 1;
                } else if (fuzzyUtilityRatioMatrix[j][i] == 0){
                    membershipMatrix[j][i] = 0;
                }else {
                    double sum = 0.0;
                    for (int k : cfList){
                        if(fuzzyUtilityRatioMatrix[j][k] != 0) {
                            sum += Math.pow(fuzzyUtilityRatioMatrix[j][i] / fuzzyUtilityRatioMatrix[j][k], 2);
                        }
                    }
                    membershipMatrix[j][i] = sum;
                }
            }
        }
        return membershipMatrix;
    }

    //将特征对聚类中心的隶属度计算出来，其中大于最小隶属度阈值min_u的特征加入该聚类中心的簇中
    public HashMap<Integer, Set<Integer>> gen_cluster_result(HashMap<Integer, Double> cList, double min_u){
        HashMap<Integer, Set<Integer>> clusterResult = new HashMap<>(); //存储聚类结果
        double[][] membershipMatrix = computeMembershipofFeature(cList); //生成特征对簇的隶属度矩阵，一列则是每个特征对聚类中心的隶属度，即行代表所有特征，列代表所有聚类中心
        Set<Integer> cfList = cList.keySet(); //存储聚类中心代表的特征
        for (Integer c : cfList){
            HashSet<Integer> features = new HashSet<>(); //存储隶属度大于min_u的特征集
            for (int i = 0; i < fuzzyUtilityRatioMatrix.length; i++){
                if (i != c && membershipMatrix[i][c] >= min_u){
                    features.add(i);
                }
            }
            clusterResult.put(c, features);
        }

        //将聚类结果打印出来
        System.out.println("聚类结果为（数字表示）：");
        for (Integer c : cfList){
            System.out.println(c + "\t" + clusterResult.get(c));
        }
        return clusterResult;
    }

    //将聚类结果有Integer类型转换为String类型
    public void clusterResfromInttoStr(HashMap<Integer, Set<Integer>> clusterResult){
        HashMap<String, Set<String>> clusterResofStr = new HashMap<>();
        List<String> fList = new ArrayList<>();
        fList.addAll(featureCount.keySet()); //生成特征集
        System.out.println("特征集为：" + fList);
        for (Map.Entry<Integer, Set<Integer>> entry : clusterResult.entrySet()){
            HashSet<String> resSet = new HashSet<>();
            for (int i : entry.getValue()){
                resSet.add(fList.get(i)); //从fList取出第i个特征的字母表示放进resSet中
            }
            clusterResofStr.put(fList.get(entry.getKey()), resSet); //fList.get(entry.getKey())表示entry.getKey()所代表的数字在fList中的字母表示
        }
        System.out.println("聚类结果为（字母表示）：");
        for (Map.Entry<String, Set<String>> entry : clusterResofStr.entrySet()){
            System.out.println(entry.getKey() + "\t" + entry.getValue());
        }
        writeResult(clusterResofStr, "result/jOutData/out_FDP_(1000-1500).txt"); //将结果写入文件中
    }

    //将聚类结果写入文件中
    private void writeResult(HashMap<String, Set<String>> clusterResofStr, String path) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(path)));
            writer.write("聚类中心\t簇");
            writer.newLine();
            for (Map.Entry<String, Set<String>> entry : clusterResofStr.entrySet()){
                writer.write(entry.getKey() + "\t" + entry.getValue());
                writer.newLine();
            }
            writer.flush();
            writer.close();
        } catch (IOException e){
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        FuzzyDensityPeaksClustering fdp = new FuzzyDensityPeaksClustering();
        fdp.Initialization("data-45.txt", "fuzzyneighbors1(300-350).txt");
//        fdp.Initialization("poi_NYC.txt", "fuzzyneighbors_NYC_(10-15).txt");
        HashMap<Integer, Double> localDensity = fdp.computeLocalDensity(10); //计算每个特征的局部密度
        List<Integer> pfList = fdp.sort(localDensity); //按局部密度降序将特征排序
        HashMap<Integer, Double> deltaList = fdp.computeRelativeDistance(pfList, localDensity); //计算每个特征的delta值
        HashMap<Integer, Double> cList = fdp.computeScore(pfList, localDensity, deltaList, 6); //计算聚类中心
        HashMap<Integer, Set<Integer>> clusterResult = fdp.gen_cluster_result(cList, 0.1); //生成有数字表示的聚类结果
        fdp.clusterResfromInttoStr(clusterResult);
    }
}