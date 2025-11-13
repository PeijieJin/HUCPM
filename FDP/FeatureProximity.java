package FDP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FeatureProximity {
    private HashMap<String, List<String>> featureInstance; //所有特征及其实例集
    private HashMap<String, Integer> featureCount; //所有特征的实例数量集
    private HashMap<String, Map<String, Map<String, Double>>> starNeighbor; //星型实例集
    private HashMap<String, Double> fuzzymembership;//存储模糊邻近关系及其隶属度
    private HashMap<String, Map<String, Map<String, Double>>> fuzzyStarNeighbor; //模糊星型实例集
    private HashMap<String, Map<String, Double>> fuzzyStarNeighborAll; //模糊星型实例集及隶属度之和

    public FeatureProximity(){
        featureInstance = new HashMap<>();
        featureCount = new HashMap<>();
        starNeighbor = new HashMap<>();
        fuzzymembership = new HashMap<>();
        fuzzyStarNeighbor = new HashMap<>();
        fuzzyStarNeighborAll = new HashMap<>();
    }

    //读入数据，featureInstances所有特征及其实例集，featureCount所有特征的实例数量集
    //通过path读入的文件为"实例 特征 x坐标 y坐标 z坐标"
    public void readInstance(String path){
        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(path)));
            String line;

            while ((line = reader.readLine()) != null){
                if (line.equals("")){
                    continue;
                }

                String tmp[] = line.split("\t");
                String feature = tmp[1];
                String id = tmp[0];
                String ins = feature + "." + id;

                if (!featureInstance.containsKey(feature)){
                    List<String> insList = new ArrayList<>();
                    insList.add(ins);
                    featureInstance.put(feature, insList);
//                	featureInstance.put(feature, new ArrayList<>());
                }
                else {
                    featureInstance.get(feature).add(ins); //一个特征的实例只会在数据集中出现一次，故可以直接加
                }

                if (!featureCount.containsKey(feature)){
                    featureCount.put(feature, 1);
                } else {
                    featureCount.put(feature, featureCount.get(feature) + 1); //HashMap当key值相同时，前一个value值会被后一个value值覆盖
                }
            }
            reader.close();
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    //生成starNeighbor(特征1, (实例1, (特征2.实例2, 模糊隶属度)))
    //生成fuzzyStarNeighbor(特征1.实例1, (特征2, (特征2.实例2, 模糊隶属度)))
    //生成fuzzyStarNeighborAll(特征1.实例1, (特征2, 实例1与特征2的所有隶属度>0的实例的隶属度之和))
    //通过path读入的文件为计算后的"特征1.实例1,特征2.实例2,模糊邻近隶属度"
    public void gen_fuzzy_starNeighbor(String path) throws FileNotFoundException, Exception{
        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(path)));
            String line;

            while ((line = reader.readLine()) != null){
                if (line.equals("")){
                    continue;
                }

                String[] tmp = line.split(",");
                String ins1 = tmp[0];
                int index1 = ins1.indexOf(".");
                String f1 = ins1.substring(0, index1);
                String id1 = ins1.substring(index1 + 1);
                String ins2 = tmp[1];
                int index2 = ins2.indexOf(".");
                String f2 = ins2.substring(0, index2);
                String id2 = ins2.substring(index2 + 1);
                double u = Double.parseDouble(tmp[2]);

                //生成星型邻居集
                if (!starNeighbor.containsKey(f1)){
                    starNeighbor.put(f1, new HashMap<>());
                }
                if (!starNeighbor.get(f1).containsKey(ins1)){
                    starNeighbor.get(f1).put(ins1, new HashMap<>());
                }
                starNeighbor.get(f1).get(ins1).put(ins2, u);
                if (!starNeighbor.containsKey(f2)){
                    starNeighbor.put(f2, new HashMap<>());
                }
                if (!starNeighbor.get(f2).containsKey(ins2)){
                    starNeighbor.get(f2).put(ins2, new HashMap<>());
                }
                starNeighbor.get(f2).get(ins2).put(ins1, u);

                //fuzzymembership存储模糊邻近关系及其隶属度
                String member = ins1 + "," + ins2;
                if (!fuzzymembership.containsKey(member)){
                    fuzzymembership.put(member, u);
                }

                //生成模糊星型实例集fuzzyStarNeighbor、模糊星型实例集及隶属度之和fuzzyStarNeighborAll
                //存储<ins1, <f2, <ins2, u>>>
                //保存<ins2, u>
                HashMap<String, Double> ins2u = new HashMap<>();
                ins2u.put(ins2, u);
                //保存<f2, <ins2, u>>
                HashMap<String, Map<String, Double>> f2ins2u = new HashMap<>();
                f2ins2u.put(f2, ins2u);
                if (fuzzyStarNeighbor.containsKey(ins1)){
                    if (fuzzyStarNeighbor.get(ins1).containsKey(f2)){
                        if (!fuzzyStarNeighbor.get(ins1).get(f2).containsKey(ins2)){
                            fuzzyStarNeighbor.get(ins1).get(f2).put(ins2, u);

                            Double preu = fuzzyStarNeighborAll.get(ins1).get(f2);
                            fuzzyStarNeighborAll.get(ins1).put(f2, preu + u);
                        }
                    } else {
                        fuzzyStarNeighbor.get(ins1).put(f2, ins2u);

                        fuzzyStarNeighborAll.get(ins1).put(f2, u);
                    }
                } else {
                    fuzzyStarNeighbor.put(ins1, f2ins2u);

                    //存储<ins1, <f2, u>>
                    fuzzyStarNeighborAll.put(ins1, new HashMap<String, Double>());
                    fuzzyStarNeighborAll.get(ins1).put(f2, u);
                }

                //存储<ins2, <f1, <ins1, u>>>
                //保存<ins1, u>
                HashMap<String, Double> ins1u = new HashMap<>();
                ins1u.put(ins1, u);
                //保存<f1, <ins1, u>>
                HashMap<String, Map<String, Double>> f1ins1u = new HashMap<>();
                f1ins1u.put(f1, ins1u);
                if (fuzzyStarNeighbor.containsKey(ins2)){
                    if (fuzzyStarNeighbor.get(ins2).containsKey(f1)){
                        if (!fuzzyStarNeighbor.get(ins2).get(f1).containsKey(ins1)){
                            fuzzyStarNeighbor.get(ins2).get(f1).put(ins1, u);

                            Double preu = fuzzyStarNeighborAll.get(ins2).get(f1);
                            fuzzyStarNeighborAll.get(ins2).put(f1, preu + u);
                        }
                    } else {
                        fuzzyStarNeighbor.get(ins2).put(f1, ins1u);

                        fuzzyStarNeighborAll.get(ins2).put(f1, u);
                    }
                } else {
                    fuzzyStarNeighbor.put(ins2, f1ins1u);

                    //存储<ins2, <f1, u>>
                    fuzzyStarNeighborAll.put(ins2, new HashMap<String, Double>());
                    fuzzyStarNeighborAll.get(ins2).put(f1, u);
                }
            }
            reader.close();
        } catch (Throwable e){
            e.printStackTrace();
        }
    }

    //将模糊星型实例集写入文件中
    public void writeFS(HashMap<String, Map<String, Map<String, Double>>> fuzzyStarNeighbor){
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File("midj/starNeighbor_FDP.txt")));
            writer.write("中心实例\t邻近特征\t邻近特征的实例及实例间的隶属度");
            writer.newLine();
            for (Map.Entry<String, Map<String, Map<String, Double>>> entry1 : fuzzyStarNeighbor.entrySet()){
                String f = entry1.getKey();
                for (Map.Entry<String, Map<String, Double>> entry2 : entry1.getValue().entrySet()){
                    writer.write(f + "\t" + entry2.getKey() + "\t" + entry2.getValue());
                    writer.newLine();
                }
            }
            writer.flush();
            writer.close();
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    //生成模糊邻近度矩阵
    private double[][] gen_fuzzy_proximity_matrix(){
        int size = featureCount.size();
        double[][] fuzzyProximityMatrix = new double[size][size];
        List<String> fList = new ArrayList<>();
        fList.addAll(featureCount.keySet()); //keySet()只获取key值
        Collections.sort(fList); //根据元素的自然顺序对指定 fList 集合的元素按升序进行排序。

        //计算右上三角模糊邻近度
        for (int i = 0; i < fList.size(); i++){
            String f1 = fList.get(i);
            for (int j = i; j < fList.size(); j++){
                String f2 = fList.get(j);
                fuzzyProximityMatrix[i][j] = computeFuzzyProximity(f1, f2);
            }
        }
        //将模糊邻近度对称到左下三角
        for (int i = 1; i <fList.size(); i++){
            for (int j = 0; j < i; j++){
                fuzzyProximityMatrix[i][j] = fuzzyProximityMatrix[j][i];
            }
        }
//        System.out.println("特征邻近度矩阵为：");
//        printMatrix(fuzzyProximityMatrix);
        return fuzzyProximityMatrix;
    }

    //计算模糊邻近度
    private double computeFuzzyProximity(String f1, String f2) {
        double fp = 0;
        fp = (computeFuzzyProximityRatio(f1, f2) + computeFuzzyProximityRatio(f2, f1)) / 2;
        return fp;
    }

    //计算模糊邻近率
    private double computeFuzzyProximityRatio(String f1, String f2) {
        double fpr = 0.0;
        List<String> insList = new ArrayList<>();
        if(starNeighbor.containsKey(f1)) {
            insList.addAll(starNeighbor.get(f1.toString()).keySet()); //取得特征f1所有有邻近关系的f1的实例
        }
        String ins;
        double maxu = 0.0, sum = 0.0;
        for (int i = 0; i < insList.size(); i++){
            ins = insList.get(i);
            maxu = computeCumulativeShareSorce(ins, f2); //计算f1的每一个实例和f2的共享分数
            sum += maxu;
        }
        fpr = sum / featureCount.get(f1);
        return fpr;
    }

    //计算实例与特征之间的累计共享分数
    private double computeCumulativeShareSorce(String instance, String feature) {
        double cs = 0.0, sum = 0.0;
        if (fuzzyStarNeighbor.containsKey(instance)){
            if (fuzzyStarNeighbor.get(instance).containsKey(feature.toString())){
                List<String> ins = new ArrayList<>();
                ins.addAll(fuzzyStarNeighbor.get(instance).get(feature.toString()).keySet()); //keySet():返回映射中所有 key 组成的 Set 视图
                for (int i = 0; i < ins.size(); i++){
                    sum += computeShareSorce(instance, ins.get(i));
                }
            }
        }

        if (sum > 1){
            cs = 1;
        } else {
            cs = sum;
        }
        return cs;
    }

    private double computeShareSorce(String ins1, String ins2) {
        double ss = 0.0;
        int index1 = ins1.indexOf(".");
        String f1 = ins1.substring(0, index1);
        String id1 = ins1.substring(index1 + 1);
        int index2 = ins2.indexOf(".");
        String f2 = ins2.substring(0, index2);
        String id2 = ins2.substring(index2 + 1);

        String s1 = ins1 + "," + ins2;
        String s2 = ins2 + "," + ins1;

        if (fuzzymembership.containsKey(s1)){
            double sumu = fuzzyStarNeighborAll.get(ins2).get(f1);
            double u = fuzzymembership.get(s1);
            ss = u / sumu;
        } else if (fuzzymembership.containsKey(s2)){
            double sumu = fuzzyStarNeighborAll.get(ins2).get(f1);
            double u = fuzzymembership.get(s2);
            ss = u / sumu;
        }
        return ss;
    }

    //方便别的类使用模糊邻近度矩阵
    public double[][] fuzzy_proximity_matrix(){
        double[][] fuzzyProximityMatrix = gen_fuzzy_proximity_matrix();
        return fuzzyProximityMatrix;
    }

    //将矩阵打印到屏幕
    public void printMatrix(double[][] matrix){
        for (int i = 0; i < matrix.length; i++){
            for (int j = 0; j < matrix.length; j++){
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println();
        }
    }

    //初始化
    public void Initialization(String path1, String path2){
        readInstance(path1);
        try {
            gen_fuzzy_starNeighbor(path2);
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        FeatureProximity featureproximity = new FeatureProximity();
//        featureproximity.Initialization("poi_NYC.txt", "fuzzyneighbors_NYC_(10-15).txt");
        featureproximity.Initialization("data-45.txt", "fuzzyneighbors1(300-350).txt");
        featureproximity.fuzzy_proximity_matrix();
    }
}