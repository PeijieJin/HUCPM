package FDP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FuzzyUtilityRatio {
    private HashMap<String, Map<String, Double>> featureInstance; //特征及其实例与效用集
    private HashMap<String, Integer> featureCount; //特征的实例数量
    private HashMap<String, Map<String, Map<String, Double>>> fuzzyStarNeighbor; //模糊星型实例集
    private HashMap<String, Double> fuzzyMembership; //模糊成员集
    private HashMap<String, List<String>> fuzzyNeighbor; //模糊邻近关系集
    public FuzzyUtilityRatio(){
        featureInstance = new HashMap<>();
        featureCount = new HashMap<>();
        fuzzyStarNeighbor = new HashMap<>();
        fuzzyMembership = new HashMap<>();
        fuzzyNeighbor = new HashMap<>();
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

    //生成模糊星型实例集
    //path路径中存储的文件形式为“f1.id1,utility(f1.id1),f2.id2,utility(f2.id2)，u”
    private void gen_fuzzy_starNeighbor(String path){
        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(path)));
            String line = null;

            while ((line = reader.readLine()) != null){
                if (line.equals("")){
                    continue;
                }

                String[] tmp = line.split(",");
                String ins1 = tmp[0]; //取出第一个实例
                int index1 = ins1.indexOf(".");
                String f1 = ins1.substring(0, index1); //第一个实例的特征
                String ins2 = tmp[2]; //取出第二个实例
                int index2 = ins2.indexOf(".");
                String f2 = ins2.substring(0, index2); //第二个实例的特征
                double u = Double.parseDouble(tmp[4]); //模糊邻近隶属度

                //生成模糊星型实例集和模糊隶属度和集
                //存储<ins1，,<f2,<ins2,,u>>>
                //存储<f2,<ins2,,u>>
                HashMap<String, Double> ins2u = new HashMap<>();
                ins2u.put(ins2, u);
                HashMap<String, Map<String, Double>> f2ins2u = new HashMap<>();
                f2ins2u.put(f2, ins2u);
                if (fuzzyStarNeighbor.containsKey(ins1)){
                    if (fuzzyStarNeighbor.get(ins1).containsKey(f2)){
                        if (!fuzzyStarNeighbor.get(ins1).get(f2).containsKey(ins2)){
                            fuzzyStarNeighbor.get(ins1).get(f2).put(ins2, u);
                        }
                    } else {
                        fuzzyStarNeighbor.get(ins1).put(f2, ins2u);
                    }
                } else {
                    fuzzyStarNeighbor.put(ins1, f2ins2u);
                }

                //存储<ins2,<f1,<ins1,u>>>
                //存储<f1,<ins1,u>>
                HashMap<String, Double> ins1u = new HashMap<>();
                ins1u.put(ins1, u);
                HashMap<String, Map<String, Double>> f1ins1u = new HashMap<>();
                f1ins1u.put(f1, ins1u);
                if (fuzzyStarNeighbor.containsKey(ins2)){
                    if (fuzzyStarNeighbor.get(ins2).containsKey(f1)){
                        if (!fuzzyStarNeighbor.get(ins2).get(f1).containsKey(ins1)){
                            fuzzyStarNeighbor.get(ins2).get(f1).put(ins1, u);
                        }
                    } else {
                        fuzzyStarNeighbor.get(ins2).put(f1, ins1u);
                    }
                } else {
                    fuzzyStarNeighbor.put(ins2, f1ins1u);
                }

                //模糊成员集，将邻近关系及其隶属度存储起来
                String member = ins1 + "," + ins2;
                if (!fuzzyMembership.containsKey(member)){
                    fuzzyMembership.put(member, u);
                }

                //模糊邻近关系，将不同特征之间的邻近关系存储起来
                String fmember = f1 + "," + f2;
                if (!fuzzyNeighbor.containsKey(fmember)){
                    fuzzyNeighbor.put(fmember, new ArrayList<>());
                }
                fuzzyNeighbor.get(fmember).add(member);
            }
            reader.close();
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    //生成模糊效用率矩阵fuzzyUtilityRatioMatrix
    private double[][] gen_fuzzy_utility_ratio_matrix(){
        int size = featureCount.size();
        double[][] fuzzyUtilityRatioMatrix = new double[size][size];
        List<String> fList = new ArrayList<>();
        fList.addAll(featureCount.keySet()); //特征集
        for (int i = 0; i < size; i++){
            String f1 = fList.get(i);
            for (int j = i ; j < size; j++){
                String f2 = fList.get(j);
                fuzzyUtilityRatioMatrix[i][j] = computefuzzyUtilityRatio(f1, f2);
            }
        }

        for (int i = 1; i < size; i++){
            for (int j = 0; j < i; j++){
                fuzzyUtilityRatioMatrix[i][j] = fuzzyUtilityRatioMatrix[j][i];
            }
        }
//        System.out.println("模糊效用率矩阵：");
//        printMatrix(fuzzyUtilityRatioMatrix);
        return fuzzyUtilityRatioMatrix;
    }

    private double computefuzzyUtilityRatio(String f1, String f2) {
        double fuf = 0.0;
        double sumuti = 0.0, sumfu = 0.0;
        List<String> fnList = new ArrayList<>();
        String fmember = f1 + "," + f2;
        //如果f1和f2的实例之间没有模糊邻近关系fuf赋值为0
        if (!fuzzyNeighbor.containsKey(fmember)){
            return fuf;
        }
        fnList.addAll(fuzzyNeighbor.get(fmember)); //取出fi和fj所有的模糊邻近关系
        //计算模糊效用率
        for (int i = 0; i < fnList.size(); i++){
            String member = fnList.get(i);
            String[] tmp = member.split(",");
            String ins1 = tmp[0];
            String ins2 = tmp[1];

            //fi和fj所有模糊邻近关系实际效用总和
            double uti1 = featureInstance.get(f1).get(ins1);
            double uti2 = featureInstance.get(f2).get(ins2);
            sumuti += uti1 + uti2;

            //fi和fj所有模糊邻近关系的实例间的模糊效用总和
            sumfu += computeFuzzyUtility(ins1, ins2);
        }
        fuf = sumfu / sumuti;
        if (fuf > 1){
            fuf = 1;
        }

        return fuf;
    }

    private double computeFuzzyUtility(String ins1, String ins2) {
        double fu = 0.0;
        int index1 = ins1.indexOf(".");
        String f1 = ins1.substring(0, index1); //ins1的特征
        int index2 = ins2.indexOf(".");
        String f2 = ins2.substring(0, index2); //ins2的特征
        String member = ins1 + "," + ins2;

        double uti1 = featureInstance.get(f1).get(ins1); //取出ins1的效用
        double uti2 = featureInstance.get(f2).get(ins2); //取出ins2的效用

        fu = computeContriScore(ins1, ins2) * uti1 + computeContriScore(ins2, ins1) * uti2;

        return fu;
    }

    private double computeContriScore(String ins1, String ins2) {
        double cu = 0.0;
        double sumcu = 0.0;
        int index = ins2.indexOf(".");
        String f = ins2.substring(0, index); //取出ins2的特征
        HashMap<String, Double> mapList = new HashMap<>();
        mapList.putAll(fuzzyStarNeighbor.get(ins1).get(f)); //取出和ins2一同共享ins1的所有实例及隶属度
        if (mapList.size() == 1){
            cu = mapList.get(ins2);
        } else {
            for (Map.Entry<String, Double> entry : mapList.entrySet()){
                sumcu += entry.getValue();
            }
            cu = mapList.get(ins2) / sumcu;
        }

        return cu;
    }

    //方便别的类使用模糊效用率矩阵
    public double[][] fuzzyUtilityRatioMatrix(){
        double[][] fuzzyUtilityRatioMatrix = gen_fuzzy_utility_ratio_matrix();
        return fuzzyUtilityRatioMatrix;
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
        FuzzyUtilityRatio fur = new FuzzyUtilityRatio();
        fur.Initialization("gen/jInitialData/data-20-8000(utility).txt", "gen/jInitialData/fuzzyneighbor(50)_data-20-8000(utility).txt");
        fur.fuzzyUtilityRatioMatrix();
    }
}
