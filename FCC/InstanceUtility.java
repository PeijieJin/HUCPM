package FC;

public class InstanceUtility {
    private String instance;
    private double utility;
    public InstanceUtility(String instance, double utility){
        this.instance = instance;
        this.utility = utility;
    }

    public String getInstance() {
        return instance;
    }

    public double getUtility() {
        return utility;
    }

    @Override
    public int hashCode() {
        return instance!=null&&utility!=0?instance.hashCode()+((Double)utility).hashCode():0;
    }


    @Override
    public boolean equals(Object object){
        if(object == this){
            return true;
        }
        if(object == null || getClass()!= object.getClass()){
            return false;
        }
        InstanceUtility ins  = (InstanceUtility) object;
        if(instance != null? !instance.equals(ins.instance):ins.instance!=null){
            return false;
        }
        if(utility!=0?utility!=(ins.utility):ins.utility!=0){
            return  false;
        }
        return true;
    }

}
