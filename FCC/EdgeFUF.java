package FC;

public class EdgeFUF {
    private String edge;
    private double fuf;

    public EdgeFUF(String edge, double fuf){
        this.edge = edge;
        this.fuf = fuf;
    }

    public String getEdge() {
        return edge;
    }

    public double getFuf() {
        return fuf;
    }

    public void setEdge(String edge) {
        this.edge = edge;
    }

    public void setFuf(double fuf) {
        this.fuf = fuf;
    }
}
