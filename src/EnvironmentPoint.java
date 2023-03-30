/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 04/12/15
 * Time: 12:32
 * To change this template use File | Settings | File Templates.
 */
public class EnvironmentPoint {

    public double c = 0.0;
    public double c_m1 = 0.0;
    public double c_m2 = 0.0;

    public double c2 = 0.0;
    public double c2_m1 = 0.0;
    public double c2_m2 = 0.0;

    public double c3 = 0.0;
    public double c3_m1 = 0.0;
    public double c3_m2 = 0.0;

    public boolean open = true;
    public boolean fixed = false;
    public boolean start = false;
    public boolean wrong = false;
    public boolean right = false;

    public double x;
    public double y;

    public EnvironmentPoint xp1;
    public EnvironmentPoint xm1;
    public EnvironmentPoint yp1;
    public EnvironmentPoint ym1;

    public EnvironmentPoint(int x, int y) {

        this.x = x * ChemicalEnvironment.grain;
        this.y = y * ChemicalEnvironment.grain;

    }

}




