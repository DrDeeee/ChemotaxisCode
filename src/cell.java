/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:46
 * To change this template use File | Settings | File Templates.
 */
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Cell {


    public boolean degradeAttractant=false;
    public boolean degradedAttractantBecomesInhibitor=false;
    public boolean degradeInhibitor=false;
    public boolean degradeChemical=false;           //this affects the mechanics of the third chemical condition
    public boolean recycleReceptors=false;
    public static double recycleTime=0.3;           //this is the time it takes for bouond receptors to be internalised
    public static double noReceptors=100000;       //amount of receptors per cell
    public static double smallCellRadius=5;
    public static double clusterRadius=50;
    public static double fractionCluster=0;
    public static boolean clusterDegradesFromSurface=false;
    public static double HC = 1;
    public static double k1 = 1;
    public static double k2 = 1;
    public static double k3 =1;
    public static double speed = 10;
    public static double CIbMax =300;       //this affects cell sensitivity
    public static double a=1;               //this is fraction of attractant that conveys signal
    public static double b=0.2;               //this is fraction of inhibitor that conveys signal
    public static double c=1;               //this is the fractionl efficiency of an additional chemical that conveys signal
    public static double minOccupDiffForChem=0.001;



    public double minrad;
    public double firstPosX;
    public double[] position;
    public double[] oldposition;
    public double[] posLastPrint;
    public double[] force;
    public double[] oldforce;
    public double[] oldoldforce;
    public MigrationSimulation sim;
    public static double rspeed = Math.sqrt(speed);  // um/min
    public static double kBIRTH = 0.00031;
    public static double kDEATH = 0.00002;
    long iTicks;
    public static int nReceptors = 10000000;
    public static final double AVAGADRO =  6.022140857E23;
    public static Random RG = new Random();
    public boolean original;
    public DoubleRectangle box;
    public double neuralrad = 16.0;  // Above which pulling in
    public double maxrad = 16.0;  // Maximum, beyond which no forces
    public static boolean inducible = false;
    public boolean weakDegrader = false;
    public static double degR = 0.2; // Basal degradation
    public double deg = 0.0;
    public static double kRepel = 0.2;
    public static double kAttract = 0.1;
    public double ld  = 10.0;
    public double CIb = CIbMax;
    public double oF = 0;
    public double oB = 0;
    public double oU = 0;
    public double oL = 0;
    public boolean active = true;
    public boolean sticky = true;
    public boolean dgr = true;
    public ArrayList<EnvironmentPoint> points = new ArrayList<EnvironmentPoint>();




    //public double[] velocity;

    public double x(){ return position[0];}
    public double y(){ return position[1];}
    public double fx(){ return force[0];}
    public double fy(){ return force[1];}
    public double f2x(){ return oldforce[0];}
    public double f2y(){ return oldforce[1];}
    public double f3x(){ return oldoldforce[0];}
    public double f3y(){ return oldoldforce[1];}



    public Cell(double[] position, MigrationSimulation sim){

        this.minrad=Math.random()>fractionCluster? smallCellRadius : clusterRadius;
        this.position = position;
        this.firstPosX=position[0];
        this.oldposition=new double[]{0.0,0.0};
        this.force = new double[]{0.0,0.0};
        this.oldforce = new double[]{0.0,0.0};
        this.oldoldforce = new double[]{0.0,0.0};
        this.box = new DoubleRectangle(position[0]-maxrad,
                                       position[1]-maxrad,
                2.0*maxrad,
                2.0*maxrad   //width and height the same.
                                );

        this.posLastPrint=new double []{0.0,0.0};
        this.sim = sim;
        this.original = true;
    }

    public Cell(double[] position, MigrationSimulation sim, boolean definitelyCluster){

        if (definitelyCluster) this.minrad=clusterRadius;
        if (!definitelyCluster) this.minrad=smallCellRadius;
        this.position = position;
        this.oldposition=new double[]{0.0,0.0};
        this.force = new double[]{0.0,0.0};
        this.oldforce = new double[]{0.0,0.0};
        this.oldoldforce = new double[]{0.0,0.0};
        this.box = new DoubleRectangle(position[0]-maxrad,
                position[1]-maxrad,
                2.0*maxrad,
                2.0*maxrad   //width and height the same.
        );

        this.posLastPrint=new double []{0.0,0.0};
        this.sim = sim;
        this.original = true;
    }




    public void clear(){
        iTicks++;

        if(iTicks%200==0) {

            oldoldforce[0] = oldforce[0];
            oldoldforce[1] = oldforce[1];

            oldforce[0] = force[0];
            oldforce[1] = force[1];

        }
        force[0] = 0.0;
        force[1] = 0.0;

    }

    public void updatePosition(){

        if(!active) return;

        double npx = this.position[0]+force[0]* AgentBasedSimulation.dt;
        double npy = this.position[1]+force[1]* AgentBasedSimulation.dt;

        if(!sim.environment.GetIsOpen(npx, npy)){

            force[0] *= -0.3;
            force[1] *= -0.3;

        }

        this.oldposition[0]=this.position[0];
        this.oldposition[1]=this.position[1];

        this.position[0]+=force[0]* AgentBasedSimulation.dt;
        this.position[1]+=force[1]* AgentBasedSimulation.dt;

        this.box.x = (position[0]-minrad);
        this.box.y = (position[1]-minrad);

        //if(this.position[1]<MigrationSimulation.padding||this.position[1]> sim.yMax-MigrationSimulation.padding) {
        //    this.position[1] = MyMaths.bounded(MigrationSimulation.padding, sim.yMax - MigrationSimulation.padding, this.position[1]);
        //}
        //if(this.position[1]<MigrationSimulation.padding) this.position[1] = 2.0*MigrationSimulation.padding-this.position[1];
        //else if(this.position[1]> sim.yMax-MigrationSimulation.padding) this.position[1] = 2.0* (sim.yMax-MigrationSimulation.padding) - this.position[1];


        //if(this.position[0]<MigrationSimulation.padding) this.position[0] = 2.0*MigrationSimulation.padding-this.position[0];
        //else if(this.position[0]> sim.xMax-MigrationSimulation.padding){
        //    this.position[0] = 2.0* (sim.xMax-MigrationSimulation.padding) - this.position[0];
            //System.out.println((sim.Ttotal/60));
        //}

    }

    public void addForce(double dx, double dy){
        force[0]+=dx;
        force[1]+=dy;
    }

    public void GetEnvironmentPointsInCell(ChemicalEnvironment env){
        points.clear();

        double dx = ChemicalEnvironment.grain;
        int iX = Math.max(0,(int) Math.floor((position[0]-minrad)/dx));
        int iY = Math.max(0,(int) Math.floor((position[1]-minrad)/dx));

        int iW = (int) Math.ceil(2.0*minrad/dx);

        for(int i = iX; i<=iX+iW; i++){
            for(int j = iY; j<=iY+iW; j++){
                if((((i*dx-position[0])*(i*dx-position[0])) + ((j*dx-position[1])*(j*dx-position[1])))<=minrad*minrad){
                    int i2 = Math.max(0,Math.min(env.profile.size()-1, i));
                    int j2 = Math.max(0,Math.min(env.profile.get(i2).size()-1, j));
                    EnvironmentPoint ep = env.profile.get(i2).get(j2);
                    points.add(ep);
                }
            }
        }

    }


    public void DegradeFromEnvironment(ChemicalEnvironment env){
        points.clear();
        GetEnvironmentPointsInCell(env);


        double km = ChemicalEnvironment.kM;
        double dt = AgentBasedSimulation.dt;
        double k1, k2, k3, k4;
        double fr;

        if (minrad==smallCellRadius){
            fr=(ChemicalEnvironment.sMax)/points.size();                // A small cells degrading power split among points in cells cross sectional area
        }else {
            if (clusterDegradesFromSurface){
                fr=((2*Math.PI*Math.pow(clusterRadius,2))/(Math.PI*Math.pow(smallCellRadius,2)))*(ChemicalEnvironment.sMax/points.size());   //A clusters degrading power (based upon a hemisphere with degrading small cells on its surface) split among points
            }else{
                fr=(((2.0/3)*Math.PI*Math.pow(clusterRadius,3))/((2.0/3)*Math.PI*Math.pow(smallCellRadius,3)))*(ChemicalEnvironment.sMax/points.size());  //A clusters degrading power (based upon a hemisphere with degrading small cells in its volume) split among points
            }
        }


        for(EnvironmentPoint e : points) {

            if (degradeInhibitor) {

                k1 = fr * Math.pow(e.c2, HC) / (Math.pow(e.c2, HC) + km);

                k2 = (e.c2 + 0.5 * dt * k1);
                k2 = fr * Math.pow(k2, HC) / (Math.pow(k2, HC) + km);

                k3 = (e.c2 + 0.5 * dt * k2);
                k3 = fr * Math.pow(k3, HC) / (Math.pow(k3, HC) + km);

                k4 = e.c2 + dt * k3;
                k4 = fr * Math.pow(k4, HC) / (Math.pow(k4, HC) + km);

                double inhibDegraded = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1 + 2 * k2 + 2 * k3 + k4);

                if (!e.fixed && e.open)
                    e.c2 -= inhibDegraded;
            }

            //e.c+=ChemicalEnvironment.baseConcentration*AgentBasedSimulation.dt/points.size();
            // If the point allows free, unrestricted diffusion, degrade.

            //RK4 integration.

            if (degradeAttractant) {
                k1 = fr * Math.pow(e.c, HC) / (Math.pow(e.c, HC) + km);

                k2 = (e.c + 0.5 * dt * k1);
                k2 = fr * Math.pow(k2, HC) / (Math.pow(k2, HC) + km);

                k3 = (e.c + 0.5 * dt * k2);
                k3 = fr * Math.pow(k3, HC) / (Math.pow(k3, HC) + km);

                k4 = e.c + dt * k3;
                k4 = fr * Math.pow(k4, HC) / (Math.pow(k4, HC) + km);

                double attrDegraded = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1 + 2 * k2 + 2 * k3 + k4);


                if (!e.fixed && e.open) {
                    e.c -= attrDegraded;


                    if (degradedAttractantBecomesInhibitor) {
                        e.c2 += attrDegraded;
                        //e.c2_m1 += degraded;
                    }
                }
            }



            if (degradeChemical) {

                k1 = fr * Math.pow(e.c3, HC) / (Math.pow(e.c3, HC) + km);

                k2 = (e.c3 + 0.5 * dt * k1);
                k2 = fr * Math.pow(k2, HC) / (Math.pow(k2, HC) + km);

                k3 = (e.c3 + 0.5 * dt * k2);
                k3 = fr * Math.pow(k3, HC) / (Math.pow(k3, HC) + km);

                k4 = e.c3 + dt * k3;
                k4 = fr * Math.pow(k4, HC) / (Math.pow(k4, HC) + km);

                double inhibDegraded = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1 + 2 * k2 + 2 * k3 + k4);

                if (!e.fixed && e.open)
                    e.c3 -= inhibDegraded;
            }



        }
    }

    public void RecycleReceptors (ChemicalEnvironment env){

        double dt = AgentBasedSimulation.dt;
        ArrayList<Double> averageChemBinding = this.averageChemBinding(env);

        double noMoleculesChem1= averageChemBinding.get(0)*noReceptors;
        double noMoleculesChem2= averageChemBinding.get(1)*noReceptors;
        double noMoleculesChem3= averageChemBinding.get(2)*noReceptors;

        double microMolesChem1Bound=(noMoleculesChem1/AVAGADRO)*1000000;
        double microMolesChem2Bound=(noMoleculesChem2/AVAGADRO)*1000000;
        double microMolesChem3Bound=(noMoleculesChem3/AVAGADRO)*1000000;

        double volumeCell = (4/3)*Math.PI*Math.pow(smallCellRadius,3);
        double volumeCellInLitres = volumeCell*1E-15;

        double concChem1Bound= microMolesChem1Bound/volumeCellInLitres;
        double concChem2Bound= microMolesChem2Bound/volumeCellInLitres;
        double concChem3Bound= microMolesChem3Bound/volumeCellInLitres;

        double amountDegradedChem1=concChem1Bound*(dt/recycleTime);
        double amountDegradedChem2=concChem2Bound*(dt/recycleTime);
        double amountDegradedChem3=concChem3Bound*(dt/recycleTime);

        points.clear();
        GetEnvironmentPointsInCell(env);

        if (recycleReceptors) {
            for(int i=0; i<points.size(); i++) {
                if (!points.get(i).fixed && points.get(i).open)
                    points.get(i).c -= amountDegradedChem1/points.size();
                    points.get(i).c2 -= amountDegradedChem2/points.size();
                    points.get(i).c3 -= amountDegradedChem3/points.size();
            }
        }
    }











    public double[] OccupancyDiffAcrossXY (ChemicalEnvironment env) {

        GetEnvironmentPointsInCell(env);

        ArrayList<Double> XinCell = new ArrayList<Double>();
        ArrayList<Double> YinCell = new ArrayList<Double>();

        double xMax;
        double xMin;
        double yMax;
        double yMin;
        double occupMaxX;
        double occupMinX;
        double occupMaxY;
        double occupMinY;




            for (int i = 0; i < points.size(); i++) {
                XinCell.add(points.get(i).x);
                YinCell.add(points.get(i).y);
            }

            Collections.sort(XinCell);
            Collections.sort(YinCell);

        if (XinCell.size()!=0 && YinCell.size()!=0) {
            xMax = Collections.max(XinCell);
            xMin = Collections.min(XinCell);
            yMax = Collections.max(YinCell);
            yMin = Collections.min(YinCell);
        }else{
            return new double [] {0.0,0.0};
        }


            double xMid=XinCell.get(XinCell.size()/2);
            double yMid=YinCell.get(YinCell.size()/2);



            EnvironmentPoint maxX= env.profile.get((int)(xMax/ChemicalEnvironment.grain)).get((int)(yMid/ChemicalEnvironment.grain));
            EnvironmentPoint minX= env.profile.get((int)(xMin/ChemicalEnvironment.grain)).get((int)(yMid/ChemicalEnvironment.grain));
            EnvironmentPoint maxY= env.profile.get((int)(xMid/ChemicalEnvironment.grain)).get((int)(yMax/ChemicalEnvironment.grain));
            EnvironmentPoint minY= env.profile.get((int)(xMid/ChemicalEnvironment.grain)).get((int)(yMin/ChemicalEnvironment.grain));




        occupMaxX = (a*maxX.c+(b*maxX.c2*(k1 / k2))) / ((maxX.c + k1 * (1 + maxX.c2 / k2)));
        occupMinX = (a*minX.c+(b*minX.c2*(k1 / k2))) / ((minX.c + k1 * (1 + minX.c2 / k2)));
        occupMaxY = (a*maxY.c+(b*maxY.c2*(k1 / k2))) / ((maxY.c + k1 * (1 + maxY.c2 / k2)));
        occupMinY = (a*minY.c+(b*minY.c2*(k1 / k2))) / ((minY.c + k1 * (1 + minY.c2 / k2)));


            return new double[]{Math.abs(occupMaxX-occupMinX),Math.abs(occupMaxY-occupMinY)};


    }








    public double cosTheta (){

        if (AgentBasedSimulation.introduceCellsNow = true){
            return 0;
        }

        if (posLastPrint[0]==0.0){
        posLastPrint[0]=this.position[0];
        posLastPrint[1]=this.position[1];
        return 0;}

        double changeInX=position[0]-posLastPrint[0];
        double changeinY=position[1]-posLastPrint[1];

        double theta=Math.atan2(changeinY,changeInX);
        double cosTh=Math.cos(theta);

        posLastPrint[0]=position[0];
        posLastPrint[1]=position[1];

        return cosTh;
    }

    public double xDirectionVelocity () {

        if (posLastPrint[0]==0){
            posLastPrint[0]=this.position[0];
            posLastPrint[1]=this.position[1];
            return 0;}

        double changeInX=this.position[0]-posLastPrint[0];
        double xVelocity=(changeInX)/AgentBasedSimulation.outInt;

        posLastPrint[0]=position[0];
        posLastPrint[1]=position[1];

        return xVelocity;
    }


    public double averageOccupancy (ChemicalEnvironment env){

        GetEnvironmentPointsInCell(env);
        double sum=0;
        double occupancy;

        if (points.size()!=0){
            for (int i=0; i<points.size(); i++){
              occupancy=(a*points.get(i).c +(b*points.get(i).c2*(k1 / k2)))/ ((points.get(i).c + k1 * (1 + points.get(i).c2 / k2)));
              sum += occupancy;
            }
        }else {
            return 0.0;
        }

        return sum/points.size();
    }


    public ArrayList<Double> averageChemBinding (ChemicalEnvironment env){

        ArrayList<Double> averageBindingEachChem = new ArrayList<Double>();
        points.clear();
        GetEnvironmentPointsInCell(env);
        double sum1=0;
        double sum2=0;
        double sum3=0;
        double bound1;
        double bound2;
        double bound3;

        if (points.size() !=0) {
            for (int i = 0; i < points.size(); i++) {
                bound1 = (points.get(i).c / k1) / (1 + (points.get(i).c / k1) + (points.get(i).c2 / k2) + (points.get(i).c3 / k3));
                bound2 = (points.get(i).c2 / k2) / (1 + (points.get(i).c / k1) + (points.get(i).c2 / k2) + (points.get(i).c3 / k3));
                bound3 = (points.get(i).c3 / k3) / (1 + (points.get(i).c / k1) + (points.get(i).c2 / k2) + (points.get(i).c3 / k3));
                sum1 += bound1;
                sum2 += bound2;
                sum3 += bound3;
            }
        }else {
            ArrayList<Double> fail = new ArrayList<Double>();
            fail.add(0.0);
            fail.add(0.0);
            fail.add(0.0);
                return fail;
            }

            averageBindingEachChem.add(sum1/points.size());
            averageBindingEachChem.add(sum2/points.size());
            averageBindingEachChem.add(sum3/points.size());


            return averageBindingEachChem;

    }






    public double[] EstimateGradientDirection(ChemicalEnvironment env){

        GetEnvironmentPointsInCell(env);




        double estX = 0;
        double estY = 0;
        double cX = 0;
        double cY = 0;

        double xMin = 0;
        double xMax = 0;
        double yMin = 0;
        double yMax = 0;


        for(EnvironmentPoint ep : points){
            cX += ep.x;
            cY += ep.y;
        }

        cX/=points.size();
        cY/=points.size();

        double occupancy;
        double meanoccupancy = 0;
        double nT;
        double nR;
        double dC;

        for(EnvironmentPoint ep : points) {
            // Calculate directional projection of occupancy across the cell
            occupancy           = ((a*ep.c/k1)+(b*ep.c2/k2)+(c*ep.c3/k3))/(1+((ep.c/k1)+(ep.c2/k2)+(ep.c3/k3)));               //this calculates agonistsic occupancy for a point
            meanoccupancy      += occupancy/points.size();
            estX               += occupancy * (ep.x - cX) / points.size();
            estY               += occupancy * (ep.y - cY) / points.size();
        }


        if(inducible) {

            double x = meanoccupancy;
            //x = k1 * x / (1-x);                                 // These run if
            //x = (x*x)/(k1+(x*x));                               // ultrasensitive.

            deg = x;                                            // If no delay
            //deg += degR * (x - deg) * MelaMigration.dt;       // If delay
        }
        else deg = 1;

        // Remove environmental attractant based on occupancy, area, depth, dt.

        //number of molecules in area
        /*    if(ep.open&&!ep.fixed) {

                nT = ChemicalEnvironment.environment_depth * ChemicalEnvironment.grain
                        * ChemicalEnvironment.grain * 1E-15 * ChemicalEnvironment.AVAGADRO * ep.c * 1E-6;
                // Number of occupied receptors in THIS part of the cell.
                nR = (nReceptors / points.size()) * occupancy;

                dC = Math.min(1, nR / nT) * ep.c * MelaMigration.dt;

                //System.out.println("nT -> "+nT+",    nR-> "+nR+",   dC -> "+dC);
                // Removed fraction of total concentration involved in gathering this information!
                ep.c -= dC;
            } */

        return new double[] {estX,estY};


        /*double cMax = 0;
        EnvironmentPoint epo = points.get(0);

        for(EnvironmentPoint ep : points){
            double cr = ep.c+Math.sqrt(ep.c)*RG.nextGaussian();

            if(cr>cMax){
                cMax = cr;
                epo = ep;
            }
        }

        if(cMax<0.01) return new double[]{0.,0.};
        else return new double[]{epo.x-this.x(),epo.y-this.y()};   */
    }

    public boolean checkExit(){
        for(EnvironmentPoint ep : points){
            if(ep.fixed) return true;
        }
        return false;
    }

    public void updateGrowth(){

        if(sim.dieoff)      this.cull(0.0,0.0,4.0);
        if(sim.proliferate) this.split(0.0,0.0, AgentBasedSimulation.dimensions[0]);

    }


    public synchronized void cull(double a, double b, double c){

        double t1 = a*x()*x();
        double t2 = b*x();
        double t3 = c;

        if((kDEATH* AgentBasedSimulation.dt>Math.random())&&((t1+t2+t3)<Math.random()* AgentBasedSimulation.dimensions[0])){
            sim.cells.remove(this);
        }
    }

    public synchronized void split(double a, double b, double c){

        if((kBIRTH* AgentBasedSimulation.dt>Math.random())/*&&((t1+t2+t3)>Math.random()*MelaMigration.dimensions[0])*/){
            Cell clone = new Cell(new double[] {position[0] + 2.0*minrad*RG.nextGaussian(), position[1] + 2.0*minrad*RG.nextGaussian()},this.sim);
            clone.position[1] = MyMaths.bounded(0.0, AgentBasedSimulation.dimensions[1]-0.0, clone.position[1]);
            if(clone.position[0]<0.1) clone.position[0] = 0.2-clone.position[0];
            else if(clone.position[0]> AgentBasedSimulation.dimensions[0]-0.1) clone.position[0] = 2.0* (AgentBasedSimulation.dimensions[0]-0.1) - clone.position[0];


            sim.newCells.add(clone);
            clone.original = false;
        }
    }
}
