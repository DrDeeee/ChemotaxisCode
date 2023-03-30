//import com.amd.aparapi.Kernel;
//import com.amd.aparapi.Range;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 19/02/2014
 * Time: 08:23
 * To change this template use File | Settings | File Templates.
 */
public class ChemicalEnvironment {

    public static boolean diffuse =true;
    public static boolean diffuseAndBackground=true;                //gives a background for each chemical defined by the lowest value out of the source and sink (otherwise 0)
    public static boolean noSinkOrSource=false;                      //eliminated the fixed concentration at either end of the bridge and gives the area limited resources


    public boolean sourceSinkDegradation=false;                     //only degrades conc1 source and sink
    public boolean degradationInSourceGivesInhibitor=false;         //adds to conc2 source and sink
    public boolean fluxIntoOutOfSourceAndSink=false;                //adds/minus from conc1/conc2 source/sink based on diffusive flux

    public static double concSink1 =1;                              //sink means LHS, source means RHS
    public static double concSource1 =5;
    public static double concSink2 =0;
    public static double concSource2 =0;
    public static double concSink3 =0;
    public static double concSource3 =0;

    public static double gam1=1;            //if chemical 1 has a shape that can be approximated by a power then gam1 defines its steepness (2=quadratic, 3=cubic etc...)
    public static double gam2=1;
    public static double gam3=1;

    public static double r1=0.01;              //if chemical 1 is subject to degradation and C<<KM then its shape can be approximated to an exponential where r1 contains information about km, smax and D
    public static double r2=0.000193;


    public static double q1=0.000005;                //if chem1 subject to degradation and C>>kM then its shape can  be approximated by a quadratic modified by q (which contains information on smax and D)




    public static boolean sourceAtTopAndBottom=false;            //sink means top, source means bottom
    public static double concSink12 =0;
    public static double concSource12 =0.03;
    public static double concSink22 =0;
    public static double concSource22 =0;
    public static double concSink32 =0;
    public static double concSource32 =0;


    public static double bridgeDepth=20;      //need this to put degradation in uMoles per min (not uMolar!)
    public static double resevoirDepth=2000;
    public static int noCellsInSource=10;
    public static double sMax   = 1000;           //remember difference between continuum and discrete vmax
    public static double kM     = 2;
    public static double  DiffC  =20000;
    public static double  grain = 3.5;





    public static final double AVAGADRO =  6.022140857E23;
    public static double  L = 1600.;      //
    public static double  W = 15.0;       //
    public static int     k = 1;       // Num branch levels
    public static double  shrinkage = 1;
    public static double  skew  = -2.0;       // (70.7 <-> 2262)
    public static double  trunk = 250;
    public static double  asym = 1;
    public static double  padding = 100;



    public double[][] chemvironment = new double[(int)(AgentBasedSimulation.dimensions[2]/grain)][(int) (AgentBasedSimulation.dimensions[3]/grain)];
    public ArrayList<ArrayList<EnvironmentPoint>> profile;
    public ArrayList<EnvironmentPoint> freepoints;


    public static double production  = 0.0;
    public static double degradation = 0.0; // Equilibrium values of 20.

    private ChemicalEnvironment(){

        profile = new ArrayList<ArrayList<EnvironmentPoint>>();

    }


    public static ChemicalEnvironment SetUpEnvironment(){

        ChemicalEnvironment ce = new ChemicalEnvironment();

        int iMax = ((int) (AgentBasedSimulation.dimensions[2]/grain));
        int jMax = ((int) (AgentBasedSimulation.dimensions[3]/grain));

        for(int i = 0; i<iMax; i++){

            ArrayList<EnvironmentPoint>  alENV = new ArrayList<EnvironmentPoint>();
            for(int j = 0; j<jMax; j++){
                EnvironmentPoint ep = new EnvironmentPoint(i,j);
                //if(i*grain<100) ep.start = true;




            if (diffuse && !diffuseAndBackground) {

                if (!sourceAtTopAndBottom) {
                    if (i == iMax - 1) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSource1;
                        ep.c2 = ep.c2_m1 = concSource2;
                        ep.c3 = ep.c3_m1 = concSource3;
                    } else if (i == 0) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSink1;
                        ep.c2 = ep.c2_m1 = concSink2;
                        ep.c3 = ep.c3_m1 = concSink3;
                    } else {
                        ep.open = true;
                        ep.c = ep.c_m1 = 0;
                        ep.c2 = ep.c2_m1 = 0;
                        ep.c3 = ep.c3_m1 = 0;
                    }
                }


                if (sourceAtTopAndBottom) {
                    if (j == jMax - 1) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSource12;
                        ep.c2 = ep.c2_m1 = concSource22;
                        ep.c3 = ep.c3_m1 = concSource32;
                    } else if (j == 0) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSink12;
                        ep.c2 = ep.c2_m1 = concSink22;
                        ep.c3 = ep.c3_m1 = concSink32;
                    } else if (i == iMax - 1) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSource1;
                        ep.c2 = ep.c2_m1 = concSource2;
                        ep.c3 = ep.c3_m1 = concSource3;
                    } else if (i == 0) {
                        ep.fixed = true;
                        ep.c = ep.c_m1 = concSink1;
                        ep.c2 = ep.c2_m1 = concSink2;
                        ep.c3 = ep.c3_m1 = concSink3;
                    } else {
                        ep.open = true;
                        ep.c = ep.c_m1 = 0;
                        ep.c2 = ep.c2_m1 = 0;
                        ep.c3 = ep.c3_m1 = 0;
                    }
                }



            }

            if (diffuse && noSinkOrSource){
                ep.open = true;
                ep.c = ep.c_m1 = concSink1 +(Math.pow((1.0*i/iMax),gam1)*(concSource1 - concSink1));
                ep.c2 =ep.c2_m1 = concSink2 +(Math.pow((1.0*i/iMax),gam2)*(concSource2 - concSink2));
                ep.c3 = ep.c3_m1 =concSource3;
                }



            else if (diffuse && diffuseAndBackground) {

                if (i == iMax - 1) {
                    ep.fixed = true;
                    ep.c = ep.c_m1 = concSource1;
                    ep.c2 = ep.c2_m1 = concSource2;
                    ep.c3 =ep.c3_m1 =concSource3;
                } else if (i == 0) {
                    ep.fixed = true;
                    ep.c = ep.c_m1 = concSink1;
                    ep.c2 = ep.c2_m1 = concSink2;
                    ep.c3 =ep.c3_m1 =concSink3;
                } else{
                    ep.open = true;
                    ep.c = ep.c_m1 = Math.min(concSink1, concSource1);
                    ep.c2 = ep.c2_m1 = Math.min(concSink2, concSource2);
                    ep.c3 = ep.c3_m1 = Math.min(concSink3, concSource3);
                }







            }else {

                //ep.c = ep.c_m1 =concSource1;
                //ep.c =ep.c_m1 = concSink1+((1.0*i/iMax)*(concSource1-concSink1));
                //ep.c = ep.c_m1 = concSink1 +(Math.pow((1.0*i/iMax),gam1)*(concSource1 - concSink1));               //power law approximation to a chemical gradient with steepness gam1
                //ep.c=ep.c_m1=(concSource1/(Math.exp(r1*iMax*grain)-Math.exp(-r1*iMax*grain)))*(Math.exp(r1*i*grain)-Math.exp(-r1*i*grain));

                                                                                                                    //exponential approximation to a gradient subject to degradation where C<<kM
                                                                                                                    //explicitely starting from 0
                                                                                                                    //r1 defines shape of curve according to degradation mechanics and so contains
                                                                                                                    //information about smax, km and D. specifically r1=smax/(km*D)

                ep.c=ep.c_m1=(q1*Math.pow(i*grain,2))+(((i*grain)/(iMax*grain))*(concSource1-(q1*Math.pow(iMax*grain,2))));
                                                                                                                    //quadratic approximation to gradient subject to degradation in limit c>>km
                                                                                                                    //gradient explicitely starting from 0 and going to concSource1




                //ep.c2 = ep.c2_m1 =concSource2;
                ep.c2 =ep.c2_m1 = concSink2 +(Math.pow((1.0*i/iMax),gam2)*(concSource2 - concSink2));
                //ep.c2 = ep.c2_m1 =concSink2+(Math.pow((1.0*i/iMax),2)*(concSource2-concSink2));
                //ep.c2=ep.c2_m1=(concSource2/(Math.exp(r2*iMax*grain)-Math.exp(-r2*iMax*grain)))*(Math.exp(r2*i*grain)-Math.exp(-r2*i*grain));

                //ep.c3 = ep.c3_m1 =concSource3;
                ep.c3 =ep.c3_m1 = concSink3 +(Math.pow((1.0*i/iMax),gam3)*(concSource3 - concSink3));
                //ep.c3 = ep.c3_m1 =concSink3+(Math.pow((1.0*i/iMax),2)*(concSource3-concSink3));

            }




                if(AgentBasedSimulation.pinned) {
                    //if(i==0||i==(iMax-1)) ep.fixed = true;
                }
                alENV.add(ep);
            }
            ce.profile.add(alENV);

        }
        ce.SetUpDiffusion();
        return ce;
    }

    public static ChemicalEnvironment EnvironmentFromFile(String imageFile){

        ChemicalEnvironment ce = new ChemicalEnvironment();

        File f = new File(imageFile);
        BufferedImage img;

        if(!f.exists()){
            System.out.println("No such file as "+imageFile);
            ce = SetUpEnvironment();
            return ce;
        }
        try{
            img =  ImageIO.read(f);
        }
        catch(IOException e){
            System.out.println("Could not open "+imageFile);
            ce = SetUpEnvironment();
            return ce;
        }
        for(int i = 0; i < img.getWidth(); i++){

            ArrayList<EnvironmentPoint>  alENV = new ArrayList<EnvironmentPoint>();

            for(int j = 0; j < img.getHeight(); j++){

                EnvironmentPoint ep = new EnvironmentPoint(i,j);

                int clr = img.getRGB(i,j);

                int blue = clr & 0xff;
                int green = (clr & 0xff00) >> 8;
                int red = (clr & 0xff0000) >> 16;



                if(red<20&&green<20&&blue<20)  ep.open  = false;
                if(red>200&&green<20&&blue<20) ep.fixed = true;
                if(red<20&&green<20&&blue>200) {
                    ep.start = true;
                    if(AgentBasedSimulation.pinned){
                        ep.fixed = true;
                        ep.c = 0;
                    }

                }
                if(ep.open && !(AgentBasedSimulation.pinned && ep.start)) ep.c = ep.c_m1 = ep.c_m2    = concSource1;

                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }

        ce.SetUpDiffusion();
        return ce;
    }

    public void DrawLine(double[] start, double[] end, double width, boolean entry, boolean fixed, boolean wrong, boolean right){     // Only works for u/d l/r

        int[] iStart = new int[] {(int) (start[0]/grain),(int) (start[1]/grain)};          // Convert coordinates to
        int[] iEnd   = new int[] {(int) (end[0]/grain),  (int) (end[1]/grain)};

        boolean bH = true;

        int iW = (int) Math.round(width/grain);

        if(iStart[0]==iEnd[0]) bH = false;

        if(bH){

            //int i1 = iStart[0]< iEnd[0] ? iStart[0] : iEnd[1];
            //int i2 = iStart[0]>=iEnd[0] ? iStart[0] : iEnd[1];

            for(int i=iStart[0]-iW; i<=iEnd[0]+iW; i++){
                for(int j=iStart[1]-iW; j<=iEnd[1]+iW; j++){
                    EnvironmentPoint ep = profile.get(i).get(j);
                    ep.open = true;
                    ep.start = entry;
                    ep.fixed = fixed;
                    ep.wrong = wrong;
                    ep.right = right;
                    ep.c = ep.c_m1    = concSource1;
                }
            }

        }
        else{
            for(int i=iStart[0]-iW; i<=iEnd[0]+iW; i++){

                //int j1 = iStart[0]< iEnd[0] ? iStart[0] : iEnd[0];
                //int j2 = iStart[0]>=iEnd[0] ? iStart[0] : iEnd[0];

                for(int j=iStart[1]-iW; j<=iEnd[1]+iW; j++){
                    ArrayList<EnvironmentPoint> l = profile.get(i);
                    EnvironmentPoint ep = l.get(j);
                    ep.open = true;
                    ep.start = entry;
                    ep.fixed = fixed;
                    ep.wrong = wrong;
                    ep.right = right;
                    ep.c = ep.c_m1    = concSource1;
                }
            }
        }

    }

    public static ChemicalEnvironment FoolsMazeEnvironment(){

        ArrayList<double[]> centrePoints = new ArrayList<double[]>();
        ArrayList<double[]> newCentrePoints = new ArrayList<double[]>();
        ChemicalEnvironment ce = new ChemicalEnvironment();

        // Make space in the environment.
        int LengthScale = (int) Math.ceil((2 * (L + 2.0 * padding)) / grain);

        for (int i = 0; i <= 1.5*LengthScale+W; i++) {

            ArrayList<EnvironmentPoint> alENV = new ArrayList<EnvironmentPoint>();

            for (int j = 0; j <=1.2*LengthScale+W; j++) {

                EnvironmentPoint ep = new EnvironmentPoint(i, j);
                ep.open = false;
                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }

        double cX = grain * (ce.profile.size() / 3.5);
        ce.DrawLine(new double[]{cX, W + padding+120}, new double[]{cX, W + padding + 220}, W, false, false, false, false); // Trunk
        ce.DrawLine(new double[]{cX, W + padding + 120}, new double[]{cX, W + padding + 120}, 30, true, false, false, false); // ORIGIN

        centrePoints.add(new double[]{cX, W + padding + 220});

        double[] ncp;
        double[] nep;

        double lW = 1.05*L;
        double lH = 2*L;

        double[] ctr = centrePoints.get(0);

        ncp = new double[]{ctr[0] - lW, ctr[1]};
        nep = new double[]{ctr[0] - lW, ctr[1] + 1.5*lH};

        ce.DrawLine(ncp, ctr, W, false, false, false, false);
        ce.DrawLine(ncp, nep, W, false, false, false, true);
        ce.DrawLine(nep, new double[]{nep[0], nep[1]+2}, W, false, true, false, false);

        ncp = new double[]{ctr[0] + lW, ctr[1]};
        nep = new double[]{ctr[0] + lW, ctr[1] + lH/8};

        ce.DrawLine(ctr, ncp, W, false, false, false, false);
        ce.DrawLine(ncp, nep, W, false, false, true, false);

        centrePoints.clear();
        centrePoints.add(nep);

        ce.DrawLine(ncp, new double[]{nep[0], nep[1] + 1.5 * L}, W, false, false, true, false);
        centrePoints.add(new double[]{nep[0], nep[1] + 1.5*L});

        for (int i = 2; i <= k; i++) {

            lW = Math.pow(0.8, 3*(i-2))*L;

            lH = L - lW;

            for (int j = 0; j < centrePoints.size(); j++) {
                ctr = centrePoints.get(j);

                nep = new double[]{ctr[0] + lW, ctr[1] + lH};

                ncp = new double[]{ctr[0] + lW, ctr[1]};

                newCentrePoints.add(nep);

                ce.DrawLine(ncp, nep, W, false, false, false, true);
                ce.DrawLine(ctr, ncp, W, false, false, false, false);

                ncp = new double[]{ctr[0] - lW, ctr[1]};
                nep = new double[]{ctr[0] - lW, ctr[1] + lH};
                newCentrePoints.add(nep);

                ce.DrawLine(ncp, nep, W, false, false, true, false);
                ce.DrawLine(ncp, ctr, W, false, false, false, false);

            }

            centrePoints.clear();
            centrePoints.addAll(newCentrePoints);
            newCentrePoints.clear();
        }

        ce.SetUpDiffusion();
        return ce;
    }

    public static ChemicalEnvironment TreeMazeEnvironment() {

        ArrayList<double[]> centrePoints = new ArrayList<double[]>();
        ArrayList<double[]> newCentrePoints = new ArrayList<double[]>();
        ChemicalEnvironment ce = new ChemicalEnvironment();

        // Make space in the environment.
        int LengthScale = (int) Math.ceil((2 * (L + 2.0 * padding)) / grain);
        if (asym >= 1.0) LengthScale *= (0.9 * asym);

        // Make the environment and fill a grid with closed EnvironmentPoint objects


        for (int i = 0; i <= 40+2*W; i++) {

            ArrayList<EnvironmentPoint> alENV = new ArrayList<EnvironmentPoint>();

            for (int j = 0; j < ((2*W+200+trunk)/grain+(2.1*L*Math.max(1, 1*Math.pow(2,skew)))/grain); j++) {

                EnvironmentPoint ep = new EnvironmentPoint(i, j);
                ep.open = false;
                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }

        // Draw the environment by opening EnvironmentPoints;

        double cX = grain * (ce.profile.size() / 2.0);

        //ce.DrawLine(new double[]{cX, W+padding+200}, new double[]{cX, W+padding+200}, 200, false, false);
        //ce.DrawLine(new double[]{cX, W + padding + grain}, new double[]{cX, W + padding + 100}, 50, true, false); // ORIGIN

        ce.DrawLine(new double[]{cX, W + padding+20}, new double[]{cX, W + padding + 20 + 30 + trunk}, W, false, false, false, false); // Trunk
        ce.DrawLine(new double[]{cX, W + padding + 20}, new double[]{cX, W + padding + 20}, 30, true, false, false, false); // ORIGIN

        centrePoints.add(new double[]{cX, W + padding + 20 + 30 + trunk});

        double[] ncp;
        double[] nep;

        for (int i = 1; i <= k; i++) {

            //double lW = Math.pow(1.8, (1.*(k-i)))*L/2;
            double lW = 120;
            double lH = 2*L - lW;


            for (int j = 0; j < centrePoints.size(); j++) {
                double[] ctr = centrePoints.get(j);

                nep = new double[]{ctr[0] + lW, ctr[1] + lH};
                ncp = new double[]{ctr[0] + lW, ctr[1]};

                newCentrePoints.add(nep);

                ce.DrawLine(ncp, nep, W, false, false, false, true);
                ce.DrawLine(ctr, ncp, W, false, false, false, false);

                ncp = new double[]{ctr[0] - lW, ctr[1]};
                nep = new double[]{ctr[0] - lW, ctr[1] + lH*Math.pow(2,skew)};
                newCentrePoints.add(nep);

                ce.DrawLine(ncp, nep, W, false, false, true, false);
                ce.DrawLine(ncp, ctr, W, false, false, false, false);

            }

            centrePoints.clear();
            centrePoints.addAll(newCentrePoints);
            newCentrePoints.clear();

            if(i==k){
                //Collections.shuffle(centrePoints);
                ce.DrawLine(new double[]{centrePoints.get(0)[0], centrePoints.get(0)[1]+2*W-1}, new double[]{centrePoints.get(0)[0], centrePoints.get(0)[1]+2*W}, W, false, true, false, false);

            }
        }

        ce.SetUpDiffusion();
        return  ce;
    }

    public static ChemicalEnvironment GridMazeEnvironment() {

        ChemicalEnvironment ce = new ChemicalEnvironment();
        int LengthScale = (int) Math.ceil((5.0*L+2.0*padding)/grain);

        for(int i = 0; i <=LengthScale/4; i++){

            ArrayList<EnvironmentPoint>  alENV = new ArrayList<EnvironmentPoint>();

            for(int j = 0; j < LengthScale/2; j++){

                EnvironmentPoint ep = new EnvironmentPoint(i,j);

                ep.open = false;

                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }

        double cX = grain*(ce.profile.size()/2.0);


        double baseline = W+padding;

        ce.DrawLine(new double[] {cX, baseline+L/3}, new double[] {cX, baseline+L/2}, W, false, false, false, false); // Trunk
        ce.DrawLine(new double[]{cX, baseline + L / 3}, new double[]{cX, baseline + 2 + L / 3}, 3 * W, true, false, false, false); // ORIGIN

        for(int i=0; i<=k; i++){

            ce.DrawLine(new double[]{cX-L/4,L/2+i*L/k + baseline},new double[]{cX+L/4,L/2+i*L/k + baseline}, W, false, false, false, false);
        }

        for(int i=0; i<k; i+=2){

            ce.DrawLine(new double[]{cX - L / 4 + (i * 2*L / (3*k)), L / 2 + baseline}, new double[]{cX - L / 4 + (i * 2*L / (3*k)), 3 * L / 2 + baseline}, W,false,false, false, false);
        }

        ce.DrawLine(new double[] {cX, baseline + 3*L/2}, new double[] {cX, baseline+ 4*L/2}, W, false, false, false, false); // Goal
        ce.DrawLine(new double[] {cX, baseline + 4 * L/2 - W}, new double[]{cX, baseline+ 4 * L/2+2},  W, false, true, false, false); // Goal

        ce.SetUpDiffusion();
        return ce;
    }

    public void SetUpDiffusion(){

        freepoints = new ArrayList<EnvironmentPoint>();

        for(int i=0; i<profile.size(); i++){
            for(int j = 0; j<profile.get(i).size(); j++){

                EnvironmentPoint ep = profile.get(i).get(j);

                if(!ep.open) continue;;

                freepoints.add(ep);

                int ixp1, ixm1, iyp1, iym1;
                ixp1 = ixm1 = i;
                iyp1 = iym1 = j;

                if(i>0)  ixm1-=1;
                if(i<profile.size()-1)  ixp1+=1;
                if(j>0)  iym1-=1;
                if(j<profile.get(i).size()-1)  iyp1+=1;


                ep.xp1 = profile.get(ixp1).get(j);
                ep.xm1 = profile.get(ixm1).get(j);
                ep.yp1 = profile.get(i).get(iyp1);
                ep.ym1 = profile.get(i).get(iym1);

                if(!ep.xp1.open) ep.xp1 = ep;
                if(!ep.xm1.open) ep.xm1 = ep;
                if(!ep.yp1.open) ep.yp1 = ep;
                if(!ep.ym1.open) ep.ym1 = ep;
            }
        }

        /*kernel = new Kernel(){
            @Override public void run() {
                int gid = getGlobalId();

                EnvironmentPoint ep = freepoints.get(gid);

                if(ep.fixed) return;

                double cpx = ep.xp1.open ? ep.xp1.c_m1 : ep.c_m1;
                double cpy = ep.yp1.open ? ep.yp1.c_m1 : ep.c_m1;
                double cmx = ep.xm1.open ? ep.xm1.c_m1 : ep.c_m1;
                double cmy = ep.ym1.open ? ep.ym1.c_m1 : ep.c_m1;

                ep.c += (MigrationSimulation.DiffC*MelaMigration.dt/(grain*grain))*(cpx+cpy+cmx+cmy-4.0*ep.c_m1);
            }
        }; */
    }

    public boolean GetIsOpen(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.open;
    }

    public boolean GetIsFixed(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.fixed;
    }

    public boolean GetIsStart(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.start;
    }

    public boolean GetIsLost(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.wrong;
    }

    public boolean GetIsCorrect(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.right;
    }

    private double[] ConvertCoordinatesToDouble(double cx, double cy){

        double x = cx/grain;
        double y = cy/grain;

        x = Math.max(1,Math.min(x,profile.size()-1)-1);
        y = Math.max(1,Math.min(y,profile.get(0).size()-1)-1);

        return new double[]{x,y};
    }

    public double litresFromUMVolume (double length, double width, double depth){

        double umCubedVolume= width*depth*length;
        double cmCubedVolume=umCubedVolume/(Math.pow(10,12));
        double litres= cmCubedVolume/1000;

        return litres;

    }

    public void degradeFromSourceOrSink (){

        //all following in micro-moles

        double sMaxMolesPerMin = litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], bridgeDepth) * sMax;
        double kMMoles = litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], bridgeDepth) * kM;


        double sourceMoles = concSource1 *litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);
        double sinkMoles = concSink1 *litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);


        double volumeSourceUM=AgentBasedSimulation.dimensions[0]* AgentBasedSimulation.dimensions[1]* resevoirDepth;
        double volumeCellUM=(4.0/3)*Math.PI*Math.pow(Cell.smallCellRadius,3);
        double molesPerCellSource=volumeCellUM/volumeSourceUM*sourceMoles;
        double molesPerCellSink=volumeCellUM/volumeSourceUM*sinkMoles;


        double degradationRateSourceInMolesPerMin=noCellsInSource*sMaxMolesPerMin*(molesPerCellSource/(molesPerCellSource+kMMoles));
        double degradationRateSinkInMolesPerMin=noCellsInSource*sMaxMolesPerMin*(molesPerCellSink/(molesPerCellSink+kMMoles));

        double amountDegradedInMolesSource=degradationRateSourceInMolesPerMin*AgentBasedSimulation.dt;
        double amountDegradedInMolesSink=degradationRateSinkInMolesPerMin*AgentBasedSimulation.dt;

      //  System.out.println(sourceMoles+"    "+amountDegradedInMolesSource);

            concSource1 -= amountDegradedInMolesSource / litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);
            if (degradationInSourceGivesInhibitor) concSource2 += amountDegradedInMolesSource / litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);

            concSink1 -= amountDegradedInMolesSink / litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);
            if (degradationInSourceGivesInhibitor) concSink2 += amountDegradedInMolesSink / litresFromUMVolume(AgentBasedSimulation.dimensions[0], AgentBasedSimulation.dimensions[1], resevoirDepth);

   //     System.out.println(concSource1+"    "+concSource2+"    "+concSink1+"    "+concSink2+"    "+MigrationSimulation.Ttotal);


    }



    public void diffuseIntoSourceOrSink(){              //note that mechanics of the code already deal with diffusion out of a source or sink (down the gradient) but dont deal with diffusion in to the source or sink!


        double sumAttr1=0;
        double sumInhib1=0;

        double sumAttrMaxMinus=0;
        double sumInhibMaxMinus=0;

        int iMax = ((int) (AgentBasedSimulation.dimensions[0] / grain));

        for (int i=0;i<freepoints.size();i++){

            if (freepoints.get(i).x / ChemicalEnvironment.grain == 1){
                sumAttr1+=freepoints.get(i).c;
                sumInhib1+=freepoints.get(i).c2;
            }
            if (freepoints.get(i).x / ChemicalEnvironment.grain == iMax-2){
                sumAttrMaxMinus+=freepoints.get(i).c;
                sumInhibMaxMinus+=freepoints.get(i).c2;
            }
        }

        double avAttr1=sumAttr1/(AgentBasedSimulation.dimensions[1]/ChemicalEnvironment.grain);
        double avInhib1=sumInhib1/(AgentBasedSimulation.dimensions[1]/ChemicalEnvironment.grain);
        double avAttrMaxMinus=sumAttrMaxMinus/(AgentBasedSimulation.dimensions[1]/ChemicalEnvironment.grain);
        double avInhibMaxMinus=sumInhibMaxMinus/(AgentBasedSimulation.dimensions[1]/ChemicalEnvironment.grain);

        double chemGradAttrStart=((avAttr1-freepoints.get(0).c)/ChemicalEnvironment.grain)/Math.pow(10,15);                                   //divison by 1x10^15 puts uMoles per litre into uMoles per um^3
        double chemGradAttrEnd= ((freepoints.get(freepoints.size()-1).c-avAttrMaxMinus)/ChemicalEnvironment.grain)/Math.pow(10,15);
        double chemGradInhibStart=((avInhib1-freepoints.get(0).c2)/ChemicalEnvironment.grain)/Math.pow(10,15);
        double chemGradInhibEnd= ((freepoints.get(freepoints.size()-1).c2-avInhibMaxMinus)/ChemicalEnvironment.grain)/Math.pow(10,15);

        double fluxAttrStart=DiffC*chemGradAttrStart;
        double fluxInhibStart=DiffC*chemGradInhibStart;
        double fluxAttrEnd=DiffC*chemGradAttrEnd;
        double fluxInhibEnd=DiffC*chemGradInhibEnd;

        if (fluxAttrStart>0){
            double molesAttrIntoSink=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxAttrStart);
            concSink1 +=molesAttrIntoSink/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
 //           System.out.println("uM Attr into sink: "+molesAttrIntoSink/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }else{
            double molesAttrOutOfSink=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxAttrStart);
            concSink1 -=molesAttrOutOfSink/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
 //           System.out.println("uM Attr out of sink: "+molesAttrOutOfSink/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }


        if (fluxInhibStart>0){
            double molesInhibIntoSink=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxInhibStart);
            concSink2 +=molesInhibIntoSink/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
  //          System.out.println("uM Inhib into sink: "+molesInhibIntoSink/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }else {
            double molesInhibOutOfSink=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxInhibStart);
            concSink2 -=molesInhibOutOfSink/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
 //           System.out.println("uM Inhib out of sink: "+ molesInhibOutOfSink/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }


        if (fluxAttrEnd<0){
            double molesAttrIntoSource=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxAttrEnd);
            concSource1 +=molesAttrIntoSource/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
  //          System.out.println("uM Attr into source: "+molesAttrIntoSource/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }else{
            double molesAttrOutOfSource=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxAttrEnd);
            concSource1 -=molesAttrOutOfSource/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
  //          System.out.println("uM Attr out of source: "+molesAttrOutOfSource/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }


        if (fluxInhibEnd<0){
            double molesInhibIntoSource=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxInhibEnd);
            concSource2 +=molesInhibIntoSource/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
 //          System.out.println("uM inhib into source: "+molesInhibIntoSource/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }else{
            double molesInhibOutOfSource=bridgeDepth*AgentBasedSimulation.dimensions[1]*AgentBasedSimulation.dt*Math.abs(fluxInhibEnd);
            concSource2 -=molesInhibOutOfSource/litresFromUMVolume(AgentBasedSimulation.dimensions[0],AgentBasedSimulation.dimensions[1],resevoirDepth);
 //           System.out.println("uM inhib out of source: "+molesInhibOutOfSource/litresFromUMVolume(resevoirLength,resevoirWidth,resevoirDepth));
        }

    }



    public void Diffuse() {
        //Explicit method

        double dx = grain;
        double dt = AgentBasedSimulation.dt;

        for(EnvironmentPoint ep : freepoints) ep.c_m2 = ep.c_m1;
        for(EnvironmentPoint ep : freepoints) ep.c_m1 = ep.c;

        for(EnvironmentPoint ep : freepoints) ep.c2_m2 = ep.c2_m1;
        for(EnvironmentPoint ep : freepoints) ep.c2_m1 = ep.c2;

        for(EnvironmentPoint ep : freepoints) ep.c3_m2 = ep.c3_m1;
        for(EnvironmentPoint ep : freepoints) ep.c3_m1 = ep.c3;


        double k = 2* DiffC*dt/(dx*dx);


        if (sourceSinkDegradation){
            degradeFromSourceOrSink();
        }
        if (fluxIntoOutOfSourceAndSink){
            diffuseIntoSourceOrSink();
        }




        //freepoints.stream().parallel().forEach(ep -> {

        for (EnvironmentPoint ep : freepoints) {

            if (sourceSinkDegradation) {

                int iMax = ((int) (AgentBasedSimulation.dimensions[0] / grain));

                if (ep.x/ChemicalEnvironment.grain==0){
                    ep.c= concSink1;
                    ep.c2= concSink2;
                }

                if (ep.x / ChemicalEnvironment.grain == iMax - 1) {
                    ep.c = concSource1;
                    ep.c2 = concSource2;
                }

            }


            double cpx, cpy, cmx, cmy;


            if (ep.fixed) continue;


            else {

                cpx = ep.xp1.c_m1;
                cpy = ep.yp1.c_m1;
                cmx = ep.xm1.c_m1;
                cmy = ep.ym1.c_m1;

                // Eulerian FTCS
                //ep.c += (MigrationSimulation.DiffC * MelaMigration.dt / (grain * grain)) * (cpx + cpy + cmx + cmy - 4.0 * ep.c_m1);

                //DuFort Frankel
                ep.c = ((1.0 - 2.0 * k) / (1.0 + 2.0 * k)) * ep.c_m2 + (k / (1.0 + 2.0 * k)) * (cpx + cmx + cpy + cmy);
            }

        }

            for(EnvironmentPoint ep : freepoints){
                if(ep.fixed) continue;
                double cpx, cpy, cmx, cmy;

                cpx = ep.xp1.c2_m1;
                cpy = ep.yp1.c2_m1;
                cmx = ep.xm1.c2_m1;
                cmy = ep.ym1.c2_m1;

                //ep.c2 += (MigrationSimulation.DiffC*MelaMigration.dt/(grain*grain))*(cpx+cpy+cmx+cmy-4.0*ep.c2_m1);

                ep.c2 = ((1.0-2.0*k)/(1.0+2.0*k))*ep.c2_m2 + (k/(1.0+2.0*k))*(cpx + cmx + cpy + cmy);
            }

        for(EnvironmentPoint ep : freepoints){
            if(ep.fixed) continue;
            double cpx, cpy, cmx, cmy;

            cpx = ep.xp1.c3_m1;
            cpy = ep.yp1.c3_m1;
            cmx = ep.xm1.c3_m1;
            cmy = ep.ym1.c3_m1;

            //ep.c2 += (MigrationSimulation.DiffC*MelaMigration.dt/(grain*grain))*(cpx+cpy+cmx+cmy-4.0*ep.c2_m1);

            ep.c3 = ((1.0-2.0*k)/(1.0+2.0*k))*ep.c3_m2 + (k/(1.0+2.0*k))*(cpx + cmx + cpy + cmy);
        }

    }

    ArrayList<EnvironmentPoint> GetStartingPoints(){

        ArrayList<EnvironmentPoint> points = new ArrayList<EnvironmentPoint>();

        for(ArrayList<EnvironmentPoint> eps : profile){
           for(EnvironmentPoint ep : eps){
               if(ep.start) points.add(ep);
           }
        }
        return points;
    }
}
