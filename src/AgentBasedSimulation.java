import java.io.File;
import java.io.IOException;
import java.util.Random;


/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:41
 * To change this template use File | Settings | File Templates.
 */
public class AgentBasedSimulation
{

    public static boolean introduceNewCells=false;
    public static boolean introduceNewClusters=false;
    public static double whenToIntroduceNewCellsInMins =1;
    public static double whenToIntroduceNewClustersInMins =0.5;
    static int pop =100;                    // Initial population
    static double T=60;                // Days, Hours, Minutes
    public static double dt = 0.0005;         // 0.05~3.0s
    public static double alpha = 0.6;
    public static double outInt = 3;
    public static String directory =System.getProperty("user.home")+"/UASimulationFOrVariableBreakdown/4/";
    public static double  dimensions[] = {600, 500, 1000, 500}; // um




    public static boolean visualise =true;
    public static boolean record = true;
    public static boolean exdeg  = true;
    public static int     yPos = 40;
    public static boolean abs = true;
    public static boolean pinned = true;
    public static boolean ctc = false;
    public static boolean competivite = false;
    public static double rdt = Math.sqrt(dt);   // Sqrt of dt for brownian motion

    public static boolean introduceCellsNow=false;              //keep these false as they will switch to true automatically when new cells to be introduced/mechanics to occur
    public static boolean introduceClustersNow=false;








            //"C:\\Users\\adowdell\\Google Drive\\PhD\\Beatson Work\\Science\\UASimulationFOrVariableBreakdown\\4\\";

    public static boolean finished = true;

    public Thread visualiser;
    //public ExecutorService threadpool;
    //ArrayList<Callable<Object>> tasks;

    public MigrationSimulation RWS;

    private static void RandomiseParameters() {

        Random r = new Random();

        ChemicalEnvironment.DiffC = 1000+2*Math.pow(10,(r.nextDouble()*4d));                                 // r.nextDouble() gives random number between 0 and 1
        ChemicalEnvironment.concSource1 =Cell.k1 *(Math.pow(10,(r.nextInt(4))));

        ChemicalEnvironment.sMax = Cell.k1 *(Math.pow(10,(r.nextDouble()*4d)));
        ChemicalEnvironment.kM = Cell.k1 *(Math.pow(10,(r.nextDouble()*4d)));

        ChemicalEnvironment.resevoirDepth=ChemicalEnvironment.bridgeDepth*Math.pow(10,(r.nextDouble()*2d)+2d);
        ChemicalEnvironment.noCellsInSource=(int)(AgentBasedSimulation.pop*(Math.pow(10,(r.nextDouble()*4d))));

        Cell.fractionCluster=0.1+0.8*r.nextDouble();

        whenToIntroduceNewCellsInMins=r.nextDouble()*T;

        double rand=r.nextDouble();

        if (rand<0.5){
            Cell.clusterDegradesFromSurface=true;
        }else {
            Cell.clusterDegradesFromSurface=false;
        }

    }


    private AgentBasedSimulation(boolean skipguis){

        if(record){
            boolean bDir = new File(directory).mkdirs();
        }
        if(visualise && !skipguis){
           makeSimGUIs();
        }
        else{
            RWS = new MigrationSimulation(false,false,ctc,abs,alpha, Cell.speed,dt,ChemicalEnvironment.grain, ChemicalEnvironment.DiffC, Cell.k1, ChemicalEnvironment.kM, ChemicalEnvironment.sMax);
        }
        Thread th1 = new Thread() {
            public void run() {
                long startTime = 0;
                long endTime   = 0;

                int j = 0;
                int k=0;
                int l=0;
                for (int i = 0; i <= T / dt; i++) {
                    if(i==0)    startTime = System.currentTimeMillis();
                    else if(i%5000 == 0) {
                        endTime = System.currentTimeMillis() - startTime;
                        System.out.println("TIME 1:: "+endTime);
                        startTime = System.currentTimeMillis();
                    }
                    RWS.Ttotal += dt;

                    if(RWS.complete) break;
                    if (RWS.paused) {
                        RWS.Ttotal-=dt;
                        i--;
                    }
                    else{
                        RWS.Iterate();
                        if (i % 1000 == 1 && visualise) RWS.controlPanel.t_time.setText(Double.toString(dt * i / 60.0));

                        if (j * dt >= outInt) {
                            j = 0;
                            if (visualise) RWS.SnapImage("con");
                            RWS.WriteCellData();
                            RWS.WriteEnvironmentData();
                        }
                        j++;
                        if (k*dt>=whenToIntroduceNewCellsInMins){
                            k=0;
                            if (introduceNewCells) introduceCellsNow=true;
                        }
                        k++;
                        if (l*dt>=whenToIntroduceNewClustersInMins){
                            l=0;
                            if (introduceNewClusters) introduceClustersNow=true;
                        }
                        l++;
                    }
                }

                if(record) {


                    RWS.WriteCellData();
                    RWS.WriteEnvironmentData();

                    RWS.WriteShortRecord();

                    try {
                        RWS.bw.flush();
                        RWS.bw.close();
                        RWS.fw.close();
                        RWS.bw4.flush();
                        RWS.bw4.close();
                        RWS.fw4.close();
                        RWS.bw3.flush();
                        RWS.bw3.close();
                        RWS.fw3.close();
                    } catch (IOException e) {
                    }
                }
                finished = true;
                System.exit(1);
            }
            //if(visualise) while(true)    try{sleep(100);} catch(Exception e){}
        };
        th1.run();
    }

    public static void main(String[] args){
        if(args.length==0) new AgentBasedSimulation(false);
        if(args.length==1){
            if(args[0].toString().equals("R")) {
                RandomiseParameters();
                //visualise = false;
                System.out.println("Randomised");
                new AgentBasedSimulation(true);
            }
        }
        else if(args.length%2==0){
            parseArgs(args);
            visualise = false;
            new AgentBasedSimulation(false);
        }
        else System.exit(-1);
    }

    private void makeSimGUIs(){

        SimGUIPanel s1 = new SimGUIPanel();

        s1.create();

        RWS = s1.SetupSimulation("/1.txt");

    }

    private static void parseArgs(String[] args){

        for(int i=0; i<args.length; i+=2){
            try{
                if(args[i].equals("p"))             pop = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("D"))        ChemicalEnvironment.DiffC = Double.parseDouble(args[i+1]);
                else if(args[i].equals("k1"))       Cell.k1 = Double.parseDouble(args[i+1]);
                else if(args[i].equals("kM"))       ChemicalEnvironment.kM    = Double.parseDouble(args[i+1]);
                else if(args[i].equals("sMax"))     ChemicalEnvironment.sMax  = Double.parseDouble(args[i+1]);
                else if(args[i].equals("a"))        abs = Boolean.parseBoolean(args[i+1]);
                else if(args[i].equals("c"))        ctc = Boolean.parseBoolean(args[i+1]);
                else if(args[i].equals("in"))       MigrationSimulation.sMazePicture = args[i+1];
                else if(args[i].equals("out"))      directory = args[i+1];
                else if(args[i].equals("alpha"))    alpha = Double.parseDouble(args[i+1]);
                else if(args[i].equals("k"))        ChemicalEnvironment.k = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("L"))        ChemicalEnvironment.L = Double.parseDouble(args[i + 1]);
                else if(args[i].equals("skew"))     ChemicalEnvironment.skew = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("asym"))     ChemicalEnvironment.asym = Double.parseDouble(args[i + 1]);
            }
            catch(NumberFormatException e){}
        }

    }
}



