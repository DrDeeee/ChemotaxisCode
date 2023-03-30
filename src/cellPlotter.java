import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Ellipse2D;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 09:50
 * To change this template use File | Settings | File Templates.
 */
public class CellPlotter extends JPanel {

    static int border = 2+ (int)Math.ceil(ChemicalEnvironment.grain/2.0);  //pixel border
    JFrame frame;
    int    upCalls;
    CellPlotter THIS = this;
    List<Cell> cells;
    ChemicalEnvironment environment;
    MigrationSimulation ms;

    private Point clickPoint;
    private Point dragPoint;
    private Rectangle selectionBounds;
    private Point     ovalSelection;
    private int ovalWidth = 0;
    private int Ovalheight = 0;
    private boolean selectionIsOval = false;
    private boolean shiftDown = false;
    private Color[]  colourSchemeCNC = {Color.WHITE,  Color.GREEN, Color.YELLOW, Color.ORANGE, Color.RED};
    private double[] colourPointsCNC = {0,            5* Cell.k1,     10* Cell.k1,   25* Cell.k1,   50* Cell.k1};

    private Color[]  colourSchemeOCC = {Color.WHITE,  Color.GREEN, Color.YELLOW, Color.ORANGE, Color.RED};
    private double[] colourPointsOCC = {0.0,          0.4,         0.6,          0.8,          0.95};

    private Color[]  colourSchemeOCD = {Color.RED, Color.ORANGE, Color.WHITE, Color.WHITE, Color.GREEN, Color.BLUE};
    private double[] colourPointsOCD = {-0.02,     -0.005,        0.0,        0.002,     0.006,        0.015};

    private Color[]  colourSchemeOCA = {Color.WHITE, Color.RED, Color.ORANGE,        Color.GREEN,    Color.BLUE};
    private double[] colourPointsOCA = {0.00000002,   0.0000002,     0.000002,       0.00002,        0.0002};


    private Color[]  colourSchemeOC2 = {new Color(200,0,220), Color.BLUE, Color.WHITE, Color.WHITE, Color.RED, Color.YELLOW, Color.ORANGE, Color.GREEN};
    private double[] colourPointsOC2 = {-0.01,                -0.003,     -0.00001,        0.002,       0.004,     0.006,        0.012,        0.015};



    private Color [] colourSchemePosChem = {Color.BLACK, Color.RED, Color.ORANGE, Color.YELLOW};
    private double [] colorPointsPosChem = {0.0, 0.003, 0.006, 0.01};

    private Color [] colourSchemeNegChem = {Color.WHITE, Color.CYAN, Color.BLUE, Color.BLACK};
    private double [] colorPointsNegChem = {-0.01,-0.006, -0.003, 0.0};

    private Color [] colourScheme2DBridge={Color.WHITE, Color.GREEN, Color.YELLOW,Color.ORANGE, Color.RED};
    private double [] colorPoints2DBridge ={0,0.001,0.005, 0.01, 0.05};



    private static BasicStroke pseudoStroke;

    public static enum GRAPHTYPE {CONCENTRATION, CONC2, OCCUPANCY, OCCUPANCYGRADIENT, OC2, };
    public static GRAPHTYPE gt = GRAPHTYPE.CONCENTRATION;

    public CellPlotter(List<Cell> cells, ChemicalEnvironment environment, MigrationSimulation ms, JFrame frame){
        this.ms = ms;
        this.environment = environment;
        this.frame = frame;
        this.cells = cells;
        this.setFocusable(true);

        selectionBounds = new Rectangle(1950, 1375, 250, 250);

        this.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent keyEvent) {
                if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT) shiftDown = true;

            }

            @Override
            public void keyReleased(KeyEvent keyEvent) {
                if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT) {
                    shiftDown = false;

                }
            }
        });

        MouseAdapter handler = new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e){
                selectionBounds = null;
                ovalSelection = convertPointToSimulationSpace(e.getPoint());
                getParent().repaint();
            }

            @Override
            public void mousePressed(MouseEvent e){
                selectionIsOval = shiftDown;
                clickPoint = convertPointToSimulationSpace(e.getPoint());
                selectionBounds = null;
                ovalSelection = null;
            }

            @Override
            public void mouseReleased(MouseEvent e){
                clickPoint = null;
                if(selectionBounds!=null); //selectAllInBounds(selectionBounds);
                //selectionBounds = null;
                ovalSelection = null;
                getParent().repaint();
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                Point dragPoint = convertPointToSimulationSpace(e.getPoint());
                int x = Math.min(clickPoint.x, dragPoint.x);
                int y = Math.min(clickPoint.y, dragPoint.y);
                int width = Math.max(clickPoint.x - dragPoint.x, dragPoint.x - clickPoint.x);
                int height = Math.max(clickPoint.y - dragPoint.y, dragPoint.y - clickPoint.y);
                selectionBounds = new Rectangle(x, y, width, height);
                getParent().getParent().repaint();
            }
        };

        this.addMouseListener(handler);
        this.addMouseMotionListener(handler);
    }

    // This method is called whenever the contents needs to be painted
    public synchronized void paintComponent(Graphics g) {
        getParent().repaint();

        int width = this.frame.getContentPane().getWidth();
        int height = this.frame.getContentPane().getHeight();

        double hScale = (width-2.0*border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-2.0*border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Cell c;
        Graphics2D g2d = (Graphics2D)g;

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.BLACK);
        //g2d.setStroke(new BasicStroke(2F));

        //for(int i=1; i<=24; i++){
        //    g2d.drawLine((int) ((1000.0*i/ms.xMax)*hScale*ms.xMax), 0, (int) ((1000.0*i/ms.xMax)*hScale*ms.xMax), 200);
        //}

        g2d.setStroke(new BasicStroke(1F));
        for(int i = 1; i<environment.profile.size()-1; i++){
            for(int j = 1; j<environment.profile.get(i).size()-1; j++){

                EnvironmentPoint ep  = environment.profile.get(i).get(j);
                EnvironmentPoint ep2 = environment.profile.get(i+1).get(j);
                EnvironmentPoint ep2Y = environment.profile.get(i).get(j+1);
                EnvironmentPoint epb = environment.profile.get(i-1).get(j);

                double gradX    = ((Cell.a*ep2.c/Cell.k1)+(Cell.b*ep2.c2/Cell.k2)+(Cell.c*ep2.c3/Cell.k3))/(1+((ep2.c/Cell.k1)+(ep2.c2/Cell.k2)+(ep2.c3/Cell.k3))) - ((Cell.a*ep.c/Cell.k1)+(Cell.b*ep.c2/Cell.k2)+(Cell.c*ep.c3/Cell.k3))/(1+((ep.c/Cell.k1)+(ep.c2/Cell.k2)+(ep.c3/Cell.k3)));
                //System.out.println(i+"    "+grad);
                double gradY   = ((Cell.a*ep2Y.c/Cell.k1)+(Cell.b*ep2Y.c2/Cell.k2)+(Cell.c*ep2Y.c3/Cell.k3))/(1+((ep2Y.c/Cell.k1)+(ep2Y.c2/Cell.k2)+(ep2Y.c3/Cell.k3))) - ((Cell.a*ep.c/Cell.k1)+(Cell.b*ep.c2/Cell.k2)+(Cell.c*ep.c3/Cell.k3))/(1+((ep.c/Cell.k1)+(ep.c2/Cell.k2)+(ep.c3/Cell.k3)));
                double gradb = (ep.c/(ep.c+ Cell.k1))-(epb.c/(epb.c+ Cell.k1));
                double attrGrad= ep2.c-ep.c;
                double occupDiffX=(2*Cell.smallCellRadius/ChemicalEnvironment.grain)*gradX;
                double occupDiffY=(2*Cell.smallCellRadius/ChemicalEnvironment.grain)*gradY;




                Color clr = Color.BLACK;
                if(ep.fixed) clr = Color.BLUE;
                else if(ep.open) {
                    if(     gt==GRAPHTYPE.OCCUPANCYGRADIENT)   {
                        if (gradX>0){
                            clr = ColorBlend(colourSchemePosChem, colorPointsPosChem, occupDiffX);
                        }else {
                            clr = ColorBlend(colourSchemeNegChem, colorPointsNegChem, occupDiffX);
                        }
                    }
                    else if(gt==GRAPHTYPE.OCCUPANCY)            clr = ColorBlend(colourSchemeOCC, colourPointsOCC, ep.c/(ep.c+ Cell.k1));
                    else if(gt==GRAPHTYPE.OC2)                  {
                        if (gradY>0){
                            clr = ColorBlend(colourSchemeNegChem, colorPointsNegChem, -occupDiffY);
                        }else {
                            clr = ColorBlend(colourSchemePosChem, colorPointsPosChem, -occupDiffY);
                        }
                    }
                    else if(gt==GRAPHTYPE.CONCENTRATION)        clr = ColorBlend(colourSchemeCNC, colourPointsCNC, ep.c);
                    else if(gt==GRAPHTYPE.CONC2)                clr = ColorBlend(colourSchemeCNC, colourPointsCNC, ep.c2*(Cell.k1 / Cell.k2));
                }
                g2d.setColor(clr);
                g2d.fillRect(border + (int) ((hScale*ChemicalEnvironment.grain*(i))-3),
                             border+(int) ((vScale*ChemicalEnvironment.grain*(j))-3),
                             (int) ((hScale*ChemicalEnvironment.grain))+6,
                             (int) ((vScale*ChemicalEnvironment.grain))+6);
            }
        }

        g2d.setStroke(new BasicStroke(1F));

        float cm = 18;

        int off = -1;

        for(int i = cells.size()-1; i>=0; i--) {
            c = cells.get(i);

            if (c.minrad == Cell.smallCellRadius) {

                Stroke outerStroke = new BasicStroke((cm/1f)*(float)(hScale*vScale),
                        BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

                Stroke innerStroke = new BasicStroke((cm/1.35f)*(float)(hScale*vScale),
                        BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

                g2d.setColor(Color.WHITE);
                //g2d.drawOval(border + (int) (hScale * (c.x() - (-2 + c.minrad))), border + (int) (vScale * (c.y() - (-2 + c.minrad))), (int) Math.max(3.0, 2.0 * hScale * c.minrad - 4), (int) Math.max(3.0, 2.0 * vScale * c.minrad - 4));

                if (c.original) g2d.setColor(new Color(60, 120, 200, 120));
                //else g2d.setColor(new Color(50,200,150,120));
                //g2d.fillOval(border + (int) (hScale * (c.x() - (-2 + c.minrad))), border + (int) (vScale * (c.y() - (-2 + c.minrad))), (int) Math.max(3.0, 2.0 * hScale * c.minrad - 4), (int) Math.max(3.0, 2.0 * vScale * c.minrad - 4));

                g2d.setColor(new Color(200, 200, 200, 200));
                g2d.setStroke(outerStroke);


                //g2d.drawLine((int) (hScale*c.x()), (int)  (vScale*c.y()), (int) (hScale*(c.x() + c.fx())), (int) (vScale*(c.y() + c.fy())));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 60f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 60f)), (int) (hScale * (c.x() + off + c.f2x() * cm / 12f)), (int) (vScale * (c.y() + off + c.f2y() * cm / 12f)));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 60f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 60f)), (int) (hScale * (c.x() + off + c.f3x() * cm / 24f)), (int) (vScale * (c.y() + off + c.f3y() * cm / 24f)));

                //g2d.fillOval((int) (hScale*(c.x()-1.2*c.minrad)), (int) (vScale*(c.y()-1.2*c.minrad)), (int) (2.4*hScale*c.minrad), (int) (2.4*vScale*c.minrad));

                if (c.weakDegrader) g2d.setColor(new Color(10, 160, 60, 100));
                else g2d.setColor(new Color(10, 60, 160, 100));

                g2d.setStroke(innerStroke);
                //g2d.drawLine((int) (hScale*c.x()), (int)  (vScale*c.y()), (int) (hScale*(c.x() + c.fx())), (int) (vScale*(c.y() + c.fy())));
                //g2d.drawLine((int) (hScale*(c.x()+off)), (int)  (vScale*(c.y()+off)), (int) (hScale*(c.x() + off + c.f2x()*cm/15f)), (int) (vScale*(c.y() + off + c.f2y()*cm/15f)));
                //g2d.drawLine((int) (hScale*(c.x()+off)), (int)  (vScale*(c.y()+off)), (int) (hScale*(c.x() + off + c.f3x()*cm/24f)), (int) (vScale*(c.y() + off + c.f3y()*cm/24f)));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 60f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 60f)), (int) (hScale * (c.x() + off + c.f2x() * cm / 12f)), (int) (vScale * (c.y() + off + c.f2y() * cm / 12f)));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 60f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 60f)), (int) (hScale * (c.x() + off + c.f3x() * cm / 24f)), (int) (vScale * (c.y() + off + c.f3y() * cm / 24f)));

            }else{

                Stroke outerStroke = new BasicStroke((cm/0.5f)*(float)(hScale*vScale),
                        BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

                Stroke innerStroke = new BasicStroke((cm/0.7f)*(float)(hScale*vScale),
                        BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

                g2d.setColor(new Color(200, 200, 200, 200));
                g2d.setStroke(outerStroke);


                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 30f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 30f)), (int) (hScale * (c.x() + off + c.f2x() * cm / 6f)), (int) (vScale * (c.y() + off + c.f2y() * cm / 6f)));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 30f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 60f)), (int) (hScale * (c.x() + off + c.f3x() * cm / 12f)), (int) (vScale * (c.y() + off + c.f3y() * cm / 12f)));

                g2d.setStroke(innerStroke);
                g2d.setColor(new Color(10, 160, 60, 100));

                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 30f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 30f)), (int) (hScale * (c.x() + off + c.f2x() * cm / 6f)), (int) (vScale * (c.y() + off + c.f2y() * cm / 6f)));
                g2d.drawLine((int) (hScale * (c.x() + off - c.f2x() * cm / 30f)), (int) (vScale * (c.y() + off - c.f2y() * cm / 30f)), (int) (hScale * (c.x() + off + c.f3x() * cm / 12f)), (int) (vScale * (c.y() + off + c.f3y() * cm / 12f)));
            }
        }

        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(1.0F));
        if(selectionBounds!=null){

            Point p1 = convertPointToDisplaySpace(new Point(selectionBounds.x,selectionBounds.y));
            Point p2 = convertPointToDisplaySpace(new Point(selectionBounds.width,selectionBounds.height));


            if(selectionIsOval) g2d.drawOval(p1.x,p1.y,p2.x,p2.y);
            else g2d.drawRect(p1.x, p1.y, p2.x, p2.y);
        }

        g2d.dispose();
    }

    public void setEnvironment(double c){

        boolean ps = ms.paused;
        if(!ps) ms.paused = true;

        for(int i=0; i<environment.profile.size(); i++){
            for(int j=0; j<environment.profile.get(i).size(); j++){
                System.out.println(selectionBounds.x);
                if(i*environment.grain >= selectionBounds.x && i*environment.grain <= selectionBounds.x+selectionBounds.width
                        && j*environment.grain >= selectionBounds.y && j*environment.grain <= selectionBounds.y+selectionBounds.height){


                    if(!selectionIsOval) environment.profile.get(i).get(j).c = c;
                    else{
                        Ellipse2D e2d = new Ellipse2D.Double(selectionBounds.x,selectionBounds.y,selectionBounds.width,selectionBounds.height);
                        if(e2d.contains(i*environment.grain,j*environment.grain)){

                            environment.profile.get(i).get(j).c_m2 = c;
                            environment.profile.get(i).get(j).c_m1 = c;
                            environment.profile.get(i).get(j).c = c;


                        }
                    }
                    // environment.profile[i][j] = c;
                }
            }
        }
        if(!ps) ms.paused = false;
    }

    Point convertPointToSimulationSpace(Point p){
        int width = this.frame.getContentPane().getWidth();
        int height = this.frame.getContentPane().getHeight();

        double hScale = (width-2.0*border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-2.0*border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Point p2 = p;

        p2.x = (int) ((p.x-border)/hScale);
        p2.y = (int) ((p.y-border)/vScale);

        return p2;
    }

    Point convertPointToDisplaySpace(Point p){
        int width = this.frame.getContentPane().getWidth();
        int height = this.frame.getContentPane().getHeight();

        double hScale = (width-2.0*border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-2.0*border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Point p2 = p;

        p2.x = (int) ((p.x*hScale)+border);
        p2.y = (int) ((p.y*vScale)+border);

        return p2;

    }

    public static Color ColorBlend(Color[] colours, double[] positions, double value){
        if(colours.length!=positions.length||colours.length<2) return new Color(0,0,0,0);
        if(value<positions[0]) value = positions[0];
        if(value>positions[positions.length-1]) value = positions[positions.length-1];

        int i = 0;
        for(i=0; i<colours.length-1; i++){

            if(value>=positions[i]&&value<=positions[i+1]) break;
        }

        double distance = positions[i+1]-positions[i];
        double d2       = (value - positions[i])/distance;
        double d1       = (positions[i+1]-value)/distance;

        double r = 0.5*((colours[i].getRed()*d1)+(colours[i+1].getRed()*d2));
        double g = 0.5*((colours[i].getGreen()*d1)+(colours[i+1].getGreen()*d2));
        double b = 0.5*((colours[i].getBlue()*d1)+(colours[i+1].getBlue()*d2));
        double a = 0.5*((colours[i].getAlpha()*d1)+(colours[i+1].getAlpha()*d2));

        return new Color((int) r, (int) g, (int) b, (int) a);
    }
}
