/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author John R. Brandon <jbrandon@gmail.com> or <jbrandon@lgl.com>
 */

import java.awt.*;

//import Jama.*;
//import Jama.LUDecomposition;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.*;

import java.util.Hashtable;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.*;

//import java.io.FileReader;
//import java.io.FileWriter;
//import java.io.IOException;

import java.awt.Color;

import javax.swing.JPanel;
import javax.swing.BorderFactory;
import javax.swing.plaf.ColorUIResource;
import javax.swing.border.Border;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.data.xy.YIntervalSeries;
import org.jfree.data.xy.YIntervalSeriesCollection;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.chart.annotations.XYTextAnnotation;

import org.jfree.ui.RefineryUtilities;
import org.jfree.ui.TextAnchor;

import org.jfree.chart.plot.DatasetRenderingOrder;

public class AppletGraph extends JApplet {
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    
    private javax.swing.JPanel panel_1;     // panel 1 holds the chart
    private javax.swing.JPanel panel_2;     // panel 2 holds the slider bars
      
    private YIntervalSeriesCollection N_data = new YIntervalSeriesCollection();
    private YIntervalSeries s1 = new YIntervalSeries("Abundance estimates");

    private DefaultXYDataset defaultXYDataset = new DefaultXYDataset();
    private DefaultXYDataset calfXYDataset = new DefaultXYDataset();
    private XYSeriesCollection catchComm = new XYSeriesCollection();
    private XYSeries catchSeries_tmp = new XYSeries("Catch");

    private double prop_mat_a[];	// 'equilibrium' proportion mature at age, based on assuming knife-edged maturity at age 7

    private double Nage_tot_0[], Nage_tot_1[];        // total numbers-at-age
    private double Nage_imm_0[], Nage_imm_1[];	// numbers-at-age that are immature
    private double Nage_mat_0[], Nage_mat_1[];	// numbers-at-age that are mature
    private double Nage_recpt_0[], Nage_recpt_1[];      // number-at-age receptive
    private double Nage_calvn_0[], Nage_calvn_1[];      // number-at-age with calf

    private double sum_NPR;             // sum over ages to get total NPR
    private double NPR_age[];           // Numbers-at-age-per-recruit

    private double initial_oneplus;     // size of female 1+ component on per-recruit basis given Finit
    private double temp_1plus;
    private double temp_mat;
    private double NPR_oneplus;         // size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing
    private double sum_Nage;
    private double b_eq;                // equilibrium birth rate at carrying capacity
    private double b_1;                 // Initial birth rate - corresponding with N_init/K depletion, for burning-in stable age/stage structure
    private double Finit;               // initial F
    private double f1_newt;             // Parameters used in Newton's algorithm to find value of Finit which results in same population size in 1900 as 1930
    private double objf_newt;
    private double df_newt;
    private double dx_newt;
    private double rtnewt;
    private double delta_Finit;
    private double N_init;              // Population size in first year of projections (e.g., 1846 or 1856)

    private double depl_init;           // Depletion in first year of trajectories
    private double depl[][];            // Depletion level each year

    private double select_init[];       // Initial selectivity associated w/historical aboriginal hunt (assumed uniform on 5+)
    private double select_histComm[];   // Selectivity associated w/historical commercial hunt (estimated??)
    private double select_modern[];     // Selectivity associated w/modern Russian aboriginal hunt (assumed uniform on 5+)

    private double b_t;                 // birth rate at year t
    private double b_t_yr[][];          // keep track of birth rate through time
    private double Nage[];		// numbers-at-age
    private double Male_age[];          // Vector with numbers-at-age for males
    private double Female_age[];        // "" for females
    private double Sage[];              // Vector with survival-at-age
    private double Page[];              // vector with rate of maturation at age - constant through time

    private double S_star;              // survival rate multiplier during mortality event of '99-'00
    private double s_pre;               // 1+ survival rate pre-1900
    private double s_post;              // 1+ survival rate post-1900 (if dictated by radio button selection)

    public int num1, num2, num3, num4, num5, num6;
    public int a_m;                                       // age at maturity
    public double delt_s, s_a, f_max, f_eq, z;
    public double b_max;                                  // max. birth rate

    private double Ntot[][];              // Total (0+) population size each year, indexed to last year +1 - final pop size is at start of last year (2005?) + 1
    private double Nplus[][];             // Vector of 1+ population size for all years
    private double N_calf[][];            // Vector of calf production for all years
    private double N_dead[][];            // Vector of natural dead each year
    private double NAll[][][];              // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males

    private double N_m_yr_1[][];            // Numbers of males by year
    private double N_calvn_yr_1[][];        // Number of calving females by year
    private double N_recpt_yr_1[][];        // Number of receptive females by year
    private double N_imm_yr_1[][];          // Number of immature females by year
    private double N_matfem_yr[][];         // Number of mature females by year

    public int rr_scaler, catchComm_scaler;
    public double catchMult;
    public double N_low, N_high;

    private Double xLowerBound = null;
    private Double xUpperBound = null;

    private double rr, KK;          // r and K for surplus production model(s)
    private double N_t[][];         // Numbers at year t (using for logistic model)
    private double K_t[][];         // Time series of Carrying Capacity
    private double K_mult;          // Carrying Capacity multiplier

    private double N_obs[][];           // Abundance estimates
    private double CV_N[][];            // CVs of the abundance estimates (ignoring Covariance)
    private double CV_Add;              // Additional variation about abundance estimates, not taken into account by survey methods
    private double VarCov[][];          // Variance covariance matrix of the abundance estimates
    public int N_n, dim1;               // number of years with data [N_n] and number of yrs of projection [dim1]
    public int first_yr, last_yr;       // number of years: eventually probably want to read this from a file

    private double Catch[][];                    // Time series of total catches
    private double catchMaleComm[][];           // Time series of historical commercial male catches
    private double catchFemaleComm[][];         // Time series of "" female catches
    private double catchMaleAborig[][];         // Time series of aboriginal male catches
    private double catchFemaleAborig[][];       // Time series of "" female catches
    private double catch_CaliforniaShore[][];   // Time series of commercial catches from CA shore stations
    private double catch_Mex[][];               // Time series of commercial catches from Mexico waters (incl. lagoons)
    private double catchMale_Mex[][];           // Males in Mexico catches
    private double catchFemale_Mex[][];         // Females in Mexico catches
    private double catchMale_Other[][];         // Time series of all male catches excluding Mexico (in/out lagoons) and CA shore stations
    private double catchFemale_Other[][];         // Time series of all male catches excluding Mexico (in/out lagoons) and CA shore stations

    private double phi_fem;                    // Sex ratio of historic commercial catch (i.e. 1.0 = 100%females)
    private double beta_cow;                    // Proportion (mature) female catch: cows w/calves

    private double eps_Mex;                    // Mex catch magnitude multipler for landings (takes into account SE calculated by Reeves et al 2010)
    private double LRF_Mex;                    // Loss Rate Factor (LRF) for Mex Catches
    private double eps_CA;                     // CA shore based magnitude multipler for landings (takes into account SE calculated by Reeves and Smith 2010)
    private double LRF_CA;                     // Loss Rate Factor (LRF) for CA shore based catches
    private double C_mult;                      // Generic catch multiplier appled to all fisheries prior to 1900

    private double C_out[][];		// Ship based catches in Mexico which were outside lagoons; dim = (1846,1873)
    private double C_in[][];		// All catches in lagoons during lagoon whaling period; dim = (1846,1873)
    private double C_mex_tmp[][];		// "" taking into account uncertainty in magnitude of catch numbers (i.e. after applying catch multiplier eps_mex);; dim = (1846,1873)

    private double C_fem_mod_out[][];	// catches by sex, after magnitude uncertainty multipliers, all areas other than lagoons during that period; dim = (1846,1873)
    private double C_male_mod_out[][];	// dim = (1846,1873)
    private double C_fem_mod_in[][];    //
    private double C_male_mod_in[][];	// catches by sex after applying estimates of catch coefficients (sex ratio etc.) in lagoons
    private double C_male[][];          // total number of male catches in a given year (all areas)
    private double C_fem[][];          // total number of female catches in a given year (all areas)

    private double C_c_mod_in[][];		// catches of cows w/calves, after beta (proportion cows w/calves) uncertainty multipliers
    private double C_noncow_mod_in[][];         // catches of non cow w/calves (immature and receptive) females, ""

    private double E_c_in[][];		// exploitation rate of cows w/calves in lagoons
    private double E_noncow_in[][];	// exploitation rate of non cow w/calves (assumed equal for vulnerable immature and receptive) females in lagoons
    private double E_out[][];		// exploitation rate outside lagoons (assumed equal for all female reproductive stages)

    private double E_i_yr[][];          // annual exploitation rate on immature animals
    private double E_r_yr[][];          // " " receptive females
    private double E_c_yr[][];          // " " calving females

    private double C_init;              // calculate aboriginal catch in 1845 implied by Finit and population size in 1845

    private double matFem_geneticN;
    private double matFem_genetic[][];          // genetic estimate of mature females at bottleneck

    private double Rec_f_yr[];                  //
    private double Rec_m_yr[];                  //
    private double Rec_f_imm_yr[];               // recruited immature females
    private double Rec_f_recp_yr[];              // recruited mature receptive females
    private double Rec_f_calv_yr[];              // recruited mature calving females

    private double F_f_imm;
    private double F_f_recp;
    private double F_f_calvn;

    private double depl_star;                   // depletion level below which depensation kicks-in

    private double NegLogLike_Abun;             // Negative Log-likelihood and variables used in its calculations (below)
    private double negLL_tmp;
    private double det_varcov;
    private double tmp_VarCov[][];            //(1,N_n,1,N_n);
    private double inv_VarCov[][];
    private double cvadd_matrix[][];
    private double cvadd_matrix2[][];
    private double log_resids[];
    private double tmp_product[];

    private boolean initialized;               // flag to see if applet has been initialized yet

    boolean Surv_const;     // If survival is constant through time = T; else (pre- and post-1900) = F

//        JSlider scroll_1 = new JSlider(JSlider.HORIZONTAL, 75, 100, (int)(depl_init*100));  // Depletion initial year (currently 1846)
    private JSlider scroll_2;

//        JSlider scroll_3 = new JSlider(JSlider.HORIZONTAL, 995, 1010, 1000);                // K multiplier
//        JSlider scroll_4 = new JSlider(JSlider.HORIZONTAL, 0, 40, (int)(depl_star*100));    // Depensation (depletion level below which depensation occurs)
//        JSlider scroll_5 = new JSlider(JSlider.HORIZONTAL, 0, 100, (int)(phi_fem*100));     // Sex ratio of historic commercial catch (i.e. 1.0 = 100%females)
//        JSlider scroll_6 = new JSlider(JSlider.HORIZONTAL, 0, 100, (int)(beta_cow*100));    // Proportion female commercial catch: cows w/calves (i.e. 1.0 = 100% of females = cows w/calves)
//        final JSlider scroll_7 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_a*1000));      // 1+ survival (constant or pre-1900 as dictated by radio button)
//        final JSlider scroll_77 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_post*1000));      // 1+ survival post-1900
//        final JSlider scroll_777 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_a*1000));      // 1+ survival post-1900
//        JSlider scroll_8 = new JSlider(JSlider.HORIZONTAL, 1, 500, (int)(delt_s*1000));     // Delta survival, difference between adult and calf survival
//        JSlider scroll_9 = new JSlider(JSlider.HORIZONTAL, 6, 12, a_m);                     // Age at maturity
//        JSlider scroll_10 = new JSlider(JSlider.HORIZONTAL, 100, 1000, (int)(z*100));       // Density dependent shape parameter (Z)
//
//        JSlider scroll_11 = new JSlider(JSlider.HORIZONTAL, -35, 35, (int)(eps_Mex*100));   // Mex catch magnitude unertainty multipler for landings
//        JSlider scroll_12 = new JSlider(JSlider.HORIZONTAL, 110, 152, (int)(LRF_Mex*100));  // LRF for Mex Catches
//
//        JSlider scroll_13 = new JSlider(JSlider.HORIZONTAL, 120, 210, (int)(LRF_CA*100));   // LRF for Mex Catches
//        JSlider scroll_14 = new JSlider(JSlider.HORIZONTAL, -52, 52, (int)(eps_CA*100));    // CA catch magnitude unertainty multipler for landings
//
//        JSlider scroll_15 = new JSlider(JSlider.HORIZONTAL, 0, 35, (int)(CV_Add*100));      // CV_add for abundance estimates
//
//        JSlider scroll_16 = new JSlider(JSlider.HORIZONTAL, 20, 100, (int)(S_star*100));      // Survival rate multiplier during mortality event of '99 and '00
//
//        JSlider scroll_17 = new JSlider(JSlider.HORIZONTAL, 50, 200, (int)(C_mult*100));      // Catch (removals) multiplier ; generic and applying it to all fisheries prior to 1900
//        final JSlider scroll_18 = new JSlider(JSlider.HORIZONTAL, 0, 2000, (int)(C_init));      // Catch (removals) multiplier ; generic and applying it to all fisheries prior to 1900

//    public JRadioButton rb1;
//    public JRadioButton rb2;

/**
 *
 * @Override
 */
    @Override

//#######################
    public void init() {
        int Year;
        Border empty;
// START INITIALIZATIONS ////////////////////////////////////////////////////////

        setSize(1200,800);              // determines the dimensions of the applet
//        setSize(1024,768);              // determines the dimensions of the applet
//        setSize(1200,900);              // determines the dimensions of the applet

//        KK = 21500;                     // initial value for carrying capacity
//        K_mult = 1.0;                   // Carrying Capacity multiplier
//
//        phi_fem = 0.80;                 // Sex ratio of historic commercial catch (ie. 1.0 = 100% females)
//        beta_cow = 0.75;                // Proportion (mature) females in commercial catch that were cows w/calves
//
//        depl_init = 0.95;                // depletion (in terms of 1+ ?) during first year of trajectories
//
//        matFem_geneticN = 1000;         //
//
//        a_m = 10;                        // Age-at-maturity (equal to the age-at-first partuition less one)
//        delt_s = 0.360;                  // Difference between calf and 1+ survival
//        s_a = 0.960;                    // Survival of animals age 1+
//        s_post = 0.960;                 // Survival of animals age 1+ after 1900
//        Surv_const=true;                // Initially survival is assumed constant through time
//        z = 2.39;                       // Density dependent shape parameter
//        b_max = 0.99;                   // maximum birth rate, i.e. 100% of mature receptive could calf next year
//        S_star = 0.86;                    // Survival rate multiplier during mortality event of '99 and '00
//
//        N_low = 0;                      //
//        N_high = 40000;                 // y-limit for trajectory plot
//
//        N_n = 23;                       // number of abundance estimates
//        N_obs = new double[2][N_n];     // N_obs[0][] = YEar of estimate; N_obs[1][] = Estimate;
//        CV_N = new double[2][N_n];      // CV_obs[0][] = Year of estimate; CV_obs[1][] = CV of Estimate;
//        CV_Add = 0.00;                   // additional variation about abundance estimates, not taken into account by survey methodology
//
//        first_yr = 1846;                // first year of projection
//        last_yr = 2010;                 // last year of projection
//        dim1 = last_yr - first_yr + 1;  // for indexing
//
//        Catch = new double[2][dim1];    // Year and combined commercial catch (both sexes)
//        catchMaleComm = new double[2][dim1];    // male commercial catches
//        catchFemaleComm = new double[2][dim1];  // female commercial catches
//        eps_Mex = 0.00;                     // Mex catch magnitude multipler for landings (takes into account SE calculated by Reeves et al 2010)
//        LRF_Mex = 1.24;                    // Loss Rate Factor (LRF) for Mex Catches
//        eps_CA = 0.0;                      // CA shore based magnitude multipler for landings (takes into account SE calculated by Reeves and Smith 2010)
//        LRF_CA = 1.20;                     // Loss Rate Factor (LRF) for CA shore based catches proposed by Reeves and Smith (2010)
//        C_mult = 1.0;                       // Generic catch multiplier applied to all fisheries prior to 1900
//
//        depl_star = 0.00;       // depletion level at which depensation kicks-in
//
//        rr_scaler = 1000;       // scale population growth rate (slider goes from 0 to 1000)
//        catchComm_scaler = 100;       // scale population growth rate (slider goes from 0 to 1000)
//        catchMult = 1.0;
//
//        Nplus = new double[2][dim1+1];                // Vector of total age 1+ population
//        N_dead = new double[2][dim1+1];               // number of dead animals per year
//        Ntot = new double[2][dim1+1];
//        b_t_yr = new double[2][dim1+1];               // birth rate through time
//        N_imm_yr_1 = new double[2][dim1+1];
//        N_m_yr_1 = new double[2][dim1+1];
//        N_calf = new double[2][dim1+1];
//        N_calvn_yr_1 = new double[2][dim1+1];
//        N_recpt_yr_1 = new double[2][dim1+1];
//        N_matfem_yr = new double[2][dim1+1];
//
//        Rec_f_yr = new double[dim1+1];
//        Rec_m_yr = new double[dim1+1];
//
//        NAll = new double[5][dim1+1][a_m+1];               // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males
//        depl = new double[2][dim1+1];                     // depletion each year
//        N_t = new double[2][dim1+1];
//        K_t = new double[2][dim1+1];
//        matFem_genetic = new double[2][dim1+1];                // genetic estimate of mature females at bottleneck

//        catchSeries_tmp = readCatchData();                               // catch data is hard-coded
//        catchComm.addSeries(catchSeries_tmp);
// END INITIALIZATIONS ////////////////////////////////////////////////////////
        initialized = false;
        Default_Values();

        for (Year=0 ;Year < dim1 + 1 ;Year++)           // Loop over years for projections
        {
           Nplus[0][Year] = Year + first_yr;
           Ntot[0][Year] = Nplus[0][Year];
           N_calf[0][Year] = Nplus[0][Year];
           K_t[0][Year] = Nplus[0][Year];
           N_matfem_yr[0][Year] = Nplus[0][Year];
           matFem_genetic[0][Year] = Nplus[0][Year];
           matFem_genetic[1][Year] = matFem_geneticN;
        }
//
//// END INITIALIZATIONS ////////////////////////////////////////////////////////
        catchSeries_tmp = readCatchData();                               // catch data is hard-coded
        catchComm.addSeries(catchSeries_tmp);
//
//// CONSIDER MOVING THIS CODE; READ IN THE CATCH DATA DURING INITIALIZATION SECTION OF PROGRAM
        readData();             // read in abundance estimates and associated CVs (and yrs)

        updateModel();           //   update logistic model if sliders are changed

        JFreeChart jFreeChart = generateGraph();                    // Get the graph (generateGraph will create the JFreeChart graph and add the red and blue point on it).

        ChartPanel chartPanel = new ChartPanel(jFreeChart);         //Put the jFreeChart in a chartPanel
        chartPanel.setPreferredSize(new Dimension(300,400));        // dimension for the chart

        chartPanel.setPopupMenu(null);
        Container content = this.getContentPane();          //add the chartPanel to the container (getContentPane is inherited from JApplet which AppletGraph extends).

        jFreeChart.setTitle("Eastern Gray Whale Population Dynamics Model");  // N[t+1] = N[t] + N[t] * r * { 1 - ( N[t] / K ) }

 //       JTabbedPane tabbedPane1 = new JTabbedPane();       // going to try and make chart panel have tabs

//        JPanel panel_0 = new JPanel();
        JPanel panel_1 = new JPanel();

//        tabbedPane1.add("Population Trajectories and Catch History", chartPanel);
//        tabbedPane1.add("Depensation", panel_1);
//        content.add(tabbedPane1);

        empty = BorderFactory.createEmptyBorder();
        chartPanel.setBorder(empty);

        content.add(panel_1);


        content.add(chartPanel);                          // adding the chart using this method makes it fill the panel

////////////////////////////////////////////////////////////////////////////////

        addSliders(content);                     // calls a function which adds the sliders below the JFreeChart object

//        testJama();

        //        System.out.println(Surv_const);

    }

//#######################
   public void Default_Values(){
        int Year;

        KK = 21500;                     // initial value for carrying capacity
        K_mult = 1.0;                   // Carrying Capacity multiplier

        phi_fem = 0.80;                 // Sex ratio of historic commercial catch (ie. 1.0 = 100% females)
        beta_cow = 0.75;                // Proportion (mature) females in commercial catch that were cows w/calves

        depl_init = 0.95;                // depletion (in terms of 1+ ?) during first year of trajectories

        matFem_geneticN = 1000;         //

        a_m = 10;                        // Age-at-maturity (equal to the age-at-first partuition less one)
        delt_s = 0.360;                  // Difference between calf and 1+ survival
        s_a = 0.960;                    // Survival of animals age 1+
        s_post = 0.960;                 // Survival of animals age 1+ after 1900
        Surv_const=true;                // Initially survival is assumed constant through time
        z = 2.39;                       // Density dependent shape parameter
        b_max = 0.99;                   // maximum birth rate, i.e. 100% of mature receptive could calf next year
        S_star = 0.86;                    // Survival rate multiplier during mortality event of '99 and '00

        N_low = 0;                      //
        N_high = 40000;                 // y-limit for trajectory plot

        N_n = 23;                       // number of abundance estimates
        N_obs = new double[2][N_n];     // N_obs[0][] = YEar of estimate; N_obs[1][] = Estimate;
        CV_N = new double[2][N_n];      // CV_obs[0][] = Year of estimate; CV_obs[1][] = CV of Estimate;
        CV_Add = 0.00;                   // additional variation about abundance estimates, not taken into account by survey methodology

        first_yr = 1846;                // first year of projection
        last_yr = 2010;                 // last year of projection
        dim1 = last_yr - first_yr + 1;  // for indexing

        Catch = new double[2][dim1];    // Year and combined commercial catch (both sexes)
        catchMaleComm = new double[2][dim1];    // male commercial catches
        catchFemaleComm = new double[2][dim1];  // female commercial catches
        eps_Mex = 0.00;                     // Mex catch magnitude multipler for landings (takes into account SE calculated by Reeves et al 2010)
        LRF_Mex = 1.24;                    // Loss Rate Factor (LRF) for Mex Catches
        eps_CA = 0.0;                      // CA shore based magnitude multipler for landings (takes into account SE calculated by Reeves and Smith 2010)
        LRF_CA = 1.20;                     // Loss Rate Factor (LRF) for CA shore based catches proposed by Reeves and Smith (2010)
        C_mult = 1.0;                       // Generic catch multiplier applied to all fisheries prior to 1900

        depl_star = 0.00;       // depletion level at which depensation kicks-in

        rr_scaler = 1000;       // scale population growth rate (slider goes from 0 to 1000)
        catchComm_scaler = 100;       // scale population growth rate (slider goes from 0 to 1000)
        catchMult = 1.0;

        Nplus = new double[2][dim1+1];                // Vector of total age 1+ population
        N_dead = new double[2][dim1+1];               // number of dead animals per year
        Ntot = new double[2][dim1+1];
        b_t_yr = new double[2][dim1+1];               // birth rate through time
        N_imm_yr_1 = new double[2][dim1+1];
        N_m_yr_1 = new double[2][dim1+1];
        N_calf = new double[2][dim1+1];
        N_calvn_yr_1 = new double[2][dim1+1];
        N_recpt_yr_1 = new double[2][dim1+1];
        N_matfem_yr = new double[2][dim1+1];

        Rec_f_yr = new double[dim1+1];
        Rec_m_yr = new double[dim1+1];

        NAll = new double[5][dim1+1][a_m+1];               // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males
        depl = new double[2][dim1+1];                     // depletion each year
        N_t = new double[2][dim1+1];
        K_t = new double[2][dim1+1];
        matFem_genetic = new double[2][dim1+1];                // genetic estimate of mature females at bottleneck

        for (Year=0 ;Year < dim1 + 1 ;Year++)           // Loop over years for projections
        {
           Nplus[0][Year] = Year + first_yr;
           Ntot[0][Year] = Nplus[0][Year];
           N_calf[0][Year] = Nplus[0][Year];
           K_t[0][Year] = Nplus[0][Year];
           N_matfem_yr[0][Year] = Nplus[0][Year];
           matFem_genetic[0][Year] = Nplus[0][Year];
           matFem_genetic[1][Year] = matFem_geneticN;
        }

        if(initialized=true){   // if slider bars exist; reset their values to defaults
//            scroll_2.setValue((int)KK);
//            addSliders(this);

        }

// END INITIALIZATIONS ////////////////////////////////////////////////////////
//        catchSeries_tmp = readCatchData();                               // catch data is hard-coded
//        catchComm.addSeries(catchSeries_tmp);

// CONSIDER MOVING THIS CODE; READ IN THE CATCH DATA DURING INITIALIZATION SECTION OF PROGRAM
//        readData();             // read in abundance estimates and associated CVs (and yrs)
//        updateModel();           //   update logistic model if sliders are changed
//        JFreeChart jFreeChart = generateGraph();                    // Get the graph (generateGraph will create the JFreeChart graph and add the red and blue point on it).

   initialized = true;      // set flag
   }
//#######################
/**
 *
 * @Override
 */
    @Override
    public void start(){
        super.start();
        this.repaint();
    }

//#######################
//public void testJama(){
//      double[][] array = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
//      Matrix A = new Matrix(array);
//      Matrix b = Matrix.random(3,1);
//      LUDecomposition LU = new LUDecomposition(b);
//      double dd = A.det();
//      Matrix x = A.solve(b);
//      Matrix Y = A.plus(A); //Matrix(array);
//    }
//#######################
    public void NegLL_Abundance(){
        int ii, jj;      // counters

        NegLogLike_Abun = 0;
        det_varcov = 0;

        tmp_product = new double[N_n];
//        cvadd_matrix = new double[N_n][N_n];
        cvadd_matrix2 = new double[N_n][N_n];
        tmp_product = new double[N_n];
        VarCov = new double[N_n][N_n];

        for (jj=1;jj<=N_n;jj++) cvadd_matrix2[jj][jj]=CV_Add*CV_Add;

// for testing : population a variance covariance matrix with just the variances on the diagonals (no covariances, yet)
        for (jj=1;jj<=N_n;jj++) VarCov[jj][jj]=CV_N[1][jj]*CV_N[1][jj];

//        tmp_VarCov = VarCov.plus(cvadd_matrix2);// + cvadd_matrix2;
// det_varcov = det(tmp_VarCov);
// det_varcov = log(det_varcov);
// det_varcov = 0.50*det_varcov;
//
//// NegLogLike_Abun=0.5*log(det_varcov);
//
// for (ii=1;ii<=N_n;ii++)
//  if(N_hat(ii)>1)
//   {
//     log_resids(ii)=log(Nplus(N_yrs(ii)))-log(N_hat(ii));
//   }
//
// inv_VarCov = inv(tmp_VarCov);
// tmp_product = log_resids * inv_VarCov;
// NegLogLike_Abun = tmp_product * log_resids;
//
// NegLogLike_Abun = det_varcov + 0.5*NegLogLike_Abun;

    }

//#######################
    public void Assign_KK(){
        int Year;

        K_t[1][0] = KK;                                  // initial carrying capacity
        for (Year = 1 ;Year < dim1 + 1; Year++){         // Loop over years for projections
            K_t[1][Year] = KK + KK*(K_mult-1)*Year;               // Linear increase through time
        }
    }

//#######################
    public void Do_projections(){

    int Year, age, ii, jj;        		// Counters
    double Shunt_f, Shunt_m;
    double tot_catch, tot_rec, Temp;
    double Fmal, Ffem;
    double Rec_f_tot, Rec_m_tot;
    double Rec_i_tot, Rec_r_tot, Rec_c_tot;
    double C_fem_tmp_out, C_male_tmp_out;       // temp mex. catches outside lagoons
    double C_fem_tmp_in, C_male_tmp_in;         // temp mex. catches inside lagoons
    double depl_tmp; 
    int AR;                                     // age-at-recruitment
    int ex_ii;                                  // index for dealing with extinctions
    double debug, b_tmp;
    double TempA[][];              // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males
    boolean surv_flag;

    TempA = new double[5][a_m+1];               // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males
    select_modern = new double[a_m+1];          // basic selectivity-at-age for modern russian whaling
    select_histComm = new double[a_m+1];        // basic selectivity-at-age for historical comm. whaling
    Rec_f_imm_yr = new double[dim1+1];          // recruited immature females
    Rec_f_recp_yr = new double[dim1+1];         // recruited mature receptive females
    Rec_f_calv_yr = new double[dim1+1];         // recruited mature calving females
    C_c_mod_in = new double[2][dim1];
    C_noncow_mod_in = new double[2][dim1];
    C_male = new double[2][dim1];               // male catches all areas
    C_fem = new double[2][dim1];                // female catches all areas

    E_c_in = new double[2][dim1];		// exploitation rate of cows w/calves in lagoons
    E_noncow_in = new double[2][dim1];	// exploitation rate of non cow w/calves (assumed equal for vulnerable immature and receptive) females in lagoons
    E_out = new double[2][dim1];		// exploitation rate outside lagoons (assumed equal for all female reproductive stages)

    E_i_yr = new double[2][dim1];
    E_r_yr = new double[2][dim1];
    E_c_yr = new double[2][dim1];

    F_f_imm = 0.0;
    F_f_recp = 0.0;
    F_f_calvn = 0.0;

    surv_flag=true;
    debug=0;

    AR = 5;                         // assumed age at recruitment into the modern russian hunt
    for(ii=0;ii<=a_m;ii++){         // this code could move to init()
       if(ii < AR){
           select_modern[ii]=0.0;
           select_histComm[ii]=0.0;
       }
       else{
           select_modern[ii]=1.0;
           select_histComm[ii]=1.0;
       }
    }

// NEED TO ADD TWO LOOPS (LIKE READ_CATCH CODE; DIFFERENT SELECTIVITIES DEPENDING ON YEAR, ETC.
    for(Year = 0; Year <= dim1-1; Year++) {          // perform modifications to catches up to 1899 (end of CA shore whaling)
////////////////////////////////////////////////////////////////////////////////////////////////////
// LAGOON (1846-1873) & CA SHORE WHALING: 1846 - 1899
////////////////////////////////////////////////////////////////////////////////////////////////////
        C_fem_tmp_out=C_out[1][Year]+catchFemale_Other[1][Year];	// catches by sex all areas other than lagoons during that period
        C_male_tmp_out=C_out[1][Year]+catchMale_Other[1][Year];      // male catches outside lagoons (note C_out is half of Mex. catches outside lagoons, hence assumption of 1:1 ratio taken into account, calculations done when reading catch data)
        C_fem_tmp_in=C_in[1][Year]*phi_fem;                          // female catches inside lagoons
        C_male_tmp_in=C_in[1][Year]*(1-phi_fem);                     // male catches inside lagoons

// apply beta (proportion cows w/calves in female catches in lagoons) uncertainty multiplier
        C_c_mod_in[1][Year]=C_fem_tmp_in*beta_cow;            // catches of cows w/calves, after beta (proportion cows w/calves) uncertainty multipliers
        C_noncow_mod_in[1][Year]=C_fem_tmp_in*(1-beta_cow);   // catches of non-cow w/calves (immature and receptive) females, ""

// sum the number of recruited males and females by reproductive stage (allowing for alternative selectivities on females by reproductive stage during the lagoon whaling period)
        Rec_i_tot = 0; Rec_r_tot = 0; Rec_c_tot = 0; Rec_m_tot = 0;

        for(age=0; age<=a_m; age++)
        {
         Rec_i_tot = Rec_i_tot + select_histComm[age]*NAll[1][Year][age]; // immature females
         Rec_r_tot = Rec_r_tot + select_histComm[age]*NAll[2][Year][age]; // mature receptive females
         Rec_c_tot = Rec_c_tot + select_histComm[age]*NAll[3][Year][age]; // cows w/calves
         Rec_m_tot = Rec_m_tot + select_histComm[age]*NAll[4][Year][age]; // males
        }
        Rec_f_yr[Year]=Rec_i_tot + Rec_r_tot + Rec_c_tot;
        Rec_m_yr[Year]=Rec_m_tot;

// calculate exploitation rates for females (each reproductive stage) and males during lagoon whaling
        E_c_in[1][Year] = C_c_mod_in[1][Year]/Rec_c_tot;	   		        // exploitation rate of cows w/calves in lagoons
        E_noncow_in[1][Year] = C_noncow_mod_in[1][Year]/(Rec_i_tot + Rec_r_tot);   // exploitation rate of non cow w/calves (assumed equal for vulnerable immature and receptive) females in lagoons
        E_out[1][Year] = C_fem_tmp_out/Rec_f_yr[Year];			        // exploitation rate of females (assumed equal for all reproductive stages) outside of lagoons

        C_male[1][Year] = C_male_tmp_out + C_male_tmp_in;
        Fmal = C_male[1][Year] / Rec_m_yr[Year];
        C_fem[1][Year] = C_fem_tmp_out + C_fem_tmp_in;       // not currently using this, but could be good to have for future reference / plotting catches by sex
        Ffem = C_fem[1][Year] / Rec_f_yr[Year];

// calculate total exploitation rates for both sexes for all whaling in a given year
        E_i_yr[1][Year] = E_noncow_in[1][Year] + E_out[1][Year];                       // total exploitation rate of sub-adult stage including catches in and out of lagoons
        E_r_yr[1][Year] = E_noncow_in[1][Year] + E_out[1][Year];                       // " " receptive stage " "
        E_c_yr[1][Year] = E_c_in[1][Year] + E_out[1][Year];                            // " " cows w/ calves " "
//        Fmal = (C_male_tmp_in(Year) + C_male_tmp_out(Year)) / Rec_m_yr(Year); // exploitation rate on males including ship mex and shore CA (after adjustments) and all other catches

       N_matfem_yr[1][Year] = 0;
        for (age=0; age<=a_m; age++)                         // Remove hunter mortality
            {
             N_matfem_yr[1][Year] += NAll[2][Year][age] + NAll[3][Year][age]; // sum mature females

//             TempA[1][age] = NAll[1][Year][age] * (1-select_histComm[age] * E_noncow_in[1][Year]) * (1-select_histComm[age] * E_out[1][Year]);    // (1) immature fem
//             TempA[2][age] = NAll[2][Year][age] * (1-select_histComm[age] * E_noncow_in[1][Year]) * (1-select_histComm[age] * E_out[1][Year]);    // (2) mature receptive fem
//             TempA[3][age] = NAll[3][Year][age] * (1-select_histComm[age] * E_c_in[1][Year]) * (1-select_histComm[age] * E_out[1][Year]);         // (3) mature calving fem, and
             TempA[1][age] = NAll[1][Year][age] * (1-select_histComm[age] * E_i_yr[1][Year]);    // (1) immature fem
             TempA[2][age] = NAll[2][Year][age] * (1-select_histComm[age] * E_r_yr[1][Year]);    // (2) mature receptive fem
             TempA[3][age] = NAll[3][Year][age] * (1-select_histComm[age] * E_c_yr[1][Year]);         // (3) mature calving fem, and

             TempA[4][age] = NAll[4][Year][age] * (1-select_histComm[age] * Fmal);                                                                 // (4) males
            }

// REMOVE CALF MORTALITY ASSOCIATED WITH KILLING A COW W/CALF if whaling is before 1900
        if (Year <= 53){
            TempA[1][0] = NAll[1][Year][0] * (1 - E_c_yr[1][Year]);
            TempA[4][0] = NAll[4][Year][0] * (1 - E_c_yr[1][Year]);
        }
//        else{
//            TempA[1][1] = NAll[1][Year][1];     // after 1900, majority of landings from northern waters, so most of young of the year at age of weaning; hence not assuming calf mortality associated with killing cow w/calf
//            TempA[4][1] = NAll[4][Year][1];
//        }

// model survival rates post-1900 (or at least mortality event)
        if(Year>=54){
// If 1999 or 2000, apply survival multiplier for mortality event
            if (Year == 153){                                       // 
                for (age=0; age<=a_m; age++){                        // Remove natural mortality and count deaths
                    Sage[age]=Sage[age]*S_star;
                }
            }else if (Year == 154){
// do nothing
            }else if (Surv_const==true) {           // Survival is constant after 1900, so just bring it back to normal after mortality event
                Sage[0] = s_a - delt_s;             // need this else statement to re-assign 'normal' rates post-mortality event
                for(ii = 1; ii <= a_m; ii++){        // assign calf survival
                    Sage[ii] = s_a;                    // fill survival-at-age vector
                }
            }else if (Surv_const==false){           // survival rate is different after 1900, so assign new rate here
                Sage[0]=s_post-delt_s;
                for (age=1; age<=a_m; age++){                         // Remove hunter mortality
                    Sage[age]=s_post;
                }
                Do_NPR(s_post);   // DEBUGGGING
             }
        }      

// Remove natural mortality and count deaths
        for (age=0; age<=a_m; age++){                        
            for (jj=1; jj<=4; jj++)
                {
                  N_dead[1][Year] += TempA[jj][age] * (1.0 - Sage[age]);
                  TempA[jj][age] = Sage[age] * TempA[jj][age];
                }
        }

// tally 1+ numbers
        Nplus[1][Year+1] = 0.0;                             
        for (age = 1; age <= a_m; age++){
            Nplus[1][Year+1] += TempA[1][age] + TempA[2][age] + TempA[3][age] + TempA[4][age];
        }

//        if(Nplus[1][Year+1]>=K_t[1][Year+1]){               // debugging
//            depl[1][Year+1]=depl[1][Year+1];
//        }

// Start each year as if January 1st, with population down in the breeding lagoons of baja - so first day of projections would be New Year's 1930.
        depl[1][Year+1] = Nplus[1][Year+1] / K_t[1][Year+1];      // Allow for time-varying carrying capacity
        depl_tmp = depl[1][Year+1];

        b_t = b_eq+(b_max-b_eq)*(1-Math.pow(depl[1][Year],z));        // calculate effect of density dependence on fecundity
        b_tmp = b_t;

// DEBUGGING
        if(b_t>1.0) {
          debug = 1;
          b_t = b_max;
          if(Year==25){
              debug=1;
          }
        }
        if(b_t<0.0){
          if(depl[1][Year]<1.0){
            debug = 1;
            b_t = b_max;
            if(Year==25){
              debug=1;
            }
          }else{
              debug = 1;
              b_t = b_eq;
          }
        }

// DEPENSATION TO DENSITY DEPENDENCE TERM
       if (depl[1][Year] <= depl_star)
           b_t = b_t * (b_eq + (1-b_eq) * depl[1][Year] / depl_star);

// IN ADMB CODE THIS WILL BE THE BIRTH RATE PENALTY; IN JAVA JUST RESETTING BIRTH RATES TO BE IN BOUNDS TO PREVENT PHOENIX RISING
        if(b_t>1.0) {
          debug = 1;
          b_t = b_max;
          if(Year==25){
              debug=1;
          }
        }                                        // impose biological constraints on birth rates

        if(b_t<0.0){
          if(depl[1][Year]<1.0){
              b_t = b_max;
          }else{                                // b_t can go negative at really low depletion levels evidently
              b_t = b_eq;
          }
          
        }

       b_t_yr[1][Year] = b_t;

// Now deal with maturation and reproduction - keeping track of PLUS GROUP dynamics for receptive and calving stages
       for (age=1; age<=a_m-1; age++)
        {
         NAll[1][Year+1][age] = TempA[1][age-1] * (1-Page[age-1]);	// immature females
         NAll[4][Year+1][age] = TempA[4][age-1];			// not dealing with male maturation - assuming they are not a limiting factor in female reproductive success
        }

        NAll[2][Year+1][a_m] = (1.0-b_t)*(TempA[1][a_m-1]*Page[a_m-1]+TempA[2][a_m])+TempA[3][a_m]; // mature receptive
        NAll[3][Year+1][a_m] = b_t*(TempA[1][a_m-1]*Page[a_m-1]+TempA[2][a_m]);                     // mature calving
        NAll[4][Year+1][a_m] = TempA[4][a_m-1] + TempA[4][a_m];                                     // males
        
        debug = NAll[1][Year+1][1];
        debug = NAll[2][Year+1][a_m];
        debug = NAll[3][Year+1][a_m];
        debug = NAll[4][Year+1][a_m];
// Transfer population components from array to vectors
       for (age=0; age<=a_m; age++) {
           N_recpt_yr_1[1][Year+1] += NAll[2][Year+1][age];
           N_calvn_yr_1[1][Year+1] += NAll[3][Year+1][age];
           N_imm_yr_1[1][Year+1] += NAll[1][Year+1][age];
           N_m_yr_1[1][Year+1] += NAll[4][Year+1][age];
       }

// CATCH EXTINCTIONS AND EXIT
        if(N_recpt_yr_1[1][Year] <= 0.0 || N_calvn_yr_1[1][Year+1] <= 0.0 || N_imm_yr_1[1][Year+1] <= 0 || N_m_yr_1[1][Year+1] <= 0.0) {                                     // in case of EXTINCTIONS
            for(ex_ii = Year; ex_ii <= dim1-1; ex_ii++){
                Nplus[1][ex_ii]=0.0;
                N_matfem_yr[1][ex_ii]=0.0;
            }
            return;                                         // send pop trajectories zero back to graphing routines
        }

// Get calves
       NAll[1][Year+1][0] = 0.5 * N_calvn_yr_1[1][Year+1];	//female calves (i.e. immature age=0)
       NAll[4][Year+1][0] = 0.5 * N_calvn_yr_1[1][Year+1];      //male calves

       NAll[2][Year+1][1] = 0;
       NAll[3][Year+1][1] = 0;

       Ntot[1][Year+1] = Nplus[1][Year+1] + N_calvn_yr_1[1][Year+1];	//total population = N_one plus + calves (N_0)

 }  // end looping over years
 N_calf = N_calvn_yr_1;     // vector notation
 

}


//#######################
    public void updateModel(){
       int ii, jj;                             // indices for loops
       
       Do_NPR(s_a);               // calculate numbers-per-recruit

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// NOT TRUSTING ANSWERS THAT GETTING FROM calc_Finit(), giving back whatever initial guess was for Finit - regardless of life history ///////////////////////
        calc_Finit();           // calculate Finit, given selectivity of hisorical aboriginal hunt (uniform 5+) and stable population size at start of trajectories (either 1846 or 1856)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Do_Rescale();           // Rescales NPR to total population size given initial depletion

        //System.out.println(Finit);
        //System.out.println(C_init);

        Assign_KK();            // Fill in vector of (possibly time varying) Carrying Capacity

        //System.out.println(Surv_const);

        Do_projections();       //

        xLowerBound = Double.valueOf(first_yr);      // this should not be hard coded either
        xUpperBound = Double.valueOf(last_yr);     // need to convert "dim1" to string.

        String N_tLabel = "Model abundance";
        if(Nplus != null){
            defaultXYDataset.addSeries(N_tLabel, Nplus);
        }else{
            //TODO
            //throw new SolexaGraphException("pointSeries was null : you should frist call the plotGraphDrawer.setSegmentSeriesFromFile(String fileAbsolutePath()");
        }

        String K_tLabel = "Model carrying capacity";
        if(K_t != null){
            defaultXYDataset.addSeries(K_tLabel, K_t);
        }else{
            //TODO
            //throw new SolexaGraphException ("pointSeries was null : you should frist call the plotGraphDrawer.setPointSeriesFromFile(String fileAbsolutePath()");
        }

        String matFemLabel = "Mature females";
        if(N_matfem_yr != null){
            defaultXYDataset.addSeries(matFemLabel, N_matfem_yr);
        }else{
            //TODO
            //throw new SolexaGraphException("pointSeries was null : you should frist call the plotGraphDrawer.setPointSeriesFromFile(String fileAbsolutePath()");
        }

        String matGenLabel = "Genetic bottleneck estimate";
        if(matFem_genetic != null){
            defaultXYDataset.addSeries(matGenLabel, matFem_genetic);
        }else{
            //TODO
            //throw new SolexaGraphException("pointSeries was null : you should frist call the plotGraphDrawer.setPointSeriesFromFile(String fileAbsolutePath()");
        }

        catchSeries_tmp = readCatchData();                               // catch data is hard-coded
        catchComm.removeSeries(0);
        catchComm.addSeries(catchSeries_tmp);

// ADD SOME ANNOTATED TEXT (FOR IMPLIED ABORIGINAL CATCH PRIOR TO 1846) TO THE CATCH PLOTS
//        XYTextAnnotation annotation = null;
//        Font font = new Font("SansSerif", Font.PLAIN, 14);
//        at = "Aboriginal removals prior to 1846: " + (int)C_init;
//        annotation = new XYTextAnnotation(at, first_yr, 1000.0);
//        annotation.setFont(font);
//        annotation.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
//        catchPlot.addAnnotation(annotation);

        readData();

    }

//#######################
    public void Do_NPR(double S_tmp1){
//******************************************************/
// Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations        
//******************************************************/
        int ii, jj, kk;
// INITIALIZE ////////////////////////////////////
        Male_age = new double [a_m+1];            // Vector with numbers-at-age for males
        Female_age = new double [a_m+1];          // "" females

        Sage = new double [a_m+1];                // Vector of survival-at-age
        Page = new double [a_m+1];                // Vector of rate of maturation-at-age
        prop_mat_a = new double[a_m+1];           // Proportion mature at age (all zeros until a_m)
        NPR_age = new double [a_m+1];             // Numbers-at-age-per-recruit
        Nage_imm_0 = new double [a_m+1];          //
        Nage_mat_0 = new double [a_m+1];          // numbers-at-age that are mature
        temp_1plus = 0.0;
        temp_mat = 0.0;
        b_eq = 0.0;                             // equilibrium birth rate at carrying capacity
        NPR_oneplus = 0.0;                      // size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing
//////////////////////////////////////////////////
        for (ii=0;ii<=a_m;ii++){              // calculate proportion mature at age.
            if(ii < a_m)
                prop_mat_a[ii] = 0.0;
            else
                prop_mat_a[ii] = 1.0;
        }

        for (jj=1;jj<=a_m;jj++)               // calculate rate of maturation, given proportion mature at age
         {
          if(1-prop_mat_a[jj]<0.001)
           Page[jj-1] = 1.0;
          else
           Page[jj-1] = (prop_mat_a[jj]-prop_mat_a[jj-1])/(1-prop_mat_a[jj-1]);
         }

        Sage[0] = S_tmp1 - delt_s;             // assign calf survival

        for(ii = 1; ii <= a_m; ii++)
            Sage[ii] = S_tmp1;                 // fill survival-at-age vector

        NPR_age[0]=0.5;                        // Numbers of females per recruit, assuming 50:50 sex ratio at birth
        Nage_imm_0[0]=0.5;                     // All calves assumed immature

       for(kk=1;kk<=a_m-1;kk++){
         NPR_age[kk] = NPR_age[kk-1]*Sage[kk-1];             //calculate numbers-at-age per recruit
         Nage_imm_0[kk] = NPR_age[kk]*(1-prop_mat_a[kk]);    //Initialize Numbers-at-Stage, based on maturity ogive
         Nage_mat_0[kk] = prop_mat_a[kk]*NPR_age[kk];
         temp_1plus += NPR_age[kk];                        // keep track of one-plus females per recruit
         temp_mat += Nage_mat_0[kk];                       // "" mature females per recruit
       }

      NPR_age[a_m] = NPR_age[a_m-1]*Sage[a_m-1]/(1-Sage[a_m-1]);   // PLUS GROUP numbers at age
      Nage_imm_0[a_m] = NPR_age[a_m]*(1-prop_mat_a[a_m]); 		 // by definition = 0.0 given assumption that age at transition to plus group equals age at maturity
      Nage_mat_0[a_m] = prop_mat_a[a_m]*NPR_age[a_m];		 // mature females

      temp_1plus = temp_1plus + NPR_age[a_m];									// keep track of one_plus component (females)
      temp_mat = temp_mat + Nage_mat_0[a_m];									// keep track of mature female component

      b_eq = 1.0 / (S_tmp1*(temp_mat - 1));                                  // equilibrium birth rate on a per recruit basis
      b_1 = b_eq + (b_max - b_eq) * (1 - Math.pow(depl_init,z));    // Birth rate in first year of projection, given intial depletion
      NPR_oneplus = temp_1plus;

    }

//#######################
    public void calc_Finit(){
//******************************************************/
// Using the Newton-Raphson method, find the root of a function known to lie in the interval [x1,x2].
//  The root (rtnewt) will be refined until its accuracty is known within +/- (xacc). (funcd) is a user
//  supplied subrouting that retuns both the function value and the first derivative of the function at the point (x).'
//  - description from Numerical Recipes book
//  In this case going to be solving for hunting mortality during 1900-1930 which would have left the population at the same size over those years
//******************************************************/
        int jj,jmax; 	// maximum number of iterations
        int debug;
// INITIALIZE ////////////////////////////////////
        f1_newt = 0.0;
        df_newt = 0.0;
        dx_newt = 0.0;
        Finit = 0.0;
        objf_newt = 0.0;

        jmax = 20;
        delta_Finit=0.000001;
        //delta_Finit=0.0001;
        rtnewt=0.05;                //initial guess for gray whale Finit (hunting mortality at start of commercial whaling)
//////////////////////////////////////////////////

        for(jj=1;jj<=jmax;jj++)
           {
        	Finit = rtnewt;
        	Initial_Fishing();	// function call
        	f1_newt = objf_newt;
        	Finit = rtnewt+delta_Finit;
        	if(Finit>1.0) Finit=1.0;
                if(Finit<0.0) Finit=0.0;
                Initial_Fishing();	// function call
        	df_newt = (objf_newt-f1_newt)/delta_Finit;
        	dx_newt = f1_newt/df_newt;
        	rtnewt = rtnewt-dx_newt;
           }

        debug = 1;
        // System.out.println(Finit);
    }

//#######################
    public void Initial_Fishing(){
//******************************************************/
// given numbers per-recruit at equilibrium - apply fishing mortality (Finit)
//******************************************************/
         int ii,jj;
         int age_recruitment;   // assumed age at recruitment into the historical aboriginal fishery (5+, might make this 1+ though?)
         double rec_init;	//initial recruitment conditioned on Finit (in terms of females)
         double temp_mature;
         double pred_rec;
// INITIALIZE ////////////////////////////////////
         age_recruitment = 1;   // assumed age at recruitment into the historical aboriginal fishery (1+)
         temp_mature = 0.0;
         temp_1plus = 0.0;
         rec_init = 0.0;
         pred_rec = 0.0;
         initial_oneplus = 0.0;
         Nage = new double [a_m+1];           // numbers-at-age (both sexes combined)
         select_init = new double[a_m+1];     // Initial selectivity associated w/historical aboriginal hunt (assumed uniform on 5+)

// THIS IS CAUSING TROUBLE : SEEMS TO BE ASSIGNING OBJECT Nage to NPR_age, such that subsequent operations on Nage are changing values of NPR_age!!
//         Nage = NPR_age;            //initialize numbers at age given number per recruit (working with females only at this point)


         //////////////////////////////////////////////////
         for(ii=0;ii<=a_m;ii++){
            Nage[ii]=NPR_age[ii];
             if(ii < age_recruitment)
                select_init[ii]=0.0;
            else
                select_init[ii]=1.0;
         }

         for(ii=1;ii<=a_m-1;ii++)
          {
           Nage[ii]=NPR_age[ii-1]*Sage[ii-1];                   // natural mortality - Sage vector is initialized in the NPR function
           Nage[ii]=Nage[ii]*(1-Finit*select_init[ii-1]);	// fishing mortality at start of time series
           temp_mature += prop_mat_a[ii]*Nage[ii];
           temp_1plus += Nage[ii];
          }

          Nage[a_m] = NPR_age[a_m-1]*Sage[a_m-1];	// plus group per recruit
          Nage[a_m] = Nage[a_m]*(1-Finit*select_init[a_m-1]);
          Nage[a_m] = Nage[a_m]/(1-Sage[a_m]*(1-Finit*select_init[a_m]));
          temp_mature += prop_mat_a[a_m]*Nage[a_m];
          temp_1plus += Nage[a_m];

          rec_init = depl_init * NPR_oneplus / temp_1plus; 	// now, rescale to get the desired 1+ depletion conditioned on Finit
          rec_init = rec_init / 2;				// divide by two to get female recruits conditioned on Finit

// calculate actual numbers per recruit conditioned on Finit
          Nage[0]=rec_init;
          for(jj=1;jj<=a_m-1;jj++)
           {
            Nage[jj] = Nage[jj-1]*Sage[jj-1];		// natural mortality
            Nage[jj] = Nage[jj]*(1-Finit*select_init[jj-1]);	// fishing mortality at start of time series
            initial_oneplus += Nage[jj];
           }

          Nage[a_m]=Nage[a_m-1]*Sage[a_m-1];                 // plus group per recruit
          Nage[a_m]=Nage[a_m]*(1-Finit*select_init[a_m-1]);
          Nage[a_m]=Nage[a_m]/(1-Sage[a_m]*(1-Finit*select_init[a_m]));
          initial_oneplus += Nage[a_m];

          pred_rec = b_1*(temp_mature-1)-1;
          objf_newt = pred_rec*pred_rec;

    }

//#######################
    public void Do_Rescale(){
////#######################
//// Rescale initial numbers at age to 1+ population size given depletion_30
////#######################
    double scale_pop;
    int aa, ii;
// INITIALIZE ////////////////////////////////////
    scale_pop = 0.0;
    aa = 0;
//    Nplus = new double[2][dim1+1];                // Vector of total age 1+ population
    N_m_yr_1 = new double[2][dim1+1];
    N_calf = new double[2][dim1+1];
    N_calvn_yr_1 = new double[2][dim1+1];
    N_recpt_yr_1 = new double[2][dim1+1];
    N_imm_yr_1 = new double[2][dim1+1];
    Ntot = new double[2][dim1+1];
    NAll = new double[5][dim1+1][a_m+1];               // Vector of stages: (1) immature fem (2) mature receptive fem (3) mature calving fem, and (4) males

    Nplus[1][0] = 0.0;
//    N_imm_yr_1[1][0] = 0.0;
//    N_recpt_yr_1[1][0] = 0.0;
//    N_calvn_yr_1[1][0] = 0.0;
//    N_m_yr_1[1][0] = 0.0;
//    Ntot[1][0] = 0.0;
//    N_calf[1][0] = 0.0;
//////////////////////////////////////////////////
   scale_pop = 0.5 * KK * depl_init / initial_oneplus;
   
   for (ii=0;ii<=a_m;ii++)                      // rescale female numbers at age per recruit to get desired initial depletion
        Nage[ii] = Nage[ii] * scale_pop;	// this now in terms of females - so, just assign another vector equal to this one to initialize males

  for(aa=0;aa<=a_m;aa++) {
    NAll[1][0][aa] = Nage[aa] * (1-prop_mat_a[aa]); 	// immature numbers at age for females
    NAll[2][0][aa] = Nage[aa] * prop_mat_a[aa]*(1-b_1); // receptive numbers at age for females
    NAll[3][0][aa] = Nage[aa] * prop_mat_a[aa]*b_1; 	// calving numbers at age for females
    NAll[4][0][aa] = Nage[aa];   			// intial numbers at age for males
    N_imm_yr_1[1][0] += N_imm_yr_1[1][0] + NAll[1][0][aa];		// total numbers at stage
    N_recpt_yr_1[1][0] += NAll[2][0][aa];
    N_calvn_yr_1[1][0] += NAll[3][0][aa];
    N_m_yr_1[1][0]     += NAll[4][0][aa];
    Ntot[1][0]         += NAll[1][0][aa] + NAll[2][0][aa] + NAll[3][0][aa] + NAll[4][0][aa];
  }

   depl[1][0] = depl_init;
   N_init = KK * depl_init;                     //treating initial depletion as estimated parameter, so multiply by K to get 1+ population size in 1930
   N_calf[1][0] = NAll[1][0][0] * 2;            // Initial number of calves (for output)
   Nplus[1][0] = Ntot[1][0] - N_calf[1][0];

// Calculate aboriginal catch in 1845
   C_init = Finit*Nplus[1][0];
}

//#######################
    public void addSliders (Container content) {
        int gridx, gridy;
        int xpos, ypos;
        int st,sl,sr,sb;        // padding in pixels around the scroll bars
        int lt,ll,lr,lb;        // padding in pixels around the labels above the scroll bars        
        int tt,tl,tr,tb;        // padding in pixels around the text boxes displaying parameter values
        int debug;
//        Border blackline, empty;

        int initialDelay = ToolTipManager.sharedInstance().getInitialDelay();   // Get current delay before showing tool tip
        int dismissDelay = ToolTipManager.sharedInstance().getDismissDelay();   // Get current delay before hiding tool tip

        ToolTipManager.sharedInstance().setInitialDelay(0);                     // Show tool tips immediately
        dismissDelay = Integer.MAX_VALUE;
        ToolTipManager.sharedInstance().setDismissDelay(dismissDelay);          // Keep the tool tip showing
        UIManager.put("ToolTip.background", new ColorUIResource(255, 247, 200)); // The color is #fff7c8.
        UIManager.put("TabbedPane.selected", Color.yellow);

        JTabbedPane tabbedPane1 = new JTabbedPane();
        tabbedPane1.setPreferredSize(new Dimension(500,400));   // SIZING OF PANELS

//        empty = BorderFactory.createEmptyBorder();
//        tabbedPane1.setBorder(empty);

        JPanel panel_2 = new JPanel();
        panel_2.setPreferredSize(new Dimension(500,400));   // SIZING OF PANELS

//        blackline = BorderFactory.createLineBorder(Color.black);
//        panel_2.setBorder(blackline);

        JPanel panel_3 = new JPanel();
        panel_3.setPreferredSize(new Dimension(500,400));   // SIZING OF PANELS

        JPanel panel_4 = new JPanel();
        panel_3.setPreferredSize(new Dimension(500,400));   // SIZING OF PANELS

        GridBagLayout gbl = new GridBagLayout();
        GridBagConstraints constraints = new GridBagConstraints();

        JSlider scroll_1 = new JSlider(JSlider.HORIZONTAL, 75, 100, (int)(depl_init*100));  // Depletion initial year (currently 1846)
        JSlider scroll_2 = new JSlider(JSlider.HORIZONTAL, 10000, 50000, (int)KK);          // (K) Carrying Capacity
        JSlider scroll_3 = new JSlider(JSlider.HORIZONTAL, 995, 1010, 1000);                // K multiplier
        JSlider scroll_4 = new JSlider(JSlider.HORIZONTAL, 0, 40, (int)(depl_star*100));    // Depensation (depletion level below which depensation occurs)
        JSlider scroll_5 = new JSlider(JSlider.HORIZONTAL, 0, 100, (int)(phi_fem*100));     // Sex ratio of historic commercial catch (i.e. 1.0 = 100%females)
        JSlider scroll_6 = new JSlider(JSlider.HORIZONTAL, 0, 100, (int)(beta_cow*100));    // Proportion female commercial catch: cows w/calves (i.e. 1.0 = 100% of females = cows w/calves)
        final JSlider scroll_7 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_a*1000));      // 1+ survival (constant or pre-1900 as dictated by radio button)
        final JSlider scroll_77 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_post*1000));      // 1+ survival post-1900
        final JSlider scroll_777 = new JSlider(JSlider.HORIZONTAL, 950, 999, (int)(s_a*1000));      // 1+ survival post-1900
        JSlider scroll_8 = new JSlider(JSlider.HORIZONTAL, 1, 500, (int)(delt_s*1000));     // Delta survival, difference between adult and calf survival
        JSlider scroll_9 = new JSlider(JSlider.HORIZONTAL, 6, 12, a_m);                     // Age at maturity
        JSlider scroll_10 = new JSlider(JSlider.HORIZONTAL, 100, 1000, (int)(z*100));       // Density dependent shape parameter (Z)

        JSlider scroll_11 = new JSlider(JSlider.HORIZONTAL, -35, 35, (int)(eps_Mex*100));   // Mex catch magnitude unertainty multipler for landings
        JSlider scroll_12 = new JSlider(JSlider.HORIZONTAL, 110, 152, (int)(LRF_Mex*100));  // LRF for Mex Catches

        JSlider scroll_13 = new JSlider(JSlider.HORIZONTAL, 120, 210, (int)(LRF_CA*100));   // LRF for Mex Catches
        JSlider scroll_14 = new JSlider(JSlider.HORIZONTAL, -52, 52, (int)(eps_CA*100));    // CA catch magnitude unertainty multipler for landings

        JSlider scroll_15 = new JSlider(JSlider.HORIZONTAL, 0, 35, (int)(CV_Add*100));      // CV_add for abundance estimates

        JSlider scroll_16 = new JSlider(JSlider.HORIZONTAL, 20, 100, (int)(S_star*100));      // Survival rate multiplier during mortality event of '99 and '00

        JSlider scroll_17 = new JSlider(JSlider.HORIZONTAL, 50, 200, (int)(C_mult*100));      // Catch (removals) multiplier ; generic and applying it to all fisheries prior to 1900
        final JSlider scroll_18 = new JSlider(JSlider.HORIZONTAL, 0, 2000, (int)(C_init));      // Catch (removals) multiplier ; generic and applying it to all fisheries prior to 1900

        panel_2.setLayout(gbl);
        panel_3.setLayout(gbl);
        panel_4.setLayout(gbl);

////////////////////////////////////////////////////////////////////////////////
// Sizing and spacing for first column of components
////////////////////////////////////////////////////////////////////////////////
        st = 0; sl = 0; sr = 10; sb = 5;         // inset in pixels around the scroll bars//
        lt = 5; ll = 0; lr = 0; lb = 0;         // inset in pixels around the labels above the scroll bars
        tt = 5; tl = 0; tr = 0; tb = 5;         // inset in pixels around the text boxes displaying parameter values

        constraints.weightx = 0.33;
        constraints.fill = GridBagConstraints.HORIZONTAL;
////////////////////////////////////////////////////////////////////////////////
// Show implied pre-1846 aboriginal annual removals (w/disabled slider)
////////////////////////////////////////////////////////////////////////////////
        xpos = 2; ypos = 7;             // text box pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.anchor = GridBagConstraints.NORTH;
        constraints.fill = GridBagConstraints.NONE;
        constraints.insets = new Insets(tt,tl,tb,tr);
        constraints.gridheight = 1;
        final JTextField jTextField18 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField18.setText(String.valueOf((int)C_init));
        jTextField18.setEnabled(true);
        gbl.setConstraints(jTextField18, constraints);
        panel_3.add(jTextField18);

        xpos = 2; ypos = 6;             // Label pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        final JLabel label_18 = new JLabel("Implied Pre-1846 Aboriginal Annual Removal");    // test label while debugging the layout
        label_18.setToolTipText("The size of the pre-1846 aboriginal annual removals implied by the fraction of carrying capacity in 1846 and life history parameter values.");
        label_18.setEnabled(true);
        gbl.setConstraints(label_18, constraints);
        panel_3.add(label_18);             // see if can add a label here.

        xpos = 6; ypos = 13;             // Scroll pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_18, constraints);
        scroll_18.setMajorTickSpacing(500);
        scroll_18.setPaintLabels(true);
        scroll_18.setPaintTicks(true);
        scroll_18.setEnabled(false);
//        panel_3.add(scroll_18);

//         scroll_18.addChangeListener(new ChangeListener() {
//            public void stateChanged(ChangeEvent evt) {
//                JSlider scroll_18 = (JSlider) evt.getSource();
//                    scroll_18.setValue((int)C_init);
//                    jTextField18.setText(String.valueOf(C_init));
////                    updateModel();                               // updates model trajectory(s)
//            }
//         });

////////////////////////////////////////////////////////////////////////////////
// CARRYING CAPACITY
////////////////////////////////////////////////////////////////////////////////
        gridx = 0; gridy = 0;           // row / column grid layout of components

        xpos = gridx+1; ypos = gridy+1;             // text box pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(tt,tl,tb,tr);
//        constraints.fill = GridBagConstraints.HORIZONTAL;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
        final JTextField jTextField2 = new JTextField(5);   // needs to be declared as 'FINAL'
        jTextField2.setText(String.valueOf(KK));
        gbl.setConstraints(jTextField2, constraints);
        panel_2.add(jTextField2);

        xpos = gridx; ypos = gridy;             // Label pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_2 = new JLabel("Carrying Capacity (K)");    // test label while debugging the layout
        label_2.setToolTipText("Maximum number of animals that can be sustained by their habitat.");
        gbl.setConstraints(label_2, constraints);
        panel_2.add(label_2);             // see if can add a label here.

        xpos = gridx; ypos = gridy+1;             // Sroll pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        scroll_2.setMajorTickSpacing(10000);
        scroll_2.setPaintLabels(true);
        scroll_2.setPaintTicks(true);
        gbl.setConstraints(scroll_2, constraints);
        panel_2.add(scroll_2);

        scroll_2.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_2 = (JSlider) evt.getSource();
                    if(KK==0) KK = 1;                       // Carrying Capacity of 0 makes no sense (not that K=1 does)
                    KK = scroll_2.getValue();
                    jTextField2.setText(String.valueOf((int)KK));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
           }
        });

////////////////////////////////////////////////////////////////////////////////
// INITIAL POPULATION SIZE OR DEPLETION
////////////////////////////////////////////////////////////////////////////////
        gridx = 0; gridy = 0;           // row / column grid layout of components

        xpos = gridx+1; ypos = gridy+1;             // text box pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField1 = new JTextField(4);   // needs to be declared as 'FINAL'
        jTextField1.setText(String.valueOf(depl_init));
        gbl.setConstraints(jTextField1, constraints);
        panel_3.add(jTextField1);

        xpos = gridx; ypos = gridy;             // Label pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_1 = new JLabel("Fraction of Carrying Capacity in 1846");    // test label while debugging the layout
        label_1.setToolTipText("The fraction of carrying capacity at the start of commercial whaling, resulting from depletion due to aboriginal whaling prior to 1846.");
        gbl.setConstraints(label_1, constraints);
        panel_3.add(label_1);             // see if can add a label here.

        Hashtable labelTable = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable.put( new Integer( 75 ), new JLabel("0.75") );
        labelTable.put( new Integer( 80 ), new JLabel("0.80") );
        labelTable.put( new Integer( 85 ), new JLabel("0.85") );
        labelTable.put( new Integer( 90 ), new JLabel("0.90") );
        labelTable.put( new Integer( 95 ), new JLabel("0.95") );
        labelTable.put( new Integer( 100 ), new JLabel("1.00") );
        scroll_1.setLabelTable( labelTable );

        xpos = gridx; ypos = gridy+1;             // Scroll pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_1, constraints);
        scroll_1.setMajorTickSpacing(5);
        scroll_1.setMinorTickSpacing(5);
        scroll_1.setPaintLabels(true);
        scroll_1.setPaintTicks(true);
        panel_3.add(scroll_1);

        scroll_1.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_1 = (JSlider) evt.getSource();
                    depl_init = (float)scroll_1.getValue()/100;
                    jTextField1.setText(String.valueOf(String.format("%.2g%n", depl_init)));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// CARRYING CAPACITY MULTIPLIER
////////////////////////////////////////////////////////////////////////////////

        xpos = 3; ypos = 1;             // text box 1 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField3 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField3.setText(String.valueOf(K_mult));
        gbl.setConstraints(jTextField3, constraints);
        panel_3.add(jTextField3);

        xpos = 2; ypos = 0;             // Label 1 pos.
//        xpos = 0; ypos = 6;             // Label 3 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_3 = new JLabel("Carrying Capacity Multiplier");    // test label while debugging the layout
        label_3.setToolTipText("Allows for a linear annual rate of change in carrying capacity, perhpas related to some trend in environmental conditions through time.");
        gbl.setConstraints(label_3, constraints);
        panel_3.add(label_3);             // see if can add a label here.

        Hashtable labelTable_Kmult = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_Kmult.put( new Integer( 995 ), new JLabel("0.995") );
        labelTable_Kmult.put( new Integer( 1000 ), new JLabel("1.00") );
        labelTable_Kmult.put( new Integer( 1005 ), new JLabel("1.05") );
        labelTable_Kmult.put( new Integer( 1010 ), new JLabel("1.10") );
        scroll_3.setLabelTable(labelTable_Kmult);

        xpos = 2; ypos = 1;             // Sroll 1 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_3, constraints);
        scroll_3.setMajorTickSpacing(10);
        scroll_3.setMinorTickSpacing(5);
        scroll_3.setPaintLabels(true);
        scroll_3.setPaintTicks(true);
        panel_3.add(scroll_3);

        scroll_3.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_3 = (JSlider) evt.getSource();
                    K_mult = scroll_3.getValue();
                    K_mult = K_mult / 1000;
                    jTextField3.setText(String.valueOf(K_mult));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// DENSITY DEPENDENCE SHAPE PARAMETER (Z)
////////////////////////////////////////////////////////////////////////////////
//        gridx = 0; gridy = 3;           // row / column grid layout of components

         xpos = 1; ypos = 4;             // text box 4 pos.
//        xpos = gridx+1; ypos = gridy+1;             // text box pos.
         constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField10 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField10.setText(String.valueOf(z));
        gbl.setConstraints(jTextField10, constraints);
        panel_2.add(jTextField10);

//        xpos = gridx; ypos = gridy;             // Label pos.
        xpos = 0; ypos = 3;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_10 = new JLabel("Density Dependence Shape Parameter (Z)");    // test label while debugging the layout
        label_10.setToolTipText("Governs the shape of the production curve (rate of change vs. population size relative to K): Z = 1.0 corresponds approximately to maximum productivity at 50% of K; Z = 2.39 at ~60% of K etc.");
        gbl.setConstraints(label_10, constraints);
        panel_2.add(label_10);             // see if can add a label here.

        Hashtable labelTable_z = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_z.put( new Integer( 200 ), new JLabel("2.0") );
        labelTable_z.put( new Integer( 400 ), new JLabel("4.0") );
        labelTable_z.put( new Integer( 600 ), new JLabel("6.0") );
        labelTable_z.put( new Integer( 800 ), new JLabel("8.0") );
        labelTable_z.put( new Integer( 1000 ), new JLabel("10.0") );
        scroll_10.setLabelTable(labelTable_z);

//        xpos = gridx; ypos = gridy;             // Scroll pos.
        xpos = 0; ypos = 4;             // Scroll 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_10, constraints);
        scroll_10.setMajorTickSpacing(100);
//        scroll_10.setMinorTickSpacing(5000);
        scroll_10.setPaintLabels(true);
        scroll_10.setPaintTicks(true);
        panel_2.add(scroll_10);

         scroll_10.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_10 = (JSlider) evt.getSource();
                    z = scroll_10.getValue();
                    z = z / 100;
                    jTextField10.setText(String.valueOf(z));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// eps_CA
////////////////////////////////////////////////////////////////////////////////
        xpos = 1; ypos = 7;             // text box 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField14 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField14.setText(String.valueOf(eps_CA));
        gbl.setConstraints(jTextField14, constraints);
        panel_2.add(jTextField14);

        xpos = 0; ypos = 6;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_14 = new JLabel("Log-Uncertainty CA Shore landings");    // test label while debugging the layout
        label_14.setToolTipText("Uncertainty (in log-space) of the estimated landings for shore-based whaling off California (1854-1899). Corresponds to the average annual standard error estimated for these landings.");
        gbl.setConstraints(label_14, constraints);
        panel_2.add(label_14);             // see if can add a label here.

        Hashtable labelTable_eps_CA = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_eps_CA.put( new Integer( -50 ), new JLabel("-0.5") );
        labelTable_eps_CA.put( new Integer( -25 ), new JLabel("-0.25") );
        labelTable_eps_CA.put( new Integer( 0 ), new JLabel("0.0") );
        labelTable_eps_CA.put( new Integer( 25 ), new JLabel("0.25") );
        labelTable_eps_CA.put( new Integer( 50 ), new JLabel("0.5") );
        scroll_14.setLabelTable(labelTable_eps_CA);

        xpos = 0; ypos = 7;             // Scroll 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_14, constraints);
        scroll_14.setMajorTickSpacing(10);
        scroll_14.setPaintLabels(true);
        scroll_14.setPaintTicks(true);
        panel_2.add(scroll_14);

         scroll_14.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_14 = (JSlider) evt.getSource();
                    eps_CA = scroll_14.getValue();
                    eps_CA = eps_CA / 100;
                    jTextField14.setText(String.valueOf(eps_CA));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// LRF_CA: uncertainty in magnitude of removals from California shore whaling
////////////////////////////////////////////////////////////////////////////////
        xpos = 1; ypos = 10;                                 // text box 7 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField13 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField13.setText(String.valueOf(LRF_CA));
        gbl.setConstraints(jTextField13, constraints);
        panel_2.add(jTextField13);

        xpos = 0; ypos = 9;                                 // Label 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_13 = new JLabel("Loss Rate Factor: CA Shore");    // test label while debugging the layout
        label_13.setToolTipText("The Loss Rate Factor (animals which were struck and died, but not landed) for shore-based whaling off California (1854-1899).");
        gbl.setConstraints(label_13, constraints);
        panel_2.add(label_13);             // see if can add a label here.

        Hashtable labelTable_LRF_CA = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_LRF_CA.put( new Integer( 120), new JLabel("1.20") );
        labelTable_LRF_CA.put( new Integer( 140 ), new JLabel("1.40") );
        labelTable_LRF_CA.put( new Integer( 160 ), new JLabel("1.60") );
        labelTable_LRF_CA.put( new Integer( 180 ), new JLabel("1.80") );
        labelTable_LRF_CA.put( new Integer( 200 ), new JLabel("2.00") );
        scroll_13.setLabelTable(labelTable_LRF_CA);

        xpos = 0; ypos = 10;             // Scroll 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_13, constraints);
        scroll_13.setMajorTickSpacing(10);
        scroll_13.setPaintLabels(true);
        scroll_13.setPaintTicks(true);
        panel_2.add(scroll_13);

         scroll_13.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_13 = (JSlider) evt.getSource();
                    LRF_CA = scroll_13.getValue();
                    LRF_CA = LRF_CA / 100;
                    jTextField13.setText(String.valueOf(LRF_CA));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// NEW COLUMN OF CONTROLS, SO TWEAK SPACING A BIT TO GET SOME SEPERATION BETWEEN SLIDERS
////////////////////////////////////////////////////////////////////////////////
//        st = 0; sl = 0; sr = 0; sb = 5;        // padding in pixels around the scroll bars//
//        lt = 5; ll = 0; lr = 0; lb = 5;         // padding in pixels around the labels above the scroll bars
//        tt = 5; tl = 0; tr = 0; tb = 5;        // padding in pixels around the text boxes displaying parameter values
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// SEX RATIO IN LAGOON CATCH
////////////////////////////////////////////////////////////////////////////////
        xpos = 3; ypos = 1;             // text box 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField5 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField5.setText(String.valueOf(phi_fem));
        gbl.setConstraints(jTextField5, constraints);
        panel_2.add(jTextField5);

        xpos = 2; ypos = 0;             // Label 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_5 = new JLabel("Fraction Female in Lagoon Removals");    // test label while debugging the layout
        label_5.setToolTipText("The fraction of the ship-based lagoon removals in Mexico (1846-1873) that were female.");
        gbl.setConstraints(label_5, constraints);
        panel_2.add(label_5);             // see if can add a label here.

        Hashtable labelTable_phi_fem = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_phi_fem.put( new Integer( 0 ), new JLabel("0.00") );
        labelTable_phi_fem.put( new Integer( 25 ), new JLabel("0.25") );
        labelTable_phi_fem.put( new Integer( 50 ), new JLabel("0.50") );
        labelTable_phi_fem.put( new Integer( 75 ), new JLabel("0.75") );
        labelTable_phi_fem.put( new Integer( 100 ), new JLabel("1.00") );
        scroll_5.setLabelTable(labelTable_phi_fem);

        xpos = 2; ypos = 1;             // Scroll 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_5, constraints);
        scroll_5.setMajorTickSpacing(25);
        scroll_5.setPaintLabels(true);
        scroll_5.setPaintTicks(true);
        panel_2.add(scroll_5);

         scroll_5.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_5 = (JSlider) evt.getSource();
                    phi_fem = scroll_5.getValue();
                    phi_fem = phi_fem / 100;
                    jTextField5.setText(String.valueOf(phi_fem));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// PROPORTION (MATURE) FEMALES IN LAGOON CATCH: COWS W/CALVES
////////////////////////////////////////////////////////////////////////////////
        xpos = 3; ypos = 4;             // text box 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField6 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField6.setText(String.valueOf(beta_cow));
        gbl.setConstraints(jTextField6, constraints);
        panel_2.add(jTextField6);

        xpos = 2; ypos = 3;             // Label 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_6 = new JLabel("Proportion Female Lagoon Removals: Cows w/Calves");    // test label while debugging the layout
        label_6.setToolTipText("The proportion of the female removals from the lagoons that were cows with calves. 100% calf mortality is assumed if a cow w/calf is killed.");
        gbl.setConstraints(label_6, constraints);
        panel_2.add(label_6);             // see if can add a label here.

        Hashtable labelTable_beta_cow = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_beta_cow.put( new Integer( 0 ), new JLabel("0.00") );
        labelTable_beta_cow.put( new Integer( 25 ), new JLabel("0.25") );
        labelTable_beta_cow.put( new Integer( 50 ), new JLabel("0.50") );
        labelTable_beta_cow.put( new Integer( 75 ), new JLabel("0.75") );
        labelTable_beta_cow.put( new Integer( 100 ), new JLabel("1.00") );
        scroll_6.setLabelTable(labelTable_beta_cow);

        xpos = 2; ypos = 4;             // Scroll 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_6, constraints);
        scroll_6.setMajorTickSpacing(25);
        scroll_6.setPaintLabels(true);
        scroll_6.setPaintTicks(true);
        panel_2.add(scroll_6);

         scroll_6.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_6 = (JSlider) evt.getSource();
                    beta_cow = scroll_6.getValue();
                    beta_cow = beta_cow / 100;
                    jTextField6.setText(String.valueOf(beta_cow));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// eps_Mex: uncertainty in magnitude of ship based landings in Mexican waters
////////////////////////////////////////////////////////////////////////////////
        xpos = 3; ypos = 7;                                 // text box 7 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField11 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField11.setText(String.valueOf(eps_Mex));
        gbl.setConstraints(jTextField11, constraints);
        panel_2.add(jTextField11);

        xpos = 2; ypos = 6;                                 // Label 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_11 = new JLabel("Log-Uncertainty of Mex. landings");    // test label while debugging the layout
        label_11.setToolTipText("Uncertainty (in log-space) of the estimated landings for ship-based landings in Mexico (1846-1873). Corresponds to the average annual standard error estimated for these landings.");
        gbl.setConstraints(label_11, constraints);
        panel_2.add(label_11);             // see if can add a label here.

        Hashtable labelTable_eps_Mex = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_eps_Mex.put( new Integer( -35 ), new JLabel("-0.30") );
        labelTable_eps_Mex.put( new Integer( -15 ), new JLabel("-0.15") );
        labelTable_eps_Mex.put( new Integer( 0 ), new JLabel("0.00") );
        labelTable_eps_Mex.put( new Integer( 15 ), new JLabel("0.15") );
        labelTable_eps_Mex.put( new Integer( 30 ), new JLabel("0.30") );
        scroll_11.setLabelTable(labelTable_eps_Mex);

        xpos = 2; ypos = 7;             // Scroll 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_11, constraints);
        scroll_11.setMajorTickSpacing(5);
        scroll_11.setPaintLabels(true);
        scroll_11.setPaintTicks(true);
        panel_2.add(scroll_11);

         scroll_11.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_11 = (JSlider) evt.getSource();
                    //Math.log(z);
                    eps_Mex = scroll_11.getValue();
                    eps_Mex = eps_Mex / 100;
                    jTextField11.setText(String.valueOf(eps_Mex));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);

            }
         });
////////////////////////////////////////////////////////////////////////////////
// LRF_Mex: uncertainty in magnitude of ship based landings in Mexican waters
////////////////////////////////////////////////////////////////////////////////
        xpos = 3; ypos = 10;                                 // text box 7 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField12 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField12.setText(String.valueOf(LRF_Mex));
        gbl.setConstraints(jTextField12, constraints);
        panel_2.add(jTextField12);

        xpos = 2; ypos = 9;                                 // Label 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_12 = new JLabel("Loss Rate Factor: Mexico");    // test label while debugging the layout
        label_12.setToolTipText("The Loss Rate Factor (animals which were struck and died, but not landed) for ship-based whaling in Mexico (1846-1873).");
        gbl.setConstraints(label_12, constraints);
        panel_2.add(label_12);             // see if can add a label here.

        Hashtable labelTable_LRF_Mex = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_LRF_Mex.put( new Integer( 110), new JLabel("1.10") );
        labelTable_LRF_Mex.put( new Integer( 120 ), new JLabel("1.20") );
        labelTable_LRF_Mex.put( new Integer( 130 ), new JLabel("1.30") );
        labelTable_LRF_Mex.put( new Integer( 140 ), new JLabel("1.40") );
        labelTable_LRF_Mex.put( new Integer( 150 ), new JLabel("1.50") );
        scroll_12.setLabelTable(labelTable_LRF_Mex);

        xpos = 2; ypos = 10;             // Scroll 5 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_12, constraints);
        scroll_12.setMajorTickSpacing(10);
        scroll_12.setPaintLabels(true);
        scroll_12.setPaintTicks(true);
        panel_2.add(scroll_12);

         scroll_12.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_12 = (JSlider) evt.getSource();
                    LRF_Mex = scroll_12.getValue();
                    LRF_Mex = LRF_Mex / 100;
                    jTextField12.setText(String.valueOf(LRF_Mex));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });


////////////////////////////////////////////////////////////////////////////////
// DEPENSATION PARAMETER
////////////////////////////////////////////////////////////////////////////////
        xpos = 1; ypos = 4;             // text box 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField4 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField4.setText(String.valueOf(depl_star));
        gbl.setConstraints(jTextField4, constraints);
        panel_3.add(jTextField4);

        xpos = 0; ypos = 3;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_4 = new JLabel("Depletion Level for Depensation");    // test label while debugging the layout
        label_4.setToolTipText("The abundance relative to carrying capacity, at which the population experiences a hypothetical decline in birth rates, perhaps due to the disruptive nature of whaling on the calving grounds.");
        gbl.setConstraints(label_4, constraints);
        panel_3.add(label_4);             // see if can add a label here.

        Hashtable labelTable_depl_star = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_depl_star.put( new Integer( 0 ), new JLabel("0.0") );
        labelTable_depl_star.put( new Integer( 10 ), new JLabel("0.1") );
        labelTable_depl_star.put( new Integer( 20 ), new JLabel("0.2") );
        labelTable_depl_star.put( new Integer( 30 ), new JLabel("0.3") );
        labelTable_depl_star.put( new Integer( 40 ), new JLabel("0.4") );
        scroll_4.setLabelTable(labelTable_depl_star);

        xpos = 0; ypos = 4;             // Scroll 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_4, constraints);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        scroll_4.setMajorTickSpacing(10);
//        scroll_4.setMinorTickSpacing(5000);
        scroll_4.setPaintLabels(true);
        scroll_4.setPaintTicks(true);
        panel_3.add(scroll_4);

         scroll_4.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_4 = (JSlider) evt.getSource();
                    depl_star = scroll_4.getValue();
                    depl_star = depl_star / 100;
                    jTextField4.setText(String.valueOf(depl_star));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// NEW COLUMN OF CONTROLS, SO TWEAK SPACING A BIT TO GET SOME SEPERATION BETWEEN SLIDERS
////////////////////////////////////////////////////////////////////////////////
//        st = 0; sl = 0; sr = 0; sb = 5;        // padding in pixels around the scroll bars//
//        lt = 5; ll = 0; lr = 0; lb = 5;         // padding in pixels around the labels above the scroll bars
//        tt = 5; tl = 0; tr = 0; tb = 5;        // padding in pixels around the text boxes displaying parameter values
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//// S_star: Survival rate multiplier during mortality event
//////////////////////////////////////////////////////////////////////////////////
        xpos = 5; ypos = 10;             // text box
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField16 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField16.setText(String.valueOf(S_star));
        gbl.setConstraints(jTextField16, constraints);
        panel_2.add(jTextField16);

        xpos = 4; ypos = 9;             // Label .
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_16 = new JLabel("Mortality Event Survival Rate Multplier");    // test label while debugging the layout
        label_16.setToolTipText("Survival rate multplier during mortality event of 1999 and 2000. A value of 1.0 indicates 100% of normal survival (i.e. no additional mortality)");
        gbl.setConstraints(label_16, constraints);
        panel_2.add(label_16);             // see if can add a label here.

        Hashtable labelTable_S_star = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_S_star.put( new Integer( 20 ), new JLabel("0.2") );
        labelTable_S_star.put( new Integer( 40 ), new JLabel("0.4") );
        labelTable_S_star.put( new Integer( 60 ), new JLabel("0.6") );
        labelTable_S_star.put( new Integer( 80 ), new JLabel("0.8") );
        labelTable_S_star.put( new Integer( 100 ), new JLabel("1.0") );
        scroll_16.setLabelTable(labelTable_S_star);

        xpos = 4; ypos = 10;             // Scroll
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_16, constraints);
        scroll_16.setMajorTickSpacing(20);
        scroll_16.setPaintLabels(true);
        scroll_16.setPaintTicks(true);
        panel_2.add(scroll_16);

         scroll_16.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_16 = (JSlider) evt.getSource();
                    S_star = scroll_16.getValue();
                    S_star = S_star / 100;
                    jTextField16.setText(String.valueOf(S_star));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// DELTA SURVIVAL RATE (DIFFERENCE BETWEEN 1+ AND CALF SURVIVAL RATES)
////////////////////////////////////////////////////////////////////////////////
        xpos = 5; ypos = 4;             // text box 
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField8 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField8.setText(String.valueOf(delt_s));
        gbl.setConstraints(jTextField8, constraints);
        panel_2.add(jTextField8);

        xpos = 4; ypos = 3;             // Label .
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_8 = new JLabel("Difference Between 1+ and Calf Survival");    // test label while debugging the layout
        label_8.setToolTipText("The difference between the survival rate of calves and the survival rate of animals age one year and older.");
        gbl.setConstraints(label_8, constraints);
        panel_2.add(label_8);             // see if can add a label here.

        Hashtable labelTable_delts = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_delts.put( new Integer( 1 ), new JLabel("0.001") );
        labelTable_delts.put( new Integer( 100 ), new JLabel("0.10") );
        labelTable_delts.put( new Integer( 200 ), new JLabel("0.20") );
        labelTable_delts.put( new Integer( 300 ), new JLabel("0.30") );
        labelTable_delts.put( new Integer( 400 ), new JLabel("0.40") );
        labelTable_delts.put( new Integer( 500 ), new JLabel("0.50") );
        scroll_8.setLabelTable(labelTable_delts);

        xpos = 4; ypos = 4;             // Scroll
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_8, constraints);
        scroll_8.setMajorTickSpacing(50);
        scroll_8.setPaintLabels(true);
        scroll_8.setPaintTicks(true);
        panel_2.add(scroll_8);

         scroll_8.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_8 = (JSlider) evt.getSource();
                    delt_s = scroll_8.getValue();
                    delt_s = delt_s / 1000;
                    jTextField8.setText(String.valueOf(delt_s));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// AGE AT MATURITY
////////////////////////////////////////////////////////////////////////////////
        xpos = 5; ypos = 7;             // text box
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField9 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField9.setText(String.valueOf(a_m));
        gbl.setConstraints(jTextField9, constraints);
        panel_2.add(jTextField9);

        xpos = 4; ypos = 6;             // Label .
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_9 = new JLabel("Age at Maturity");    // test label while debugging the layout
        label_9.setToolTipText("The age at first possible conception. The age at first possible partuition is one year greater, given a gestation time of approx. one year.");
        gbl.setConstraints(label_9, constraints);
        panel_2.add(label_9);             // see if can add a label here.

        xpos = 4; ypos = 7;             // Scroll
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_9, constraints);
        scroll_9.setMajorTickSpacing(1);
        scroll_9.setPaintLabels(true);
        scroll_9.setPaintTicks(true);
        panel_2.add(scroll_9);

         scroll_9.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_9 = (JSlider) evt.getSource();
                    a_m = scroll_9.getValue();
                    jTextField9.setText(String.valueOf(a_m));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
//////////////////////////////////////////////////////////////////////////////////
//// CV_add
//////////////////////////////////////////////////////////////////////////////////
        xpos = 1; ypos = 7;             // text box
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
        constraints.insets = new Insets(tt,tl,tb,tr);
        final JTextField jTextField15 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField15.setText(String.valueOf(CV_Add));
        gbl.setConstraints(jTextField15, constraints);
        panel_3.add(jTextField15);

        xpos = 0; ypos = 6;             // Label .
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        JLabel label_15 = new JLabel("Additional Uncertainty in Abundance Estimates");    // test label while debugging the layout
        label_15.setToolTipText("Additional unexplained uncertainty about the means of the abundance estimates, the variance estimates of which are known to underestimate the true uncertainty about the estimated abundance.");
        gbl.setConstraints(label_15, constraints);
        panel_3.add(label_15);             // see if can add a label here.

        Hashtable labelTable_CV_add = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_CV_add.put( new Integer( 0 ), new JLabel("0.0") );
        labelTable_CV_add.put( new Integer( 10 ), new JLabel("0.1") );
        labelTable_CV_add.put( new Integer( 20 ), new JLabel("0.2") );
        labelTable_CV_add.put( new Integer( 30 ), new JLabel("0.3") );
        scroll_15.setLabelTable(labelTable_CV_add);

        xpos = 0; ypos = 7;             // Scroll
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_15, constraints);
        scroll_15.setMajorTickSpacing(5);
        scroll_15.setPaintLabels(true);
        scroll_15.setPaintTicks(true);
        panel_3.add(scroll_15);

         scroll_15.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_15 = (JSlider) evt.getSource();
                    CV_Add = scroll_15.getValue();
                    CV_Add = CV_Add / 100;
                    jTextField15.setText(String.valueOf(CV_Add));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// Post-1900 1+ SURVIVAL RATE (if dictated by radio button
////////////////////////////////////////////////////////////////////////////////
        xpos = 7; ypos = 7;             // text box 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
//        constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.insets = new Insets(tt,tl,tb,tr);
        constraints.gridheight = 1;
        final JTextField jTextField77 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField77.setText(String.valueOf(s_post));
        jTextField77.setEnabled(false);
        gbl.setConstraints(jTextField77, constraints);
        panel_3.add(jTextField77);

        xpos = 6; ypos = 6;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        final JLabel label_77 = new JLabel("Age 1+ Survival Post-1900");    // test label while debugging the layout
        label_77.setToolTipText("The annual survival rate of animals aged one year and older (i.e. all ages except calves) post-1900");
        label_77.setEnabled(false);
        gbl.setConstraints(label_77, constraints);
        panel_3.add(label_77);             // see if can add a label here.

        Hashtable labelTable_Sa_post1900 = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_Sa_post1900.put( new Integer( 950 ), new JLabel("0.95") );
        labelTable_Sa_post1900.put( new Integer( 960 ), new JLabel("0.96") );
        labelTable_Sa_post1900.put( new Integer( 970 ), new JLabel("0.97") );
        labelTable_Sa_post1900.put( new Integer( 980 ), new JLabel("0.98") );
        labelTable_Sa_post1900.put( new Integer( 990 ), new JLabel("0.99") );
        scroll_77.setLabelTable(labelTable_Sa_post1900);

        xpos = 6; ypos = 7;             // Scroll 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_77, constraints);
        scroll_77.setMajorTickSpacing(10);
        scroll_77.setPaintLabels(true);
        scroll_77.setPaintTicks(true);
        scroll_77.setEnabled(false);
        panel_3.add(scroll_77);

         scroll_77.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_77 = (JSlider) evt.getSource();
                    s_post = scroll_77.getValue();
                    s_post = s_post / 1000;
                    jTextField77.setText(String.valueOf(s_post));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// Pre-1900 1+ SURVIVAL RATE (if dictated by radio button
////////////////////////////////////////////////////////////////////////////////
        xpos = 7; ypos = 4;             // text box 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
//        constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.insets = new Insets(tt,tl,tb,tr);
        constraints.gridheight = 1;
        final JTextField jTextField777 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField777.setText(String.valueOf(s_post));
        jTextField777.setEnabled(false);
        gbl.setConstraints(jTextField777, constraints);
        panel_3.add(jTextField777);

        xpos = 6; ypos = 3;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        final JLabel label_777 = new JLabel("Age 1+ Survival Pre-1900");    // test label while debugging the layout
        label_777.setToolTipText("The annual survival rate of animals aged one year and older (i.e. all ages except calves) pre-1900");
        label_777.setEnabled(false);
        gbl.setConstraints(label_777, constraints);
        panel_3.add(label_777);             // see if can add a label here.

        Hashtable labelTable_Sa_pre1900 = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_Sa_pre1900.put( new Integer( 950 ), new JLabel("0.95") );
        labelTable_Sa_pre1900.put( new Integer( 960 ), new JLabel("0.96") );
        labelTable_Sa_pre1900.put( new Integer( 970 ), new JLabel("0.97") );
        labelTable_Sa_pre1900.put( new Integer( 980 ), new JLabel("0.98") );
        labelTable_Sa_pre1900.put( new Integer( 990 ), new JLabel("0.99") );
        scroll_777.setLabelTable(labelTable_Sa_pre1900);

        xpos = 6; ypos = 4;             // Scroll pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_777, constraints);
        scroll_777.setMajorTickSpacing(10);
        scroll_777.setPaintLabels(true);
        scroll_777.setPaintTicks(true);
        scroll_777.setEnabled(false);
        panel_3.add(scroll_777);

         scroll_777.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_777 = (JSlider) evt.getSource();
                    s_a = scroll_777.getValue();
                    s_a = s_a / 1000;
                    jTextField777.setText(String.valueOf(s_a));
                    if(Surv_const==true){
                        scroll_777.setValue((int)(s_a*1000));
                        jTextField777.setText(String.valueOf(s_a));
                    }
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

////////////////////////////////////////////////////////////////////////////////
// 1+ SURVIVAL RATE 
////////////////////////////////////////////////////////////////////////////////
        xpos = 5; ypos = 1;             // text box 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
//        constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.insets = new Insets(tt,tl,tb,tr);
        constraints.gridheight = 1;
        final JTextField jTextField7 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField7.setText(String.valueOf(s_a));
        gbl.setConstraints(jTextField7, constraints);
        panel_2.add(jTextField7);

        xpos = 4; ypos = 0;             // Label 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        final JLabel label_7 = new JLabel("Age 1+ Survival");    // test label while debugging the layout
        label_7.setToolTipText("The annual survival rate of animals aged one year and older (i.e. all ages except calves).");
        gbl.setConstraints(label_7, constraints);
        panel_2.add(label_7);             // see if can add a label here.

        Hashtable labelTable_Sa = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_Sa.put( new Integer( 950 ), new JLabel("0.95") );
        labelTable_Sa.put( new Integer( 960 ), new JLabel("0.96") );
        labelTable_Sa.put( new Integer( 970 ), new JLabel("0.97") );
        labelTable_Sa.put( new Integer( 980 ), new JLabel("0.98") );
        labelTable_Sa.put( new Integer( 990 ), new JLabel("0.99") );
        scroll_7.setLabelTable(labelTable_Sa);

        xpos = 4; ypos = 1;             // Scroll 4 pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_7, constraints);
        scroll_7.setMajorTickSpacing(10);
        scroll_7.setPaintLabels(true);
        scroll_7.setPaintTicks(true);
        panel_2.add(scroll_7);

         scroll_7.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_7 = (JSlider) evt.getSource();
                    s_a = scroll_7.getValue();
                    s_a = s_a / 1000;
                    jTextField7.setText(String.valueOf(s_a));
                    if(Surv_const==true){
                        scroll_77.setValue((int)(s_a*1000));
                        jTextField77.setText(String.valueOf(s_a));
                        scroll_777.setValue((int)(s_a*1000));
                        jTextField777.setText(String.valueOf(s_a));
                    }
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });
////////////////////////////////////////////////////////////////////////////////
// Generic Catch (Removals) Multiplier Applied to all fisheries pre-1900
////////////////////////////////////////////////////////////////////////////////
        xpos = 3; ypos = 4;             // text box pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.NORTHWEST;
//        constraints.anchor = GridBagConstraints.NORTH;
//        constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.insets = new Insets(tt,tl,tb,tr);
        constraints.gridheight = 1;
        final JTextField jTextField17 = new JTextField(3);   // needs to be declared as 'FINAL'
        jTextField17.setText(String.valueOf(C_mult));
        gbl.setConstraints(jTextField17, constraints);
        panel_3.add(jTextField17);

        xpos = 2; ypos = 3;             // Label pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(lt,ll,lb,lr);
        constraints.anchor = GridBagConstraints.CENTER;
        constraints.fill = GridBagConstraints.NONE;
        final JLabel label_17 = new JLabel("Landings Multiplier: 1846-1900");    // test label while debugging the layout
        label_17.setToolTipText("Multiplies the estimated landings (all fisheries combined) during the period 1846-1900.");
        gbl.setConstraints(label_17, constraints);
        panel_3.add(label_17);             // see if can add a label here.

        Hashtable labelTable_C_mult = new Hashtable();                 //Create the label table (re-label in units of population growth rate, sliders only take integer values)
        labelTable_C_mult.put( new Integer( 50 ), new JLabel("0.50") );
//        labelTable_C_mult.put( new Integer( 75 ), new JLabel("0.75") );
        labelTable_C_mult.put( new Integer( 100 ), new JLabel("1.0") );
//        labelTable_C_mult.put( new Integer( 125 ), new JLabel("1.25") );
        labelTable_C_mult.put( new Integer( 150 ), new JLabel("1.50") );
//        labelTable_C_mult.put( new Integer( 175 ), new JLabel("1.75") );
        labelTable_C_mult.put( new Integer( 200 ), new JLabel("2.00") );
        scroll_17.setLabelTable(labelTable_C_mult);

        xpos = 2; ypos = 4;             // Scroll pos.
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);         constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(scroll_17, constraints);
        scroll_17.setMajorTickSpacing(25);
        scroll_17.setPaintLabels(true);
        scroll_17.setPaintTicks(true);
        panel_3.add(scroll_17);

         scroll_17.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                JSlider scroll_17 = (JSlider) evt.getSource();
                    C_mult = scroll_17.getValue();
                    C_mult = C_mult / 100;
                    jTextField17.setText(String.valueOf(C_mult));
                    updateModel();                               // updates model trajectory(s)
                    jTextField18.setText(String.valueOf((int)C_init));
                    scroll_18.setValue((int)C_init);
            }
         });

//////////////////////////////////////////////////////////////////////////////////
//// Add radio buttons to allow user to choose between a constant rate of survival or two rates (pre- and post-1900)
//////////////////////////////////////////////////////////////////////////////////
        final JRadioButton rb1 = new JRadioButton("Constant Survival", true);
        final JRadioButton rb2 = new JRadioButton("Two Survival Rates (pre- and post-1900)", false);
//        rb1.addActionListener(this);
//        rb1.addActionListener(this);
//        JRadioButton rb2 = new JRadioButton("Two Survival Rates", false);
        ButtonGroup bg = new ButtonGroup();
        bg.add(rb1);                            // add each radio button to a common button group, so options are mutually exclusive
        bg.add(rb2);

        JPanel radioPanel = new JPanel();
        radioPanel.setLayout(new GridLayout(2, 0));
        radioPanel.add(rb1);
        radioPanel.add(rb2);
        radioPanel.setBorder(BorderFactory.createTitledBorder(
           BorderFactory.createEtchedBorder(), "Age 1+ Survival Rate Through Time"));

        st = 0; sl = 0; sr = 0; sb = 5;        // padding in pixels around the scroll bars//

        xpos = 5; ypos = 0;
        constraints.gridx = xpos;
        constraints.gridy = ypos;
        constraints.insets = new Insets(st,sl,sb,sr);
        constraints.gridheight = 3;
        constraints.gridwidth = 2;
        constraints.fill = GridBagConstraints.NONE;
        constraints.anchor = GridBagConstraints.CENTER;
        gbl.setConstraints(radioPanel, constraints);
        panel_3.add(radioPanel);

        rb1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
               label_7.setEnabled(true);
               jTextField7.setEnabled(true);
               scroll_7.setEnabled(true);
               scroll_7.setValue((int)(s_a*1000));
               jTextField7.setText(String.valueOf(s_a));
               label_7.setText("Age 1+ Survival");
               label_7.setToolTipText("The annual survival rate of animals aged one year and older (i.e. all ages except calves)");

               label_7.setEnabled(true);
               jTextField7.setEnabled(true);
               scroll_7.setEnabled(true);

               scroll_777.setValue((int)(s_a*1000));
               jTextField777.setText(String.valueOf(s_a));
               scroll_77.setValue((int)(s_a*1000));
               jTextField77.setText(String.valueOf(s_a));

               label_777.setEnabled(false);
               jTextField777.setEnabled(false);
               scroll_777.setEnabled(false);

               label_77.setEnabled(false);
               jTextField77.setEnabled(false);
               scroll_77.setEnabled(false);
               label_77.setToolTipText("Select 'Two Survival Rates' above to explore this scenario");
               label_777.setToolTipText("Select 'Two Survival Rates' above to explore this scenario");
               
               Surv_const=true;                             // use this flag to dictate survival rate for projections

               updateModel();                               // updates model trajectory(s)
            }
        });

        rb2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
               label_7.setText("Age 1+ Survival Not Constant");
               label_7.setToolTipText("Change the different survival rates in the Parameter Set 2, or select 'Constant Survival'");
               label_7.setEnabled(false);
               jTextField7.setEnabled(false);
               scroll_7.setEnabled(false);
               jTextField7.setText("NA");

               label_777.setEnabled(true);
               jTextField777.setEnabled(true);
               scroll_777.setEnabled(true);
               label_77.setEnabled(true);
               jTextField77.setEnabled(true);
               scroll_77.setEnabled(true);
               label_77.setToolTipText("Post-1900 annual survival rate of animals aged one year and older (i.e. all ages except calves).");
               label_77.setToolTipText("Pre-1900 annual survival rate of animals aged one year and older (i.e. all ages except calves).");
               Surv_const=false;                            // use this flag to dictate survival rate for projections
               updateModel();                               // updates model trajectory(s)
            }
        });

        JButton reset_button = new JButton("Reset Default Values");
        reset_button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                Default_Values();
                updateModel();
            }
        });

        panel_4.add(reset_button);

        tabbedPane1.addTab("Parameter Set: 1 / 2", null, panel_2, "Life history parameters and uncertainties in historical removals off California and Mexico");
        tabbedPane1.addTab("Parameter Set: 2 / 2", null, panel_3, "Additional uncertainties and possible changes in parameters through time");

//        tabbedPane1.addTab("Additional Controls", null, panel_4, "Reset Default Values");
        tabbedPane1.setTabPlacement(2);                            // put tabs on top of panel (could also try on left side '2')

        content.add(tabbedPane1, BorderLayout.SOUTH);

        //content.add(panel_2, BorderLayout.SOUTH);                // instruct frame to put sliders on bottom panel

    }

//#######################
    public JFreeChart generateGraph(){
        String at;            // string for annotated text, displaying implied aboriginal catch prior to 1846

        NumberAxis abscisse = new NumberAxis("Year");
 //       NumberAxis abscisse = (NumberAxis) domainAxis;
//        abscisse.setAutoRangeIncludesZero(false);
//        abscisse.DEFAULT MINIMUM AXIS VALUE(1846);

        ValueAxis ordonate = new NumberAxis("Abundance");
        ValueAxis ordonate_calf = new NumberAxis("Calf Production");

        CombinedDomainXYPlot combinedDomainXYPlot = new CombinedDomainXYPlot(abscisse);

        XYLineAndShapeRenderer xyLineAndShapeRenderer = new XYLineAndShapeRenderer();

        xyLineAndShapeRenderer.setSeriesPaint(0, Color.RED);    // series 0 = model abundance
        xyLineAndShapeRenderer.setSeriesPaint(1, Color.BLUE);   // series 1 = carrying capacity
        xyLineAndShapeRenderer.setSeriesPaint(2, Color.GREEN);  // series 2 = mature females
        xyLineAndShapeRenderer.setSeriesPaint(3, Color.GREEN);  // series 3 = genetic estimate of mature female bottleneck
        
        Dimension dimension = new Dimension(4,4);
        Shape shape = new Rectangle(dimension);
        xyLineAndShapeRenderer.setSeriesShape(0, shape);
        xyLineAndShapeRenderer.setSeriesShape(1, shape);
        xyLineAndShapeRenderer.setSeriesShape(2, shape);        // ABUNDANCE

        xyLineAndShapeRenderer.setSeriesLinesVisible(0, true);  // model abundance
        xyLineAndShapeRenderer.setSeriesLinesVisible(1, true);  // carrying capacity
        xyLineAndShapeRenderer.setSeriesLinesVisible(2, true);  // mature females
        xyLineAndShapeRenderer.setSeriesLinesVisible(3, true); // abundance estimates

        xyLineAndShapeRenderer.setSeriesShapesVisible(0, false);        // don't plot data points
        xyLineAndShapeRenderer.setSeriesShapesVisible(1, false);
        xyLineAndShapeRenderer.setSeriesShapesVisible(2, false);
        xyLineAndShapeRenderer.setSeriesShapesVisible(3, false);

        xyLineAndShapeRenderer.setDrawSeriesLineAsPath(true);
        
        xyLineAndShapeRenderer.setPlot(combinedDomainXYPlot);

// DEBUGGING: HAVING SOME TROUBLE SETTING X-AXIS OBJECT TO NULL FOR SECOND PLOT, SO TRYING HERE AS WELL WITH FIRST PLOT
        XYPlot xyPlot = new XYPlot(defaultXYDataset,abscisse,ordonate,xyLineAndShapeRenderer);
//        XYPlot xyPlot = new XYPlot(defaultXYDataset,null,ordonate,xyLineAndShapeRenderer);

        xyPlot.getRenderer().setSeriesStroke(0, new BasicStroke(3.0f));     // control line thickness
        xyPlot.getRenderer().setSeriesStroke(1, new BasicStroke(3.0f));     
        xyPlot.getRenderer().setSeriesStroke(2, new BasicStroke(3.0f));
        float dash[] = { 5.0f };
        xyPlot.getRenderer().setSeriesStroke(3, new BasicStroke(3.0f, BasicStroke.CAP_BUTT,
        BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f));

//           NumberAxis xAxis = new NumberAxis();
//           DecimalFormat format = (DecimalFormat)DecimalFormat.getNumberInstance(Locale.ENGLISH);
//           format.applyPattern("#");
//           xAxis.setNumberFormatOverride(format);
//           plot.setDomainAxis(xAxis);

        if(xLowerBound != null && xUpperBound != null){
            xyPlot.getDomainAxis().setLowerBound(xLowerBound);
            xyPlot.getDomainAxis().setUpperBound(xUpperBound);
            xyPlot.getRangeAxis().setLowerBound(N_low);
            xyPlot.getRangeAxis().setUpperBound(N_high);
//            xyPlot.getRangeAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
            xyPlot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());   // seems to eliminate commas in year labels (e.g., 1990 instead of 1,990)
        }

// now see about adding abundance data CIs as an additional dataset with a different renderer
        XYErrorRenderer CI_renderer = new XYErrorRenderer();       // for drawing data CIs
        CI_renderer.setSeriesPaint(0, Color.BLACK);
        CI_renderer.setSeriesShapesVisible(0, true);             // not sure what this will do?
        CI_renderer.setBaseLinesVisible(false);
        CI_renderer.setBaseShapesVisible(false);
        
//        xyPlot.datasetChanged(null)
        xyPlot.setDataset(3, N_data);
        xyPlot.setRenderer(3, CI_renderer);
        xyPlot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

// CATCHES
        XYBarRenderer catchRenderer = new XYBarRenderer();                        // would like to eventually plot catches as bars, but might have to wait
        catchRenderer.setSeriesPaint(0, Color.GRAY);
        catchRenderer.setShadowVisible(false);
        XYPlot catchPlot = new XYPlot(catchComm, abscisse, new NumberAxis("Removals"), catchRenderer);

        catchPlot.getRangeAxis().setLowerBound(0);
        catchPlot.getRangeAxis().setUpperBound(2500);

// ADD SOME ANNOTATED TEXT (FOR IMPLIED ABORIGINAL CATCH PRIOR TO 1846) TO THE CATCH PLOTS
//        XYTextAnnotation annotation = null;
//        Font font = new Font("SansSerif", Font.PLAIN, 14);
//        at = "Aboriginal removals prior to 1846: " + (int)C_init;
//        annotation = new XYTextAnnotation(at, first_yr, 1000.0);
//        annotation.setFont(font);
//        annotation.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
//        catchPlot.addAnnotation(annotation);

// Now see about adding a third plot showing calf production
//        XYPlot Calf_xyPlot = new XYPlot(defaultXYDataset,abscisse,ordonate_calf,xyLineAndShapeRenderer);

        combinedDomainXYPlot.add(xyPlot, 3);                       // add abundance plot (~three times size of catchPlot, given by weight argument)
//        combinedDomainXYPlot.add(Calf_xyPlot, 1);                       // add calf production plot
        combinedDomainXYPlot.add(catchPlot, 1);                       // add abundance plot
        combinedDomainXYPlot.setGap(10.0);                      // space between plots
//        combinedDomainXYPlot.setOrientation(PlotOrientation.VERTICAL);                      // space between plots

        JFreeChart jFreeChart = new JFreeChart(combinedDomainXYPlot);

        return jFreeChart;

    }

//#######################
    public void readData(){
        int ii, jj;                     // index
        double sig_tot;
        //        int first_yr, last_yr;      // First and last year of projection

        N_obs[0][0] = 1967;
        N_obs[0][1] = 1968;
        N_obs[0][2] = 1969;
        N_obs[0][3] = 1970;
        N_obs[0][4] = 1971;
        N_obs[0][5] = 1972;
        N_obs[0][6] = 1973;
        N_obs[0][7] = 1974;
        N_obs[0][8] = 1975;
        N_obs[0][9] = 1976;
        N_obs[0][10] = 1977;
        N_obs[0][11] = 1978;
        N_obs[0][12] = 1979;
        N_obs[0][13] = 1984;
        N_obs[0][14] = 1985;
        N_obs[0][15] = 1987;
        N_obs[0][16] = 1992;
        N_obs[0][17] = 1993;
        N_obs[0][18] = 1995;
        N_obs[0][19] = 1997;
        N_obs[0][20] = 2000;
        N_obs[0][21] = 2001;
        N_obs[0][22] = 2006;

//        N_obs[1][0] = 13776;              // latest published estimates from Rugh et al.
//        N_obs[1][1] = 12869;
//        N_obs[1][2] = 13431;
//        N_obs[1][3] = 11416;
//        N_obs[1][4] = 10406;
//        N_obs[1][5] = 16098;
//        N_obs[1][6] = 15960;
//        N_obs[1][7] = 13812;
//        N_obs[1][8] = 15481;
//        N_obs[1][9] = 16317;
//        N_obs[1][10] = 17996;
//        N_obs[1][11] = 13971;
//        N_obs[1][12] = 17447;
//        N_obs[1][13] = 22862;
//        N_obs[1][14] = 21444;
//        N_obs[1][15] = 22250;
//        N_obs[1][16] = 18844;
//        N_obs[1][17] = 24638;
//        N_obs[1][18] = 24065;
//        N_obs[1][19] = 29758;
//        N_obs[1][20] = 19448;
//        N_obs[1][21] = 18178;
//        N_obs[1][22] = 20110;

        N_obs[1][0] = 13426;              // latest estimates from Laake et al.
        N_obs[1][1] = 14548;
        N_obs[1][2] = 14553;
        N_obs[1][3] = 12771;
        N_obs[1][4] = 11079;
        N_obs[1][5] = 17365;
        N_obs[1][6] = 17375;
        N_obs[1][7] = 15290;
        N_obs[1][8] = 17564;
        N_obs[1][9] = 18377;
        N_obs[1][10] = 19538;
        N_obs[1][11] = 15384;
        N_obs[1][12] = 19763;
        N_obs[1][13] = 23499;
        N_obs[1][14] = 22921;
        N_obs[1][15] = 26916;
        N_obs[1][16] = 15762;
        N_obs[1][17] = 20103;
        N_obs[1][18] = 20944;
        N_obs[1][19] = 21135;
        N_obs[1][20] = 16369;
        N_obs[1][21] = 16033;
        N_obs[1][22] = 19126;

        for(ii = 0; ii < N_n; ii++)           //Set the x coordinates of the time series(s)
           CV_N[0][ii] = N_obs[0][ii];

        CV_N[1][0] = 0.094;              // latest estimates from Laake et al.
        CV_N[1][1] = 0.08;
        CV_N[1][2] = 0.083;
        CV_N[1][3] = 0.081;
        CV_N[1][4] = 0.092;
        CV_N[1][5] = 0.079;
        CV_N[1][6] = 0.082;
        CV_N[1][7] = 0.084;
        CV_N[1][8] = 0.086;
        CV_N[1][9] = 0.08;
        CV_N[1][10] = 0.088;
        CV_N[1][11] = 0.08;
        CV_N[1][12] = 0.083;
        CV_N[1][13] = 0.089;
        CV_N[1][14] = 0.081;
        CV_N[1][15] = 0.058;
        CV_N[1][16] = 0.067;
        CV_N[1][17] = 0.055;
        CV_N[1][18] = 0.061;
        CV_N[1][19] = 0.068;
        CV_N[1][20] = 0.061;
        CV_N[1][21] = 0.069;
        CV_N[1][22] = 0.071;

//        CV_N[1][0] = 0.0785;              // latest published estimates from Rugh et al.
//        CV_N[1][1] = 0.0550;
//        CV_N[1][2] = 0.0564;
//        CV_N[1][3] = 0.0517;
//        CV_N[1][4] = 0.0590;
//        CV_N[1][5] = 0.0518;
//        CV_N[1][6] = 0.0546;
//        CV_N[1][7] = 0.0565;
//        CV_N[1][8] = 0.0600;
//        CV_N[1][9] = 0.0501;
//        CV_N[1][10] = 0.0694;
//        CV_N[1][11] = 0.0539;
//        CV_N[1][12] = 0.0564;
//        CV_N[1][13] = 0.0603;
//        CV_N[1][14] = 0.0522;
//        CV_N[1][15] = 0.0501;
//        CV_N[1][16] = 0.0632;
//        CV_N[1][17] = 0.0599;
//        CV_N[1][18] = 0.0579;
//        CV_N[1][19] = 0.1049;
//        CV_N[1][20] = 0.0967;
//        CV_N[1][21] = 0.0979;
//        CV_N[1][22] = 0.0878;

        N_data.removeSeries(s1);
        //s1 = new YintervalSeries;
        s1.clear();
        for(jj = 1; jj < N_n; jj++) {          //note that loop starts at index = 1 (not zero), which was initialized above
//            CV_N[1][jj] = CV_N[1][jj] + CV_Add;
            sig_tot = Math.sqrt(CV_N[1][jj]*CV_N[1][jj] + CV_Add*CV_Add);
            CV_N[1][jj] = sig_tot;
//            s1.add(N_obs[0][jj], N_obs[1][jj], N_obs[1][jj]*Math.exp(-Math.sqrt(Math.log(1+CV_N[1][jj]*CV_N[1][jj]))), N_obs[1][jj]*Math.exp(Math.sqrt(Math.log(1+CV_N[1][jj]*CV_N[1][jj]))));
            s1.add(N_obs[0][jj], N_obs[1][jj], N_obs[1][jj]*Math.exp(-1.645*sig_tot), N_obs[1][jj]*Math.exp(1.645*sig_tot));
        }

        N_data.addSeries(s1);
        
//        
// catchSeries.add(Catch[0][ii], Catch[1][ii]);
// return catchSeries;
    }

//#######################
    public XYSeries readCatchData() {
//        DefaultXYDataset CatchData = new DefaultXYDataset();
        int ii;

        XYSeries catchSeries = new XYSeries("Removals");

//    private double catchMaleAborig[][];           // Time series of total catches
//    private double catchFemaleAborig[][];           // Time series of total catches

        Catch = new double[2][dim1];
//        catchMaleComm = new double[2][dim1];
//        catchFemaleComm = new double[2][dim1];
        catch_CaliforniaShore = new double[2][dim1];
        catchMale_Other = new double[2][dim1];
        catchFemale_Other = new double[2][dim1];
        catch_Mex = new double[2][dim1];
        C_mex_tmp = new double[2][dim1];            // catches in Mexico after taking into account uncertainty in magnitude of landings
        C_in = new double[2][dim1];                 // catches in the lagoons
        C_out = new double[2][dim1];                // ship based catches in mex. outside lagoons
        catchMale_Mex = new double[2][dim1];        // lagoon catches that were male
        catchFemale_Mex = new double[2][dim1];      // "" that were female
        
        catch_CaliforniaShore[1][0] = 0;
        catch_CaliforniaShore[1][1] = 0;
        catch_CaliforniaShore[1][2] = 0;
        catch_CaliforniaShore[1][3] = 0;
        catch_CaliforniaShore[1][4] = 0;
        catch_CaliforniaShore[1][5] = 0;
        catch_CaliforniaShore[1][6] = 0;
        catch_CaliforniaShore[1][7] = 0;
        catch_CaliforniaShore[1][8] = 13;
        catch_CaliforniaShore[1][9] = 13;
        catch_CaliforniaShore[1][10] = 20;
        catch_CaliforniaShore[1][11] = 26;
        catch_CaliforniaShore[1][12] = 58;
        catch_CaliforniaShore[1][13] = 54;
        catch_CaliforniaShore[1][14] = 93;
        catch_CaliforniaShore[1][15] = 84;
        catch_CaliforniaShore[1][16] = 105;
        catch_CaliforniaShore[1][17] = 107;
        catch_CaliforniaShore[1][18] = 120;
        catch_CaliforniaShore[1][19] = 109;
        catch_CaliforniaShore[1][20] = 111;
        catch_CaliforniaShore[1][21] = 123;
        catch_CaliforniaShore[1][22] = 128;
        catch_CaliforniaShore[1][23] = 119;
        catch_CaliforniaShore[1][24] = 128;
        catch_CaliforniaShore[1][25] = 127;
        catch_CaliforniaShore[1][26] = 125;
        catch_CaliforniaShore[1][27] = 115;
        catch_CaliforniaShore[1][28] = 108;
        catch_CaliforniaShore[1][29] = 95;
        catch_CaliforniaShore[1][30] = 83;
        catch_CaliforniaShore[1][31] = 99;
        catch_CaliforniaShore[1][32] = 92;
        catch_CaliforniaShore[1][33] = 119;
        catch_CaliforniaShore[1][34] = 98;
        catch_CaliforniaShore[1][35] = 94;
        catch_CaliforniaShore[1][36] = 95;
        catch_CaliforniaShore[1][37] = 91;
        catch_CaliforniaShore[1][38] = 94;
        catch_CaliforniaShore[1][39] = 66;
        catch_CaliforniaShore[1][40] = 28;
        catch_CaliforniaShore[1][41] = 30;
        catch_CaliforniaShore[1][42] = 25;
        catch_CaliforniaShore[1][43] = 27;
        catch_CaliforniaShore[1][44] = 20;
        catch_CaliforniaShore[1][45] = 20;
        catch_CaliforniaShore[1][46] = 19;
        catch_CaliforniaShore[1][47] = 19;
        catch_CaliforniaShore[1][48] = 13;
        catch_CaliforniaShore[1][49] = 13;
        catch_CaliforniaShore[1][50] = 6;
        catch_CaliforniaShore[1][51] = 6;
        catch_CaliforniaShore[1][52] = 6;
        catch_CaliforniaShore[1][53] = 2;
        catch_CaliforniaShore[1][54] = 0;
        catch_CaliforniaShore[1][55] = 0;
        catch_CaliforniaShore[1][56] = 0;
        catch_CaliforniaShore[1][57] = 0;
        catch_CaliforniaShore[1][58] = 0;
        catch_CaliforniaShore[1][59] = 0;
        catch_CaliforniaShore[1][60] = 0;
        catch_CaliforniaShore[1][61] = 0;
        catch_CaliforniaShore[1][62] = 0;
        catch_CaliforniaShore[1][63] = 0;
        catch_CaliforniaShore[1][64] = 0;
        catch_CaliforniaShore[1][65] = 0;
        catch_CaliforniaShore[1][66] = 0;
        catch_CaliforniaShore[1][67] = 0;
        catch_CaliforniaShore[1][68] = 0;
        catch_CaliforniaShore[1][69] = 0;
        catch_CaliforniaShore[1][70] = 0;
        catch_CaliforniaShore[1][71] = 0;
        catch_CaliforniaShore[1][72] = 0;
        catch_CaliforniaShore[1][73] = 0;
        catch_CaliforniaShore[1][74] = 0;
        catch_CaliforniaShore[1][75] = 0;
        catch_CaliforniaShore[1][76] = 0;
        catch_CaliforniaShore[1][77] = 0;
        catch_CaliforniaShore[1][78] = 0;
        catch_CaliforniaShore[1][79] = 0;
        catch_CaliforniaShore[1][80] = 0;
        catch_CaliforniaShore[1][81] = 0;
        catch_CaliforniaShore[1][82] = 0;
        catch_CaliforniaShore[1][83] = 0;
        catch_CaliforniaShore[1][84] = 0;
        catch_CaliforniaShore[1][85] = 0;
        catch_CaliforniaShore[1][86] = 0;
        catch_CaliforniaShore[1][87] = 0;
        catch_CaliforniaShore[1][88] = 0;
        catch_CaliforniaShore[1][89] = 0;
        catch_CaliforniaShore[1][90] = 0;
        catch_CaliforniaShore[1][91] = 0;
        catch_CaliforniaShore[1][92] = 0;
        catch_CaliforniaShore[1][93] = 0;
        catch_CaliforniaShore[1][94] = 0;
        catch_CaliforniaShore[1][95] = 0;
        catch_CaliforniaShore[1][96] = 0;
        catch_CaliforniaShore[1][97] = 0;
        catch_CaliforniaShore[1][98] = 0;
        catch_CaliforniaShore[1][99] = 0;
        catch_CaliforniaShore[1][100] = 0;
        catch_CaliforniaShore[1][101] = 0;
        catch_CaliforniaShore[1][102] = 0;
        catch_CaliforniaShore[1][103] = 0;
        catch_CaliforniaShore[1][104] = 0;
        catch_CaliforniaShore[1][105] = 0;
        catch_CaliforniaShore[1][106] = 0;
        catch_CaliforniaShore[1][107] = 0;
        catch_CaliforniaShore[1][108] = 0;
        catch_CaliforniaShore[1][109] = 0;
        catch_CaliforniaShore[1][110] = 0;
        catch_CaliforniaShore[1][111] = 0;
        catch_CaliforniaShore[1][112] = 0;
        catch_CaliforniaShore[1][113] = 0;
        catch_CaliforniaShore[1][114] = 0;
        catch_CaliforniaShore[1][115] = 0;
        catch_CaliforniaShore[1][116] = 0;
        catch_CaliforniaShore[1][117] = 0;
        catch_CaliforniaShore[1][118] = 0;
        catch_CaliforniaShore[1][119] = 0;
        catch_CaliforniaShore[1][120] = 0;
        catch_CaliforniaShore[1][121] = 0;
        catch_CaliforniaShore[1][122] = 0;
        catch_CaliforniaShore[1][123] = 0;
        catch_CaliforniaShore[1][124] = 0;
        catch_CaliforniaShore[1][125] = 0;
        catch_CaliforniaShore[1][126] = 0;
        catch_CaliforniaShore[1][127] = 0;
        catch_CaliforniaShore[1][128] = 0;
        catch_CaliforniaShore[1][129] = 0;
        catch_CaliforniaShore[1][130] = 0;
        catch_CaliforniaShore[1][131] = 0;
        catch_CaliforniaShore[1][132] = 0;
        catch_CaliforniaShore[1][133] = 0;
        catch_CaliforniaShore[1][134] = 0;
        catch_CaliforniaShore[1][135] = 0;
        catch_CaliforniaShore[1][136] = 0;
        catch_CaliforniaShore[1][137] = 0;
        catch_CaliforniaShore[1][138] = 0;
        catch_CaliforniaShore[1][139] = 0;
        catch_CaliforniaShore[1][140] = 0;
        catch_CaliforniaShore[1][141] = 0;
        catch_CaliforniaShore[1][142] = 0;
        catch_CaliforniaShore[1][143] = 0;
        catch_CaliforniaShore[1][144] = 0;
        catch_CaliforniaShore[1][145] = 0;
        catch_CaliforniaShore[1][146] = 0;
        catch_CaliforniaShore[1][147] = 0;
        catch_CaliforniaShore[1][148] = 0;
        catch_CaliforniaShore[1][149] = 0;
        catch_CaliforniaShore[1][150] = 0;
        catch_CaliforniaShore[1][151] = 0;
        catch_CaliforniaShore[1][152] = 0;
        catch_CaliforniaShore[1][153] = 0;
        catch_CaliforniaShore[1][154] = 0;
        catch_CaliforniaShore[1][155] = 0;
        catch_CaliforniaShore[1][156] = 0;
        catch_CaliforniaShore[1][157] = 0;
        catch_CaliforniaShore[1][158] = 0;
        catch_CaliforniaShore[1][159] = 0;
        catch_CaliforniaShore[1][160] = 0;
        catch_CaliforniaShore[1][161] = 0;
        catch_CaliforniaShore[1][162] = 0;

        catch_Mex[1][0] = 28;
        catch_Mex[1][1] = 105;
        catch_Mex[1][2] = 91;
        catch_Mex[1][3] = 0;
        catch_Mex[1][4] = 0;
        catch_Mex[1][5] = 0;
        catch_Mex[1][6] = 55;
        catch_Mex[1][7] = 207;
        catch_Mex[1][8] = 200;
        catch_Mex[1][9] = 141;
        catch_Mex[1][10] = 186;
        catch_Mex[1][11] = 217;
        catch_Mex[1][12] = 527;
        catch_Mex[1][13] = 568;
        catch_Mex[1][14] = 712;
        catch_Mex[1][15] = 606;
        catch_Mex[1][16] = 181;
        catch_Mex[1][17] = 186;
        catch_Mex[1][18] = 283;
        catch_Mex[1][19] = 303;
        catch_Mex[1][20] = 189;
        catch_Mex[1][21] = 252;
        catch_Mex[1][22] = 103;
        catch_Mex[1][23] = 36;
        catch_Mex[1][24] = 37;
        catch_Mex[1][25] = 48;
        catch_Mex[1][26] = 0;
        catch_Mex[1][27] = 8;
        catch_Mex[1][28] = 0;
        catch_Mex[1][29] = 0;
        catch_Mex[1][30] = 0;
        catch_Mex[1][31] = 0;
        catch_Mex[1][32] = 0;
        catch_Mex[1][33] = 0;
        catch_Mex[1][34] = 0;
        catch_Mex[1][35] = 0;
        catch_Mex[1][36] = 0;
        catch_Mex[1][37] = 0;
        catch_Mex[1][38] = 0;
        catch_Mex[1][39] = 0;
        catch_Mex[1][40] = 0;
        catch_Mex[1][41] = 0;
        catch_Mex[1][42] = 0;
        catch_Mex[1][43] = 0;
        catch_Mex[1][44] = 0;
        catch_Mex[1][45] = 0;
        catch_Mex[1][46] = 0;
        catch_Mex[1][47] = 0;
        catch_Mex[1][48] = 0;
        catch_Mex[1][49] = 0;
        catch_Mex[1][50] = 0;
        catch_Mex[1][51] = 0;
        catch_Mex[1][52] = 0;
        catch_Mex[1][53] = 0;
        catch_Mex[1][54] = 0;
        catch_Mex[1][55] = 0;
        catch_Mex[1][56] = 0;
        catch_Mex[1][57] = 0;
        catch_Mex[1][58] = 0;
        catch_Mex[1][59] = 0;
        catch_Mex[1][60] = 0;
        catch_Mex[1][61] = 0;
        catch_Mex[1][62] = 0;
        catch_Mex[1][63] = 0;
        catch_Mex[1][64] = 0;
        catch_Mex[1][65] = 0;
        catch_Mex[1][66] = 0;
        catch_Mex[1][67] = 0;
        catch_Mex[1][68] = 0;
        catch_Mex[1][69] = 0;
        catch_Mex[1][70] = 0;
        catch_Mex[1][71] = 0;
        catch_Mex[1][72] = 0;
        catch_Mex[1][73] = 0;
        catch_Mex[1][74] = 0;
        catch_Mex[1][75] = 0;
        catch_Mex[1][76] = 0;
        catch_Mex[1][77] = 0;
        catch_Mex[1][78] = 0;
        catch_Mex[1][79] = 0;
        catch_Mex[1][80] = 0;
        catch_Mex[1][81] = 0;
        catch_Mex[1][82] = 0;
        catch_Mex[1][83] = 0;
        catch_Mex[1][84] = 0;
        catch_Mex[1][85] = 0;
        catch_Mex[1][86] = 0;
        catch_Mex[1][87] = 0;
        catch_Mex[1][88] = 0;
        catch_Mex[1][89] = 0;
        catch_Mex[1][90] = 0;
        catch_Mex[1][91] = 0;
        catch_Mex[1][92] = 0;
        catch_Mex[1][93] = 0;
        catch_Mex[1][94] = 0;
        catch_Mex[1][95] = 0;
        catch_Mex[1][96] = 0;
        catch_Mex[1][97] = 0;
        catch_Mex[1][98] = 0;
        catch_Mex[1][99] = 0;
        catch_Mex[1][100] = 0;
        catch_Mex[1][101] = 0;
        catch_Mex[1][102] = 0;
        catch_Mex[1][103] = 0;
        catch_Mex[1][104] = 0;
        catch_Mex[1][105] = 0;
        catch_Mex[1][106] = 0;
        catch_Mex[1][107] = 0;
        catch_Mex[1][108] = 0;
        catch_Mex[1][109] = 0;
        catch_Mex[1][110] = 0;
        catch_Mex[1][111] = 0;
        catch_Mex[1][112] = 0;
        catch_Mex[1][113] = 0;
        catch_Mex[1][114] = 0;
        catch_Mex[1][115] = 0;
        catch_Mex[1][116] = 0;
        catch_Mex[1][117] = 0;
        catch_Mex[1][118] = 0;
        catch_Mex[1][119] = 0;
        catch_Mex[1][120] = 0;
        catch_Mex[1][121] = 0;
        catch_Mex[1][122] = 0;
        catch_Mex[1][123] = 0;
        catch_Mex[1][124] = 0;
        catch_Mex[1][125] = 0;
        catch_Mex[1][126] = 0;
        catch_Mex[1][127] = 0;
        catch_Mex[1][128] = 0;
        catch_Mex[1][129] = 0;
        catch_Mex[1][130] = 0;
        catch_Mex[1][131] = 0;
        catch_Mex[1][132] = 0;
        catch_Mex[1][133] = 0;
        catch_Mex[1][134] = 0;
        catch_Mex[1][135] = 0;
        catch_Mex[1][136] = 0;
        catch_Mex[1][137] = 0;
        catch_Mex[1][138] = 0;
        catch_Mex[1][139] = 0;
        catch_Mex[1][140] = 0;
        catch_Mex[1][141] = 0;
        catch_Mex[1][142] = 0;
        catch_Mex[1][143] = 0;
        catch_Mex[1][144] = 0;
        catch_Mex[1][145] = 0;
        catch_Mex[1][146] = 0;
        catch_Mex[1][147] = 0;
        catch_Mex[1][148] = 0;
        catch_Mex[1][149] = 0;
        catch_Mex[1][150] = 0;
        catch_Mex[1][151] = 0;
        catch_Mex[1][152] = 0;
        catch_Mex[1][153] = 0;
        catch_Mex[1][154] = 0;
        catch_Mex[1][155] = 0;
        catch_Mex[1][156] = 0;
        catch_Mex[1][157] = 0;
        catch_Mex[1][158] = 0;
        catch_Mex[1][159] = 0;
        catch_Mex[1][160] = 0;
        catch_Mex[1][161] = 0;
        catch_Mex[1][162] = 0;

        catchMale_Other[1][0] = 98.25;
        catchMale_Other[1][1] = 101.75;
        catchMale_Other[1][2] = 100.75;
        catchMale_Other[1][3] = 96.25;
        catchMale_Other[1][4] = 96.25;
        catchMale_Other[1][5] = 93.5;
        catchMale_Other[1][6] = 96.5;
        catchMale_Other[1][7] = 104;
        catchMale_Other[1][8] = 103.5;
        catchMale_Other[1][9] = 100.5;
        catchMale_Other[1][10] = 103;
        catchMale_Other[1][11] = 104.5;
        catchMale_Other[1][12] = 120.5;
        catchMale_Other[1][13] = 122.5;
        catchMale_Other[1][14] = 130;
        catchMale_Other[1][15] = 86.5;
        catchMale_Other[1][16] = 65;
        catchMale_Other[1][17] = 65;
        catchMale_Other[1][18] = 70;
        catchMale_Other[1][19] = 71;
        catchMale_Other[1][20] = 65;
        catchMale_Other[1][21] = 68.5;
        catchMale_Other[1][22] = 61;
        catchMale_Other[1][23] = 57.5;
        catchMale_Other[1][24] = 57.5;
        catchMale_Other[1][25] = 58;
        catchMale_Other[1][26] = 55.5;
        catchMale_Other[1][27] = 56;
        catchMale_Other[1][28] = 55.5;
        catchMale_Other[1][29] = 55.5;
        catchMale_Other[1][30] = 55;
        catchMale_Other[1][31] = 55;
        catchMale_Other[1][32] = 55;
        catchMale_Other[1][33] = 55;
        catchMale_Other[1][34] = 55;
        catchMale_Other[1][35] = 54;
        catchMale_Other[1][36] = 54;
        catchMale_Other[1][37] = 54;
        catchMale_Other[1][38] = 54;
        catchMale_Other[1][39] = 54;
        catchMale_Other[1][40] = 54;
        catchMale_Other[1][41] = 54;
        catchMale_Other[1][42] = 54;
        catchMale_Other[1][43] = 54;
        catchMale_Other[1][44] = 54;
        catchMale_Other[1][45] = 31;
        catchMale_Other[1][46] = 31;
        catchMale_Other[1][47] = 31;
        catchMale_Other[1][48] = 31;
        catchMale_Other[1][49] = 31;
        catchMale_Other[1][50] = 31;
        catchMale_Other[1][51] = 31;
        catchMale_Other[1][52] = 31;
        catchMale_Other[1][53] = 31;
        catchMale_Other[1][54] = 31;
        catchMale_Other[1][55] = 30.5;
        catchMale_Other[1][56] = 30.5;
        catchMale_Other[1][57] = 30.5;
        catchMale_Other[1][58] = 30.5;
        catchMale_Other[1][59] = 28.5;
        catchMale_Other[1][60] = 28.5;
        catchMale_Other[1][61] = 28.5;
        catchMale_Other[1][62] = 28.5;
        catchMale_Other[1][63] = 28.5;
        catchMale_Other[1][64] = 29;
        catchMale_Other[1][65] = 29;
        catchMale_Other[1][66] = 28;
        catchMale_Other[1][67] = 29;
        catchMale_Other[1][68] = 38;
        catchMale_Other[1][69] = 28;
        catchMale_Other[1][70] = 26;
        catchMale_Other[1][71] = 26;
        catchMale_Other[1][72] = 26;
        catchMale_Other[1][73] = 26;
        catchMale_Other[1][74] = 27;
        catchMale_Other[1][75] = 46;
        catchMale_Other[1][76] = 33;
        catchMale_Other[1][77] = 26;
        catchMale_Other[1][78] = 27;
        catchMale_Other[1][79] = 100;
        catchMale_Other[1][80] = 51;
        catchMale_Other[1][81] = 36;
        catchMale_Other[1][82] = 31;
        catchMale_Other[1][83] = 24;
        catchMale_Other[1][84] = 23;
        catchMale_Other[1][85] = 5;
        catchMale_Other[1][86] = 10;
        catchMale_Other[1][87] = 60;
        catchMale_Other[1][88] = 88;
        catchMale_Other[1][89] = 74;
        catchMale_Other[1][90] = 95;
        catchMale_Other[1][91] = 11;
        catchMale_Other[1][92] = 30;
        catchMale_Other[1][93] = 18;
        catchMale_Other[1][94] = 51;
        catchMale_Other[1][95] = 33;
        catchMale_Other[1][96] = 53;
        catchMale_Other[1][97] = 52;
        catchMale_Other[1][98] = 2;
        catchMale_Other[1][99] = 24;
        catchMale_Other[1][100] = 13;
        catchMale_Other[1][101] = 11;
        catchMale_Other[1][102] = 7;
        catchMale_Other[1][103] = 10;
        catchMale_Other[1][104] = 5;
        catchMale_Other[1][105] = 6;
        catchMale_Other[1][106] = 17;
        catchMale_Other[1][107] = 21;
        catchMale_Other[1][108] = 15;
        catchMale_Other[1][109] = 22;
        catchMale_Other[1][110] = 46;
        catchMale_Other[1][111] = 36;
        catchMale_Other[1][112] = 56;
        catchMale_Other[1][113] = 75;
        catchMale_Other[1][114] = 58;
        catchMale_Other[1][115] = 78;
        catchMale_Other[1][116] = 59;
        catchMale_Other[1][117] = 68;
        catchMale_Other[1][118] = 90;
        catchMale_Other[1][119] = 71;
        catchMale_Other[1][120] = 95;
        catchMale_Other[1][121] = 155;
        catchMale_Other[1][122] = 89;
        catchMale_Other[1][123] = 90;
        catchMale_Other[1][124] = 71;
        catchMale_Other[1][125] = 58;
        catchMale_Other[1][126] = 62;
        catchMale_Other[1][127] = 97;
        catchMale_Other[1][128] = 94;
        catchMale_Other[1][129] = 58;
        catchMale_Other[1][130] = 69;
        catchMale_Other[1][131] = 87;
        catchMale_Other[1][132] = 94;
        catchMale_Other[1][133] = 58;
        catchMale_Other[1][134] = 54;
        catchMale_Other[1][135] = 36;
        catchMale_Other[1][136] = 57;
        catchMale_Other[1][137] = 46;
        catchMale_Other[1][138] = 59;
        catchMale_Other[1][139] = 55;
        catchMale_Other[1][140] = 46;
        catchMale_Other[1][141] = 48;
        catchMale_Other[1][142] = 44;
        catchMale_Other[1][143] = 61;
        catchMale_Other[1][144] = 67;
        catchMale_Other[1][145] = 67;
        catchMale_Other[1][146] = 0;
        catchMale_Other[1][147] = 0;
        catchMale_Other[1][148] = 21;
        catchMale_Other[1][149] = 48;
        catchMale_Other[1][150] = 18;
        catchMale_Other[1][151] = 48;
        catchMale_Other[1][152] = 64;
        catchMale_Other[1][153] = 69;
        catchMale_Other[1][154] = 63;
        catchMale_Other[1][155] = 62;
        catchMale_Other[1][156] = 80;
        catchMale_Other[1][157] = 71;
        catchMale_Other[1][158] = 43;
        catchMale_Other[1][159] = 49;
        catchMale_Other[1][160] = 57;
        catchMale_Other[1][161] = 51;
        catchMale_Other[1][162] = 64;

        catchFemale_Other[1][0] = 98.25;
        catchFemale_Other[1][1] = 101.75;
        catchFemale_Other[1][2] = 100.75;
        catchFemale_Other[1][3] = 96.25;
        catchFemale_Other[1][4] = 96.25;
        catchFemale_Other[1][5] = 93.5;
        catchFemale_Other[1][6] = 96.5;
        catchFemale_Other[1][7] = 104;
        catchFemale_Other[1][8] = 103.5;
        catchFemale_Other[1][9] = 100.5;
        catchFemale_Other[1][10] = 103;
        catchFemale_Other[1][11] = 104.5;
        catchFemale_Other[1][12] = 120.5;
        catchFemale_Other[1][13] = 122.5;
        catchFemale_Other[1][14] = 130;
        catchFemale_Other[1][15] = 86.5;
        catchFemale_Other[1][16] = 65;
        catchFemale_Other[1][17] = 65;
        catchFemale_Other[1][18] = 70;
        catchFemale_Other[1][19] = 71;
        catchFemale_Other[1][20] = 65;
        catchFemale_Other[1][21] = 68.5;
        catchFemale_Other[1][22] = 61;
        catchFemale_Other[1][23] = 57.5;
        catchFemale_Other[1][24] = 57.5;
        catchFemale_Other[1][25] = 58;
        catchFemale_Other[1][26] = 55.5;
        catchFemale_Other[1][27] = 56;
        catchFemale_Other[1][28] = 55.5;
        catchFemale_Other[1][29] = 55.5;
        catchFemale_Other[1][30] = 55;
        catchFemale_Other[1][31] = 55;
        catchFemale_Other[1][32] = 55;
        catchFemale_Other[1][33] = 55;
        catchFemale_Other[1][34] = 55;
        catchFemale_Other[1][35] = 54;
        catchFemale_Other[1][36] = 54;
        catchFemale_Other[1][37] = 54;
        catchFemale_Other[1][38] = 54;
        catchFemale_Other[1][39] = 54;
        catchFemale_Other[1][40] = 54;
        catchFemale_Other[1][41] = 54;
        catchFemale_Other[1][42] = 54;
        catchFemale_Other[1][43] = 54;
        catchFemale_Other[1][44] = 54;
        catchFemale_Other[1][45] = 31;
        catchFemale_Other[1][46] = 31;
        catchFemale_Other[1][47] = 31;
        catchFemale_Other[1][48] = 31;
        catchFemale_Other[1][49] = 31;
        catchFemale_Other[1][50] = 31;
        catchFemale_Other[1][51] = 31;
        catchFemale_Other[1][52] = 31;
        catchFemale_Other[1][53] = 31;
        catchFemale_Other[1][54] = 31;
        catchFemale_Other[1][55] = 30.5;
        catchFemale_Other[1][56] = 30.5;
        catchFemale_Other[1][57] = 30.5;
        catchFemale_Other[1][58] = 30.5;
        catchFemale_Other[1][59] = 28.5;
        catchFemale_Other[1][60] = 28.5;
        catchFemale_Other[1][61] = 28.5;
        catchFemale_Other[1][62] = 28.5;
        catchFemale_Other[1][63] = 28.5;
        catchFemale_Other[1][64] = 29;
        catchFemale_Other[1][65] = 29;
        catchFemale_Other[1][66] = 29;
        catchFemale_Other[1][67] = 29;
        catchFemale_Other[1][68] = 38;
        catchFemale_Other[1][69] = 29;
        catchFemale_Other[1][70] = 26;
        catchFemale_Other[1][71] = 26;
        catchFemale_Other[1][72] = 26;
        catchFemale_Other[1][73] = 26;
        catchFemale_Other[1][74] = 27;
        catchFemale_Other[1][75] = 44;
        catchFemale_Other[1][76] = 29;
        catchFemale_Other[1][77] = 26;
        catchFemale_Other[1][78] = 26;
        catchFemale_Other[1][79] = 86;
        catchFemale_Other[1][80] = 43;
        catchFemale_Other[1][81] = 48;
        catchFemale_Other[1][82] = 33;
        catchFemale_Other[1][83] = 26;
        catchFemale_Other[1][84] = 24;
        catchFemale_Other[1][85] = 5;
        catchFemale_Other[1][86] = 10;
        catchFemale_Other[1][87] = 55;
        catchFemale_Other[1][88] = 78;
        catchFemale_Other[1][89] = 80;
        catchFemale_Other[1][90] = 103;
        catchFemale_Other[1][91] = 13;
        catchFemale_Other[1][92] = 34;
        catchFemale_Other[1][93] = 21;
        catchFemale_Other[1][94] = 74;
        catchFemale_Other[1][95] = 44;
        catchFemale_Other[1][96] = 68;
        catchFemale_Other[1][97] = 67;
        catchFemale_Other[1][98] = 4;
        catchFemale_Other[1][99] = 34;
        catchFemale_Other[1][100] = 17;
        catchFemale_Other[1][101] = 20;
        catchFemale_Other[1][102] = 12;
        catchFemale_Other[1][103] = 16;
        catchFemale_Other[1][104] = 6;
        catchFemale_Other[1][105] = 8;
        catchFemale_Other[1][106] = 27;
        catchFemale_Other[1][107] = 27;
        catchFemale_Other[1][108] = 24;
        catchFemale_Other[1][109] = 37;
        catchFemale_Other[1][110] = 76;
        catchFemale_Other[1][111] = 60;
        catchFemale_Other[1][112] = 92;
        catchFemale_Other[1][113] = 121;
        catchFemale_Other[1][114] = 98;
        catchFemale_Other[1][115] = 130;
        catchFemale_Other[1][116] = 92;
        catchFemale_Other[1][117] = 112;
        catchFemale_Other[1][118] = 129;
        catchFemale_Other[1][119] = 110;
        catchFemale_Other[1][120] = 125;
        catchFemale_Other[1][121] = 219;
        catchFemale_Other[1][122] = 112;
        catchFemale_Other[1][123] = 124;
        catchFemale_Other[1][124] = 80;
        catchFemale_Other[1][125] = 95;
        catchFemale_Other[1][126] = 120;
        catchFemale_Other[1][127] = 81;
        catchFemale_Other[1][128] = 90;
        catchFemale_Other[1][129] = 113;
        catchFemale_Other[1][130] = 96;
        catchFemale_Other[1][131] = 100;
        catchFemale_Other[1][132] = 90;
        catchFemale_Other[1][133] = 125;
        catchFemale_Other[1][134] = 128;
        catchFemale_Other[1][135] = 100;
        catchFemale_Other[1][136] = 111;
        catchFemale_Other[1][137] = 125;
        catchFemale_Other[1][138] = 110;
        catchFemale_Other[1][139] = 115;
        catchFemale_Other[1][140] = 125;
        catchFemale_Other[1][141] = 111;
        catchFemale_Other[1][142] = 107;
        catchFemale_Other[1][143] = 119;
        catchFemale_Other[1][144] = 95;
        catchFemale_Other[1][145] = 102;
        catchFemale_Other[1][146] = 0;
        catchFemale_Other[1][147] = 0;
        catchFemale_Other[1][148] = 23;
        catchFemale_Other[1][149] = 44;
        catchFemale_Other[1][150] = 25;
        catchFemale_Other[1][151] = 31;
        catchFemale_Other[1][152] = 61;
        catchFemale_Other[1][153] = 55;
        catchFemale_Other[1][154] = 52;
        catchFemale_Other[1][155] = 50;
        catchFemale_Other[1][156] = 51;
        catchFemale_Other[1][157] = 57;
        catchFemale_Other[1][158] = 68;
        catchFemale_Other[1][159] = 75;
        catchFemale_Other[1][160] = 77;
        catchFemale_Other[1][161] = 81;
        catchFemale_Other[1][162] = 66;

        catchSeries.clear();

        for(ii = 0; ii <= dim1-1; ii++) {          // perform modifications to historical commercial catches, including CA shore thru 1899
            Catch[0][ii] = first_yr+ii;        // assign years to first column
            catchMaleComm[0][ii] = Catch[0][ii];    // ""
            catchFemaleComm[0][ii] = Catch[0][ii];  // ""
            catchFemale_Other[0][ii] = Catch[0][ii];    // ""
            catchMale_Other[0][ii] = Catch[0][ii];    // ""
            catch_CaliforniaShore[0][ii] = Catch[0][ii];    // ""
            catch_Mex[0][ii] = Catch[0][ii];    // ""

            if (ii<54){                        // apply generic catch multiplier prior to 1900
                C_mex_tmp[1][ii] =  C_mex_tmp[1][ii] * C_mult;
                catch_CaliforniaShore[1][ii] = catch_CaliforniaShore[1][ii] * C_mult;
                catchMale_Other[1][ii] = catchMale_Other[1][ii] * C_mult;
                catchFemale_Other[1][ii] = catchFemale_Other[1][ii] * C_mult;
            }
            
            C_mex_tmp[1][ii] = Math.exp(eps_Mex) * catch_Mex[1][ii];	// Apply catch multiplier to ship based landings in Mex.
            C_mex_tmp[1][ii] = C_mex_tmp[1][ii] * LRF_Mex;              // Loss Rate Factor (LRF_MEX) multipler to calculate removals, given landings
            C_in[1][ii]=C_mex_tmp[1][ii]*2/3;                           // Ship based catches in Mexico which were in lagoons (assuming 2/3 inside)
            C_out[1][ii]=C_mex_tmp[1][ii]/3;                            // Ship based catches in Mexico which were outside lagoons
            C_out[1][ii]=0.5*C_out[1][ii];                              // Now, 'C_out' is the catch by sex (1:1 ratio) outside the lagoons

            catch_CaliforniaShore[1][ii]=catch_CaliforniaShore[1][ii]*Math.exp(eps_CA);
            catch_CaliforniaShore[1][ii]=catch_CaliforniaShore[1][ii]*LRF_CA;

            catchMale_Other[1][ii] = catchMale_Other[1][ii] + 0.5 * catch_CaliforniaShore[1][ii];
            catchFemale_Other[1][ii] = catchFemale_Other[1][ii] + 0.5 * catch_CaliforniaShore[1][ii];

            Catch[1][ii] = C_mex_tmp[1][ii]+catch_CaliforniaShore[1][ii]+catchMale_Other[1][ii]+catchFemale_Other[1][ii];
            catchSeries.add(Catch[0][ii], Catch[1][ii]);

        }

//        XYSeriesCollection collection_tmp = new XYSeriesCollection();
//        collection_tmp.addSeries(catchSeries);
//        return collection_tmp;
        return catchSeries;
    }

//#######################
    public void Calc_Gen_Time(){
        
    }

}
