import java.io.*;

class nwchem_Task extends Object{

  public String theory = new String("md");
  public String operation = new String("energy");

  public String system = new String("system_calc");

  public String integration = new String("leapfrog");
  public double time = 0.0;
  public double timestep = 0.001;
  public int equi = 0;
  public int data = 1000;

  public boolean isotherm = false;
  public boolean isotherm_both = false;
  public double temp = 298.15;
  public double tmprlxw = 0.1;
  public double tmprlxs = 0.1;

  public boolean isobar = false;
  public double press = 1.0125e5;
  public double prsrlx = 0.4;
  public double compr = 4.53e-10;

  public double temp_reas = 298.15;
  public int freq_reas = 0;

  public int interaction;
  public int set = 1;
  public boolean pset2 = false;
  public boolean pset3 = false;
  public int polmax = 14;
  public double poltol = 0.001;
  
  public int gridx = 8;
  public int gridy = 8;
  public int gridz = 8;
  public int order = 4;
  public int nodpme = 0;
  public int fft = 1;
  public double alpha = 0.0;

  public boolean twin = false;
  public double r_short = 0.9;
  public double r_long = 1.2;
  public boolean qmmm = false;
  public double r_qmmm = 0.4;

  public int sd_iter = 100;
  public double sd_init = 0.01;
  public double sd_xmin = 0.0001;
  public double sd_xmax = 0.1;

  public int cg_iter = 0;
  public double cg_init = 0.3;
  public double cg_xmin = 0.0001;
  public double cg_xmax = 0.1;
  public int cg_cycle = 10;

  public double ti_error = 5.0;
  public double ti_drift = 5.0;
  public double ti_factor = 0.75;
  public double ti_delta = 0.075;
  public int ti_steps = 21;
  public int ti_max = 21;
  public int ti_over = 1000;
  public boolean ti_decomp = false;
  public boolean ti_sss = false;
  public String ti_direct = new String("forward");
  public String ti_start = new String("new");

  public boolean prt_top = false;
  public boolean prt_topn = false;
  public boolean prt_topw = false;
  public boolean prt_tops = false;
  public boolean prt_extra = false;
  public boolean prt_energy = false;
  public int prt_step = 10;
  public int prt_stat = 100;
  public int prt_expect = 1;
  
  public int freq_pair = 1;
  public int freq_long = 1;
  public int freq_center = 0;
  public int freq_motion = 0;
  public int freq_rdf = 0;
  public int num_rdf = 1000;
  public double r_rdf = 0.9;
  public int freq_coo = 0;
  public int freq_sco = 0;
  public int freq_vel = 0;
  public int freq_svl = 0;
  public int freq_rst = 0;
  public int freq_prp = 0;
  public int freq_syn = 0;
  public int freq_tim = 0;
  public int freq_mnd = 0;
  public int freq_gib = 0;
  
  public boolean cnv = false;
  public boolean fet = false;
  public boolean acf = false;
  public boolean keep = false;
  
  public String format = new String("ecce");
  
  public boolean fix = false;
  public boolean unfix = false;
  public boolean fix_s = false;
  public boolean fix_w = false;
  public boolean fix_X = false;
  public boolean fix_all = false;
  public boolean fix_none = false;

  public boolean shake_w = true;
  public boolean shake_s = true;
  public int shake_w_iter = 100;
  public int shake_s_iter = 100;
  public double shake_w_toler = 0.001;
  public double shake_s_toler = 0.001;

  public int npx = 0;
  public int npy = 0;
  public int npz = 0;
  public int nbx = 0;
  public int nby = 0;
  public int nbz = 0;
  public int mad = 6;
  
  public int test = 0;
  public int debug = 0;

  public int mwm = 0;
  public int msa = 0;

  public boolean load_none = false;
  public boolean load_pair = false;
  public boolean load_size = false;
  public boolean load_reset = false;
  public double load_factor = 0.75;
  public int load_num_pair = 5;

  public boolean distar = false;
  public boolean draver = false;
  public double drsscl = 0.0;
  public int nfdrss = 0;
  public int Limit = 0;
  public int Heap = 1;
  public int Stack = 24;
  public int Global = 12;
  public boolean Verify = false;

  nwchem_Task(){
  }
}

