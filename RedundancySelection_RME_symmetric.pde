// Agent based simulation of the selection-mutation dynamics on two dimensional symmetric redundant fitness landscape
//
// Individuals have two traits p and q, the selective quantity is r= p^2 + q^2
// Fitness landscape is symmetric, therefore fitness is (negatively) proportional to the distance from the center r
// The phenotype space P(p,q) is the sphere of radius 1 centered in (p=0;q=0)
// P = { (p,q) | p^2 + q^2 <= 1 }


// Parameters

int RandomNumberSeed=3;                // set random number seed
float mu=1;                            // mutation event rate
float gamma=1;                         // competition event rate
float tot_rate=mu+gamma;               // total events rate
float Mutation=0.02;                   // width mutation step
int L=100000;                          // population size
float w=0.001;                         // selection strength
int bin=50;                            // number of bins population histogram
int n_iterations=500;                   // number of iterations within a single time step
int time_idx=0;                        // initialise time
int time_max=200;                      // number of time steps
int pre_itr_max=1000;                  // number of transient iterations that are not accounted for time averages

// Auxiliary arrays and variables

int[] pop_hist=new int[bin];                      // not normalised histogram of fitness-related phenotype trait p
float[] fitness_hist=new float[bin];             // normalised histogram of fitness-related phenotype trait p 
float[] avg_r=new float[bin];                    // temporal vector of average phenotype trait p
Table table,mean_r,hist_r;
String lbl_r,lbl_c;
float avg;

// Define directory and name of the data files

String address="";                                // here put the address of the output data files (optional)
String filename_avgr= address + day() + "-" + month() + "-" + year() + "-" + hour() + "-" + minute() + "_avgr.csv";   // temporal trajectory of average trait p
String filename_hist= address + day() + "-" + month() + "-" + year() + "-" + hour() + "-" + minute() + "_hist.csv";   // temporal average of trait p distribution  
PrintWriter parameters;                        // txt file printing parameters used

// Define new class Player

class Player{
  float p;          // trait 1
  float q;        // trait 2
  float r;        // selective quantity 
  Player (float p_, float q_){
  p=p_;
  q=q_;
  r=pow(p,2)+pow(q,2);
}

void InheritFromMut(Player parent){
   float ptemp=parent.p;
   float qtemp=parent.q;
   float rtemp=parent.r;
   p=parent.p+random(-Mutation,+Mutation);      // mutate trait 1 of random quantity
   q=parent.q+random(-Mutation,+Mutation);      // mutate trait 2 of random quantity
   r=sqrt(pow(p,2)+pow(q,2));                   // update selective quantity
   if((pow(p,2)+pow(q,2))>1) {p=ptemp; q=qtemp;  r=rtemp;}  // if mutation leads outside boundary then reject    
 }

// Define Asexual reproduction without mutation
 
void InheritFrom(Player parent){
   p=parent.p;
   q=parent.q;
   r=parent.r;
   }
}

Player[] player;          // Define population of players

// Initialisation of the system

void setup(){
  int i,l,index;
  float p0,q0;
  randomSeed(RandomNumberSeed);
  for(l=0;l<bin; l++){pop_hist[l]=0; fitness_hist[l]=float(0);}    // Initialise histograms
  
  
  parameters = createWriter(address + day() + "-" + month() + "-" + year() + "-" + hour() + "-" + minute() + "_parameters.txt");
  table = new Table();
  hist_r = new Table();
  for(i=0;i<bin;i++){
  String lbl=str(i);
  table.addColumn(lbl, Table.FLOAT);
  hist_r.addColumn(lbl, Table.FLOAT);
  }
  mean_r = new Table();
  mean_r.addColumn("<r>", Table.FLOAT);

// Initialise population with randomly uniform distributed phenotypes in the Adaptive Space and update histogram

  player =new Player[L];
  for (i=0;i<L;i++){
    float r0=sqrt(random(1.0));
    float theta0=2*PI*random(1.0);
    player[i]= new Player (p0=r0*cos(theta0), q0=r0*sin(theta0));
    player[i].r=r0;
    index=get_index_r(player[i].r);
    pop_hist[index]++;
  }
}

// Runs the algorithm

void draw(){
  
  while(time_idx<time_max){            // the algorithm is repeated time_max times
    for (int pre_itr=0;pre_itr<pre_itr_max;pre_itr++) dynamics();    // the first pre_itr_max times are considered of transient
    for (int itr=0;itr<n_iterations;itr++) dynamics();                 // for each timestep, the dynamics is simulated n_iterations times
    
    // after each timestep, the population distribution and the average statistics are recorded
    update_histogram();
    update_statistics();
    time_idx++;
    println(time_idx);
  }
  
  // Finally, the averages over time are computed and the data are saved
  
  compute_temporal_average();
  
  // Save tables data and parameters.txt file
  
  saveTable(mean_r,filename_avgr);
  saveTable(hist_r,filename_hist);
  parameters.println("RandomNumberSeed = " + str(RandomNumberSeed));
  parameters.println("mu = " + str(mu));
  parameters.println("gamma = " + str(gamma));
  parameters.println("Mutation = " + str(Mutation));
  parameters.println("w = " + str(w));
  parameters.println("bin = " + str(bin));
  parameters.println("n_iterations = " + str(n_iterations));
  parameters.println("time_max = " + str(L));
  parameters.println("pre_itr_max = " + str(pre_itr_max));
  parameters.close();
  exit(); 
}

// Define dynamics

void dynamics(){
  int i,i1,index;
  float r1;    
  for (int repeat=1; repeat<=L; repeat++){  // events are repeated L times
    r1=random(1.0);         
    i=int(random(float(L)));                 // a focal player i is sampled
    index=get_index_r(player[i].r);      
    if(r1<=mu/tot_rate){                    // i mutates with probability mu/tot_rate
      pop_hist[index]--;
      player[i].InheritFromMut(player[i]);
      index=get_index_r(player[i].r);
      pop_hist[index]++;
    }         
    else{                                  // otherwise i competes with another random player i1
      i1=int(random(float(L)));
      while(i1==i) {i1=int(random(float(L)));}
      
      // Player i1 wins the competition and reproduces with probability proportional to the difference
      // between their fitness and i's
      
      if(random(1.0)<=(1+w*(player[i1].r-player[i].r))/2){      // probability i1 wins
        pop_hist[index]--;
        player[i].InheritFrom(player[i1]);      // traits player i become those of player i1
        index=get_index_r(player[i1].r);
        pop_hist[index]++;
      }
      
      // Otherwise player i wins the competition and i1 becomes i
      
      else{
        pop_hist[index]++;
        index=get_index_r(player[i1].r);
        pop_hist[index]--;
        player[i1].InheritFrom(player[i]);
      }
    }          
  }
}

// updates histogram at the end of each timestep

void update_histogram(){    
  int i;
  for(i=0;i<bin;i++){
    fitness_hist[i]=(float(pop_hist[i])/float(L))*bin;
  }
  TableRow row = table.addRow();
  for(i=0;i<bin;i++){
    String lbl=str(i);
    row.setFloat(lbl, fitness_hist[i]);
  }
}    
        
// Updates average trait at the end of each timestep

void update_statistics() {
  int i;
  float average_r=0;
  for(i=0; i<L; i++) {average_r+= (player[i].r)/L;}
  TableRow row = mean_r.addRow();
  row.setFloat("<r>", average_r);
}

// Auxiliary function getting index of histogram    
    
int get_index_r(float x){
  int idx;  
  idx=int(x*bin);    
  if(idx==bin) idx--;
  return idx;
}

// Finally, computes average over time of histogram

void compute_temporal_average() {
 int i,j;
 int lines=table.getRowCount();
  for(j=0; j<bin; j++){
    avg=float(0);
    for(i=0;i<lines;i++){
      avg+=table.getFloat(i,str(j));
    }
    avg_r[j]=avg/float(lines);
  }  
  TableRow row = hist_r.addRow();
  for(i=0;i<bin;i++){
    row.setFloat(str(i), avg_r[i]);
  }        
}
