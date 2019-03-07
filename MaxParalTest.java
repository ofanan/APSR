/* Given n, k, budget and target_epsilon, calculate and print possible values of s, d; 
and return values for s, d which maximize parallism, while guaranteeing SLA and satisfying the budget constraint s*d < budget
*/
import java.lang.*;
import java.util.stream.*;
import java.util.Arrays; 

class MaxParal {
  double[][] P; 
  int MAX_NUM_OF_SCHEDS = 4;
     
  // C'tor. the main effort of construction is allocating and filling the table P, defined as:
  // P(a,b) = \sum_{1<=i<=b, 0<=X_i<=a, \sum X_i<=b} \prod (i^{X_i})
  // Input: maximum number of schedulers to be considered: should be at least 1.
  MaxParal (int max_num_of_scheds) {
    if (max_num_of_scheds < 2) {
        System.out.println ("MaxParal: Wrong input. Maximal number of schedulers should be at least 2.");
        System.exit(0);
    }
    this.MAX_NUM_OF_SCHEDS = max_num_of_scheds;
    this.P = new double [MAX_NUM_OF_SCHEDS+1][MAX_NUM_OF_SCHEDS+1];
    
    int cur_a, cur_b; //loops iterators
    for (cur_b = 0; cur_b <= MAX_NUM_OF_SCHEDS; cur_b++)
      P[0][cur_b] =  1;   
    for (cur_a = 0; cur_a <= MAX_NUM_OF_SCHEDS; cur_a++)
      P[cur_a][1] =  1;   
    for (cur_b=2; cur_b<=this.MAX_NUM_OF_SCHEDS; cur_b++) {
        for (cur_a=1; cur_a<=this.MAX_NUM_OF_SCHEDS; cur_a++) {
            for (int j=0; j<=cur_a; j++) {
                this.P[cur_a][cur_b] += Math.pow (cur_b, j) * P[cur_a-j][cur_b-1];
            }
        }
    }
  }
  
  // For debug only. prints the matrix P. 
  void print_P (int dummy) {
    for (int i = 0; i <= this.MAX_NUM_OF_SCHEDS; i++) {
      for (int j = 0; j <= this.MAX_NUM_OF_SCHEDS; j++) {
        System.out.print(this.P[i][j] + " ");
      }
      System.out.println();
    }  
  }

	void Debug (int dummy) {
		int line = 20;
		System.out.print("[");
		for (int j = 0; j <= this.MAX_NUM_OF_SCHEDS; j++) {
			System.out.print(this.P[line][j] + ", ");
		}
		System.out.print("]\n");
	}
  // Calculate E[Hs | Fs=f] using the formula:
  // E[Hs | Fs=f] = k^{f-1} * \sum {h=1, 2, ... f} [ h*(k-1)!/*(k-h)! P(f-h, h) ]
  double calc_E_H_s_cond_f (int k, int f) {
    if (f <0 ){
        System.out.println("Error: Cannot have negative number of potentially happy agents");
        System.exit (0);           
    }
    else if (f<2) //If the # of pot-happy agents is 0 or 1, so is the expected # of happy agents
        return f;
    else if (f > k) {
        System.out.println("Error: Cannot have more than k potentially happy agents");
        System.exit (0);
    }
    
    // Now we know that f>=2
    double E_H_s_cond_f = 0;
    double k_fact_over_k_min_h_fact;
    for (int h=1; h <=f; h++) {
        k_fact_over_k_min_h_fact = 1;
        for (int j = 1; j < h; j++)
            k_fact_over_k_min_h_fact *= (k-j);
        E_H_s_cond_f += (h * k_fact_over_k_min_h_fact * this.P[f-h][h]);
    }
    return E_H_s_cond_f * Math.pow (k, 1-f);
  }

  // Accessory function for decreasing the complexity of calculating the binomial values.
  // Code taken from: https://introcs.cs.princeton.edu/java/36inheritance/ProbStat.java.html
  public static double logFact(int n) {
    double ans = 0.0;
    for (int i = 1; i <= n; i++)
        ans += Math.log(i);
    return ans;
  }

  // Return the probability of getting exactly k heads when throwing N biased p coins
  // Code taken from: https://introcs.cs.princeton.edu/java/36inheritance/ProbStat.java.html
  public static double binomial(int N, double p, int k) {
    if (p==1) {
        if (k==N)
            return 1;
        else 
            return 0;              
    }
    return Math.exp(logFact(N) - logFact(k) - logFact(N-k) + k*Math.log(p) + (N-k)*Math.log(1-p));
  }
  
  // Calculate epsilon for a sys with given n,k,s,d		// E[Hs] = \su
  double calc_epsilon (int n, int k, int s, int d) {
    if (k > n || s > n) {
        System.out.println ("Wrong inputs to calc_epsilon: n = " + n + " k = " + k + " s = " + s + " d = " + d);
        System.exit (0);
    }
    double sigma = 1 - Math.pow ( (double)(n-k)/n, d);
    double E_H_s = 0;
    for (int f=1; f<=s; f++) {
        E_H_s += this.calc_E_H_s_cond_f (k, f) * binomial (s, sigma, f);
    }
    return 1 - E_H_s / s;
  }
  
  // // Used for debugging. Prints the pre-computed LUT P. 
  // void PrintP (int dummy) {
    // System.out.println(Arrays.deepToString(this.P));       
  // }
  
  /* Class's main function. Does the following:
   1. Finds and prints all the pairs of (s,d) which satisfies the SLA and budget constraints.
   2. Prints and returns the pair of (s,d) with the maximal parallelism (maximal s), which still satisfies the SLA and budget constraints.
  */
  int[] main (int n, int k, int Smax, int budget, double target_epsilon) {
    System.out.println ("Combinatorial analysis: n = " +n+ " k = " +k+ " budget = " +budget+ " max allowed epsilon = " +target_epsilon);
    System.out.println ("**************************************************************************************");
    
    int d = 0; // default value (indicating a failure in finding values which satisfy the SLA & budget constraints)
    int s = 0;
    double calculated_epsilon;
            int[] rtrn_s_d = new int[2];
    
    for (s=1; s <= Smax; s++) {
      d = budget / s;
      calculated_epsilon = this.calc_epsilon (n, k, s, d);
      if (calculated_epsilon > target_epsilon)
        break;
      System.out.println ("s = " +s+ ", d = " +d+ ", epsilon = " +calculated_epsilon);
        rtrn_s_d[0] = s;
        rtrn_s_d[1] = d;
    }
    System.out.println ("Maximum possible parallelism for this system: s = " +rtrn_s_d[0]+ ", d = " +rtrn_s_d[1]);
    return rtrn_s_d; 
  }
  
} 

// Class for checking the class MaxParal.
public class MaxParalTest {
    
  /*Main function, does the following:
  1. Fix parameters (n, k, max_num_of_scheds, budget, target_epsilon required by the SLA
  2. Generate a MaxParal object (named max_paral) with these parameters.
  3. Query max_paral for finding the highest value of s which still satisfies the SLA and budget requirements.
  */
  public static void main(String[] args) {
      final int n = 1000; //total number of hosts
      int k = 200;  //num of hosts with enough available resources
      int max_num_of_scheds = k; //maximal number of schedulers to check
      int budget = n; //Total number of probes of hosts. Every configuration (s, d) should satisfy s*d <= budget
      double target_epsilon = 0.1; //acceptable decline rate, defined by the SLA
      
      MaxParal max_paral = new MaxParal (max_num_of_scheds);
      max_paral.main (n, k, max_num_of_scheds, budget, target_epsilon);
      
  }
}
