package multimapper;

public class TestSnatchAssignment {

	public static void main(String[] args) {

			double snatch_read_count_t_A = 1500;
			double snatch_read_count_t_B = 2000;
			double snatch_read_count_t_C = 100000;
			
			double cap_read_count_A = 100000;
			double cap_read_count_B = 11000;
			double cap_read_count_C = 301000;
		
			double reads_multi_to_A_and_B = 10;
			
		
			// Initialize theta0 vector based on unique mapping counts 
			double theta_t_A = (snatch_read_count_t_A)/((snatch_read_count_t_A) + (snatch_read_count_t_B) + (snatch_read_count_t_C));
			double theta_t_B = (snatch_read_count_t_B)/((snatch_read_count_t_A) + (snatch_read_count_t_B) + (snatch_read_count_t_C));
			double theta_t_C = (snatch_read_count_t_C)/((snatch_read_count_t_A) + (snatch_read_count_t_B) + (snatch_read_count_t_C));
			//System.out.println(theta_t_A + "***"+ theta_t_B+"***"+theta_t_C);
			
			double snatch_pref_t_A = snatch_read_count_t_A/(snatch_read_count_t_A+cap_read_count_A);
			double snatch_pref_t_B = snatch_read_count_t_B/(snatch_read_count_t_B+cap_read_count_B);
			double snatch_pref_t_C = snatch_read_count_t_C/(snatch_read_count_t_C+cap_read_count_C);
			
			double snatch_pref_t_A_normalized = snatch_pref_t_A/(snatch_pref_t_A+snatch_pref_t_B+snatch_pref_t_C);
			double snatch_pref_t_B_normalized = snatch_pref_t_B/(snatch_pref_t_A+snatch_pref_t_B+snatch_pref_t_C);
			double snatch_pref_t_C_normalized = snatch_pref_t_C/(snatch_pref_t_A+snatch_pref_t_B+snatch_pref_t_C);

			double log_l_t = Math.log10(theta_t_A*snatch_pref_t_A_normalized*snatch_read_count_t_A + theta_t_B*snatch_pref_t_B_normalized*snatch_read_count_t_B + theta_t_C*snatch_pref_t_C_normalized*snatch_read_count_t_C) ;
			double log_l_t_plus_1 = log_l_t + 0.01;
			
			while ((log_l_t_plus_1 - log_l_t) > 0.0001)
			{
				double snatch_readcount_t_plus_1_A = snatch_read_count_t_A + reads_multi_to_A_and_B * (theta_t_A/(theta_t_A + theta_t_B))*(snatch_pref_t_A_normalized/(snatch_pref_t_A_normalized+snatch_pref_t_B_normalized));
				double snatch_readcount_t_plus_1_B = snatch_read_count_t_B + reads_multi_to_A_and_B * (theta_t_B/(theta_t_A + theta_t_B))*(snatch_pref_t_B_normalized/(snatch_pref_t_A_normalized+snatch_pref_t_B_normalized));
				double snatch_readcount_t_plus_1_C = snatch_read_count_t_C;

				double snatch_pref_t_plus_1_A = snatch_readcount_t_plus_1_A/(snatch_readcount_t_plus_1_A+cap_read_count_A);
				double snatch_pref_t_plus_1_B = snatch_readcount_t_plus_1_B/(snatch_readcount_t_plus_1_B+cap_read_count_B);
				double snatch_pref_t_plus_1_C = snatch_readcount_t_plus_1_C/(snatch_readcount_t_plus_1_C+cap_read_count_C);
				
				double snatch_pref_t_plus_1_A_normalized = snatch_pref_t_plus_1_A/(snatch_pref_t_plus_1_A+snatch_pref_t_plus_1_B+snatch_pref_t_plus_1_C);
				double snatch_pref_t_plus_1_B_normalized = snatch_pref_t_plus_1_B/(snatch_pref_t_plus_1_A+snatch_pref_t_plus_1_B+snatch_pref_t_plus_1_C);
				double snatch_pref_t_plus_1_C_normalized = snatch_pref_t_plus_1_C/(snatch_pref_t_plus_1_A+snatch_pref_t_plus_1_B+snatch_pref_t_plus_1_C);

				double theta_t_plus_1_A = (snatch_readcount_t_plus_1_A)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));
				double theta_t_plus_1_B = (snatch_readcount_t_plus_1_B)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));
				double theta_t_plus_1_C = (snatch_readcount_t_plus_1_C)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));

				log_l_t =  Math.log10(theta_t_A*snatch_pref_t_A_normalized*snatch_read_count_t_A + theta_t_B*snatch_pref_t_B_normalized*snatch_read_count_t_B + theta_t_C*snatch_pref_t_C_normalized*snatch_read_count_t_C) ;
				log_l_t_plus_1 = Math.log10(theta_t_plus_1_A*snatch_pref_t_plus_1_A_normalized*snatch_readcount_t_plus_1_A + theta_t_plus_1_B*snatch_pref_t_plus_1_B_normalized*snatch_readcount_t_plus_1_B + theta_t_plus_1_C*snatch_pref_t_plus_1_C_normalized*snatch_readcount_t_plus_1_C) ;
				
				theta_t_A = theta_t_plus_1_A;
				theta_t_B = theta_t_plus_1_B;
				theta_t_C = theta_t_plus_1_C;
				
				snatch_read_count_t_A = snatch_readcount_t_plus_1_A;
				snatch_read_count_t_B = snatch_readcount_t_plus_1_B;
				snatch_read_count_t_C = snatch_readcount_t_plus_1_C;
				
				snatch_pref_t_A_normalized = snatch_pref_t_plus_1_A_normalized;
				snatch_pref_t_B_normalized = snatch_pref_t_plus_1_B_normalized;
				snatch_pref_t_C_normalized = snatch_pref_t_plus_1_C_normalized;
				
			}
			System.out.println(theta_t_A+"***"+theta_t_B+ "***"+theta_t_C);
			System.out.println(snatch_pref_t_A_normalized+"***"+snatch_pref_t_B_normalized+ "***"+snatch_pref_t_C_normalized);
			System.out.println((snatch_read_count_t_A+cap_read_count_A)+"***"+(snatch_read_count_t_B+cap_read_count_B)+ "***"+(snatch_read_count_t_C+cap_read_count_C));

	}
	

}
