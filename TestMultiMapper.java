package multimapper;

public class TestMultiMapper {

	public static void main(String[] args) {

			double readcounts_t_A = 1000;
			double readcounts_t_B = 100;
			double readcounts_t_C = 100;
		
			double reads_multi_to_A_and_B = 90;
			
			double length_A = 200;
			double length_B = 200;
			double length_C = 200;
			
		
			// Initial theta0 vector
			double theta_t_A = (readcounts_t_A/length_A)/((readcounts_t_A/length_A) + (readcounts_t_B/length_B) + (readcounts_t_C/length_C));
			double theta_t_B = (readcounts_t_B/length_B)/((readcounts_t_A/length_A) + (readcounts_t_B/length_B) + (readcounts_t_C/length_C));
			double theta_t_C = (readcounts_t_C/length_C)/((readcounts_t_A/length_A) + (readcounts_t_B/length_B) + (readcounts_t_C/length_C));
			double log_l_t = Math.log10(theta_t_A*readcounts_t_A + theta_t_B*readcounts_t_B + theta_t_C*readcounts_t_C) ;
			double log_l_t_plus_1 = log_l_t + 1;
			
			//while ((log_l_t_plus_1 - log_l_t) > 0.0001)
			{
				double readcounts_t_plus_1_A = readcounts_t_A + reads_multi_to_A_and_B * (theta_t_A/(theta_t_A + theta_t_B));
				double readcounts_t_plus_1_B = readcounts_t_B + reads_multi_to_A_and_B * (theta_t_B/(theta_t_A + theta_t_B));
				double readcounts_t_plus_1_C = readcounts_t_C;


				double theta_t_plus_1_A = (readcounts_t_plus_1_A/length_A)/((readcounts_t_plus_1_A/length_A) + (readcounts_t_plus_1_B/length_B) + (readcounts_t_plus_1_C/length_C));
				double theta_t_plus_1_B = (readcounts_t_plus_1_B/length_B)/((readcounts_t_plus_1_A/length_A) + (readcounts_t_plus_1_B/length_B) + (readcounts_t_plus_1_C/length_C));
				double theta_t_plus_1_C = (readcounts_t_plus_1_C/length_C)/((readcounts_t_plus_1_A/length_A) + (readcounts_t_plus_1_B/length_B) + (readcounts_t_plus_1_C/length_C));

				log_l_t =  Math.log10(theta_t_A*readcounts_t_A + theta_t_B*readcounts_t_B + theta_t_C*readcounts_t_C) ;
				log_l_t_plus_1 = Math.log10(theta_t_plus_1_A*readcounts_t_plus_1_A + theta_t_plus_1_B*readcounts_t_plus_1_B + theta_t_plus_1_C*readcounts_t_plus_1_C) ;
				
				theta_t_A = theta_t_plus_1_A;
				theta_t_B = theta_t_plus_1_B;
				theta_t_C = theta_t_plus_1_C;
				readcounts_t_A = readcounts_t_plus_1_A;
				readcounts_t_B = readcounts_t_plus_1_B;
				readcounts_t_C = readcounts_t_plus_1_C;
				System.out.println(theta_t_A+"***"+theta_t_B+ "***"+theta_t_C);
				System.out.println(readcounts_t_A+"***"+readcounts_t_B+ "***"+readcounts_t_C);
			}
	}

}
