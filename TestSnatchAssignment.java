package multimapper;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.StringTokenizer;

public class TestSnatchAssignment {

	public static void main(String[] args) throws IOException {
/*		
			Map<String, Double> snatch_readCounts_from_uniqueMappers = new HashMap<String,Double>();
			Map<String, Double> capReadCounts = new HashMap<String, Double>();
			

			Map<String, Double> snatch_readiness_by_gene_unnormalized = new HashMap<String, Double>();
			Map<String, Double> snatch_readiness_by_gene_normalized = new HashMap<String, Double>();
			Map<String, Double> snatch_preference = new HashMap<String, Double>();
			
			Scanner snatch_readCounts_from_uniqueMappers_file = new Scanner(new File(args[0]));
			double sum_snatch_reads = 0.0;
			while (snatch_readCounts_from_uniqueMappers_file.hasNext())
			{
				String s = snatch_readCounts_from_uniqueMappers_file.nextLine();
				StringTokenizer st = new StringTokenizer(s, ",");
				String geneName = st.nextToken();
				Double readCount = Double.parseDouble(st.nextToken());
				snatch_readCounts_from_uniqueMappers.put(geneName, readCount);
				sum_snatch_reads = sum_snatch_reads + readCount;
			}
			
			
			Scanner cap_readCounts_file = new Scanner(new File(args[0]));
			double sum_capreadcounts = 0.0;
			while (cap_readCounts_file.hasNext())
			{
				String s = cap_readCounts_file.nextLine();
				StringTokenizer st = new StringTokenizer(s, ",");
				String geneName = st.nextToken();
				Double readCount = Double.parseDouble(st.nextToken());
				capReadCounts.put(geneName, readCount);
				sum_capreadcounts = sum_capreadcounts + readCount; 
			}
			
			Iterator<String> geneNames = snatch_readCounts_from_uniqueMappers.keySet().iterator();
			double sum_snatch_readiness_unnnormalized = 0.0;
			while (geneNames.hasNext())
			{
				String geneName = geneNames.next();
				Double capReadCounts_thisgene = capReadCounts.get(geneName);
				Double snatch_readCounts_thisgene = snatch_readCounts_from_uniqueMappers.get(geneName);
				double snatch_readiness = snatch_readCounts_thisgene/(snatch_readCounts_thisgene + capReadCounts_thisgene);
				snatch_readiness_by_gene_unnormalized.put(geneName,new Double(snatch_readiness));
				snatch_preference.put(geneName, (snatch_readCounts_thisgene/sum_snatch_reads));
				sum_snatch_readiness_unnnormalized = sum_snatch_readiness_unnnormalized + snatch_readiness;
			}
			
			geneNames = snatch_readCounts_from_uniqueMappers.keySet().iterator();
			while (geneNames.hasNext())
			{
				String geneName = geneNames.next();
				double snatch_readiness_unnormalized = snatch_readiness_by_gene_unnormalized.get(geneName);
				snatch_readiness_by_gene_normalized.put(geneName,(snatch_readiness_unnormalized/sum_snatch_readiness_unnnormalized));
			}
*/
			// Read counts from unique mappers
			double snatch_readcount_t_A = 1500;
			double snatch_readcount_t_B = 2000;
			double snatch_readcount_t_C = 10000;
			
			double cap_read_count_A = 100000;
			double cap_read_count_B = 1100067;
			double cap_read_count_C = 301000;
			
			
			double totalReads_t_A = snatch_readcount_t_A + cap_read_count_A;
			double totalReads_t_B = snatch_readcount_t_B + cap_read_count_B;
			double totalReads_t_C = snatch_readcount_t_C + cap_read_count_C;
		
			double reads_multi_to_A_and_B = 1000;
			
		
			// Initialize snatch preference vector based on unique mapping counts 
			double snatch_preference_t_A = (snatch_readcount_t_A)/(snatch_readcount_t_A + snatch_readcount_t_B + snatch_readcount_t_C);
			double snatch_preference_t_B = (snatch_readcount_t_B)/(snatch_readcount_t_A + snatch_readcount_t_B + snatch_readcount_t_C);
			double snatch_preference_t_C = (snatch_readcount_t_C)/(snatch_readcount_t_A + snatch_readcount_t_B + snatch_readcount_t_C);
			//System.out.println(snatch_preference_t_A + "***"+ snatch_preference_t_B+"***"+snatch_preference_t_C);
			
			double snatch_contribution_t_A = snatch_readcount_t_A/(snatch_readcount_t_A+cap_read_count_A);
			double snatch_contribution_t_B = snatch_readcount_t_B/(snatch_readcount_t_B+cap_read_count_B);
			double snatch_contribution_t_C = snatch_readcount_t_C/(snatch_readcount_t_C+cap_read_count_C);
			
			double snatch_contribution_t_A_normalized = snatch_contribution_t_A/(snatch_contribution_t_A+snatch_contribution_t_B+snatch_contribution_t_C);
			double snatch_contribution_t_B_normalized = snatch_contribution_t_B/(snatch_contribution_t_A+snatch_contribution_t_B+snatch_contribution_t_C);
			double snatch_contribution_t_C_normalized = snatch_contribution_t_C/(snatch_contribution_t_A+snatch_contribution_t_B+snatch_contribution_t_C);
			
			double abundance_t_A = totalReads_t_A/(totalReads_t_A + totalReads_t_B + totalReads_t_C);
			double abundance_t_B = totalReads_t_B/(totalReads_t_A + totalReads_t_B + totalReads_t_C);
			double abundance_t_C = totalReads_t_C/(totalReads_t_A + totalReads_t_B + totalReads_t_C);
			System.out.println("Initial abundances.......");
			System.out.println(abundance_t_A+ "***" + abundance_t_B + "***" + abundance_t_C);
			System.out.println("Initial snatch preferences.......");
			System.out.println(snatch_preference_t_A+ "***" + snatch_preference_t_B + "***" + snatch_preference_t_C);
			System.out.println("Initial snatch contributions.....");
			System.out.println(snatch_contribution_t_A_normalized+ "***" + snatch_contribution_t_B_normalized + "***" + snatch_contribution_t_C_normalized);
			

			double log_l_t = snatch_readcount_t_A*Math.log10(snatch_preference_t_A*snatch_contribution_t_A_normalized*abundance_t_A) + 
						 	 snatch_readcount_t_B*Math.log10(snatch_preference_t_B*snatch_contribution_t_B_normalized*abundance_t_B) + 
						 	 snatch_readcount_t_C*Math.log10(snatch_preference_t_C*snatch_contribution_t_C_normalized*abundance_t_C) ;
			
			double log_l_t_plus_1 = log_l_t + 0.01;
			
			while ((log_l_t_plus_1 - log_l_t) > 0.0001)
			{
				//System.out.println((log_l_t_plus_1 - log_l_t));
				double snatch_readcount_t_plus_1_A = snatch_contribution_t_A_normalized*snatch_preference_t_A*abundance_t_A*snatch_readcount_t_A + reads_multi_to_A_and_B * abundance_t_A * (snatch_preference_t_A/(snatch_preference_t_A + snatch_preference_t_B))*(snatch_contribution_t_A_normalized/(snatch_contribution_t_A_normalized+snatch_contribution_t_B_normalized));
				double snatch_readcount_t_plus_1_B = snatch_contribution_t_B_normalized*snatch_preference_t_B*abundance_t_B*snatch_readcount_t_B + reads_multi_to_A_and_B * abundance_t_B * (snatch_preference_t_B/(snatch_preference_t_A + snatch_preference_t_B))*(snatch_contribution_t_B_normalized/(snatch_contribution_t_A_normalized+snatch_contribution_t_B_normalized));
				double snatch_readcount_t_plus_1_C = snatch_contribution_t_C_normalized*snatch_preference_t_C*abundance_t_C*snatch_readcount_t_C; 

				double snatch_contribution_t_plus_1_A = snatch_readcount_t_plus_1_A/(snatch_readcount_t_plus_1_A+cap_read_count_A);
				double snatch_contribution_t_plus_1_B = snatch_readcount_t_plus_1_B/(snatch_readcount_t_plus_1_B+cap_read_count_B);
				double snatch_contribution_t_plus_1_C = snatch_readcount_t_plus_1_C/(snatch_readcount_t_plus_1_C+cap_read_count_C);
				
				double snatch_contribution_t_plus_1_A_normalized = snatch_contribution_t_plus_1_A/(snatch_contribution_t_plus_1_A+snatch_contribution_t_plus_1_B+snatch_contribution_t_plus_1_C);
				double snatch_contribution_t_plus_1_B_normalized = snatch_contribution_t_plus_1_B/(snatch_contribution_t_plus_1_A+snatch_contribution_t_plus_1_B+snatch_contribution_t_plus_1_C);
				double snatch_contribution_t_plus_1_C_normalized = snatch_contribution_t_plus_1_C/(snatch_contribution_t_plus_1_A+snatch_contribution_t_plus_1_B+snatch_contribution_t_plus_1_C);

				double snatch_preference_t_plus_1_A = (snatch_readcount_t_plus_1_A)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));
				double snatch_preference_t_plus_1_B = (snatch_readcount_t_plus_1_B)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));
				double snatch_preference_t_plus_1_C = (snatch_readcount_t_plus_1_C)/((snatch_readcount_t_plus_1_A) + (snatch_readcount_t_plus_1_B) + (snatch_readcount_t_plus_1_C));
				
				double abundance_t_plus_1_A = (snatch_readcount_t_plus_1_A+cap_read_count_A)/(snatch_readcount_t_plus_1_A+cap_read_count_A+snatch_readcount_t_plus_1_B+cap_read_count_B+snatch_readcount_t_plus_1_C+cap_read_count_C);
				double abundance_t_plus_1_B = (snatch_readcount_t_plus_1_B+cap_read_count_B)/(snatch_readcount_t_plus_1_A+cap_read_count_A+snatch_readcount_t_plus_1_B+cap_read_count_B+snatch_readcount_t_plus_1_C+cap_read_count_C);
				double abundance_t_plus_1_C = (snatch_readcount_t_plus_1_C+cap_read_count_C)/(snatch_readcount_t_plus_1_A+cap_read_count_A+snatch_readcount_t_plus_1_B+cap_read_count_B+snatch_readcount_t_plus_1_C+cap_read_count_C);

				log_l_t = snatch_readcount_t_A*Math.log10(snatch_preference_t_A*snatch_contribution_t_A_normalized*abundance_t_A) + 
					 	  snatch_readcount_t_B*Math.log10(snatch_preference_t_B*snatch_contribution_t_B_normalized*abundance_t_B) + 
					 	  snatch_readcount_t_C*Math.log10(snatch_preference_t_C*snatch_contribution_t_C_normalized*abundance_t_C) ;
				log_l_t_plus_1 = snatch_readcount_t_plus_1_A*Math.log10(snatch_preference_t_plus_1_A*snatch_contribution_t_plus_1_A_normalized*abundance_t_plus_1_A) + 
					 	 		 snatch_readcount_t_plus_1_B*Math.log10(snatch_preference_t_plus_1_B*snatch_contribution_t_plus_1_B_normalized*abundance_t_plus_1_B) + 
					 	 		 snatch_readcount_t_plus_1_C*Math.log10(snatch_preference_t_plus_1_C*snatch_contribution_t_plus_1_C_normalized*abundance_t_plus_1_C) ;
				
				snatch_preference_t_A = snatch_preference_t_plus_1_A;
				snatch_preference_t_B = snatch_preference_t_plus_1_B;
				snatch_preference_t_C = snatch_preference_t_plus_1_C;
								
				snatch_contribution_t_A_normalized = snatch_contribution_t_plus_1_A_normalized;
				snatch_contribution_t_B_normalized = snatch_contribution_t_plus_1_B_normalized;
				snatch_contribution_t_C_normalized = snatch_contribution_t_plus_1_C_normalized;
				
				abundance_t_A = abundance_t_plus_1_A;
				abundance_t_B = abundance_t_plus_1_B;
				abundance_t_C = abundance_t_plus_1_C;
							
				snatch_readcount_t_A = snatch_readcount_t_plus_1_A;
				snatch_readcount_t_B = snatch_readcount_t_plus_1_B;
				snatch_readcount_t_C = snatch_readcount_t_plus_1_C;
				
			}
			System.out.println("Final snatch preferences.......");
			System.out.println(snatch_preference_t_A+"***"+snatch_preference_t_B+ "***"+snatch_preference_t_C);
			System.out.println("Final snatch contributions.......");
			System.out.println(snatch_contribution_t_A_normalized+"***"+snatch_contribution_t_B_normalized+ "***"+snatch_contribution_t_C_normalized);
			System.out.println("Final abundances.......");
			System.out.println(abundance_t_A+"***"+abundance_t_B+ "***"+abundance_t_C);
			System.out.println((snatch_readcount_t_A+cap_read_count_A)+"***"+(snatch_readcount_t_B+cap_read_count_B)+ "***"+(snatch_readcount_t_C+cap_read_count_C));

	}
	

}
