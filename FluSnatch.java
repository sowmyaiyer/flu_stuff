package multimapper;

import java.awt.List;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.StringTokenizer;

public class FluSnatch {

	public static void main(String[] args) throws IOException {
		
			Map<String, Gene> genes_t = new LinkedHashMap<String, Gene>();
			
			Scanner readCounts_file = new Scanner(new File(args[0]));
			double sum_snatch_contribution = 0.0;
			double sum_read_counts = 0.0;
			double sum_snatches = 0.0;
			while (readCounts_file.hasNextLine())
			{
				String s = readCounts_file.nextLine();
				String[] st = s.split(",");
				String geneName = st[0];
				double capReadCount = Double.parseDouble(st[1]);
				double snatchReadCount = Double.parseDouble(st[2]); 
				
				Gene gene = new Gene(capReadCount,snatchReadCount);
				genes_t.put(geneName, gene);
				double snatch_contribution_this_gene = capReadCount/(capReadCount+snatchReadCount);
				sum_snatch_contribution = sum_snatch_contribution + snatch_contribution_this_gene;
				sum_read_counts = sum_read_counts + snatchReadCount + capReadCount;
				sum_snatches = sum_snatches + snatchReadCount;
			}


/*
HWI-ST560:128:C30W5ACXX:4:1101:10396:78065	2	CCDC69,TMEM105
HWI-ST560:128:C30W5ACXX:4:1101:10405:17473	2	FAM27E2,HEATR5A
HWI-ST560:128:C30W5ACXX:4:1101:10473:52214	2	HSF2,PHLDB1
HWI-ST560:128:C30W5ACXX:4:1101:10496:74555	2	CCL13,CENPA
*/
			ArrayList<MultiRead> multiReadMap = new ArrayList<MultiRead>();
			Scanner multiReads_file = new Scanner(new File(args[1]));
			while (multiReads_file.hasNextLine())
			{
				MultiRead m = new MultiRead();
				
				String s = multiReads_file.nextLine();
				String[] words = s.split("\\t");
				String readId = words[0];
				int multiplicity = Integer.parseInt(words[1]);
				String[] genenames = words[2].split(",");
				m.setMultiplicity(multiplicity);
				m.setGeneNames(genenames);
				multiReadMap.add(m);
			}
			
		
			
			double log_l_t = 0.0;
			Iterator<String> it = genes_t.keySet().iterator();
			while (it.hasNext()) {
				Gene gene = genes_t.get(it.next());
				gene.init(sum_snatch_contribution, sum_snatches, sum_read_counts);
				log_l_t = log_l_t + gene.getSnatchReadCounts()*(Math.log10(gene.getSnatchPreference()) + Math.log10(gene.getSnatchContribution()) + Math.log10(gene.getAbundance()));
			}
			double log_l_t_plus_1 = log_l_t + 0.01;			

			
			
			while ((log_l_t_plus_1 - log_l_t) > 0.0001)
			{
				Map<String, Gene> genes_t_plus_1 = new LinkedHashMap<String, Gene>();
				
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
