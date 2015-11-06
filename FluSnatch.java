package multimapper;

import java.util.List;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

public class FluSnatch {

	public static void main(String[] args) throws IOException {
		
			Map<String, Gene> genes_t = new LinkedHashMap<String, Gene>();
			List<String> all_geneNames = new ArrayList<String>();
			
			Scanner readCounts_file = new Scanner(new File(args[0]));
			double sum_snatch_contribution = 0.0;
			double sum_read_counts = 0.0;
			double sum_snatches = 0.0;
			while (readCounts_file.hasNextLine())
			{
				String s = readCounts_file.nextLine();
				String[] st = s.split("\\t");
				String geneName = st[0];
				all_geneNames.add(geneName);
				double capReadCount = Double.parseDouble(st[1]);
				double snatchReadCount = Double.parseDouble(st[2]); 
				
				Gene gene = new Gene(capReadCount,snatchReadCount);
				genes_t.put(geneName, gene);
				double snatch_contribution_this_gene = snatchReadCount/(capReadCount+snatchReadCount);
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
				List<String> genenames = Arrays.asList(words[2].split(","));
				m.setMultiplicity(multiplicity);
				m.setGeneNames(genenames);
				multiReadMap.add(m);
			}
			
		
			
			double log_l_t = 0.0;
			Iterator<String> it = genes_t.keySet().iterator();
			while (it.hasNext()) {
				Gene gene = genes_t.get(it.next());
				gene.init(sum_snatch_contribution, sum_snatches, sum_read_counts);
				log_l_t = log_l_t + gene.getSnatchReadCounts()*(Math.log10(gene.getSnatchPreference()) + Math.log10(gene.getSnatchContributionNormalized()) + Math.log10(gene.getAbundance()));
			}
			double log_l_t_plus_1 = log_l_t + 0.01;			

			
			Map<String, Gene> genes_t_plus_1 = new HashMap<String, Gene>();
			while ((log_l_t_plus_1 - log_l_t) > 0.0001)
			{
				
				Iterator<String> it_gene = all_geneNames.iterator();
				// Initialize genes_t_plus_1 with values from genes_t
				while (it_gene.hasNext())
				{
					String genename = it_gene.next();
					Gene gene_t = genes_t.get(genename);
					Gene gene_t_plus_1 = new Gene();
					gene_t_plus_1.init(gene_t);
					//double new_snatch_read_counts = gene_t.getSnatchContributionNormalized()*gene_t.getSnatchPreference()*gene_t.getAbundance()*gene_t.getSnatchReadCounts(); //snatch read counts from unique_mappers. to be added to multimappers later
					double new_snatch_read_counts = gene_t.getSnatchReadCounts();
					gene_t_plus_1.setSnatchReadCounts(new_snatch_read_counts);
					genes_t_plus_1.put(genename, gene_t_plus_1);
				}
				
				//Alter genes_t_plus_1 based on multi reads
				Iterator<MultiRead> it_multi = multiReadMap.iterator();
				while (it_multi.hasNext())
				{
					MultiRead m = it_multi.next();
					List<String> genes = m.getGeneNames();
					Iterator<String> genenames = genes.iterator();
					int multiplicity = m.getMultiplicity();
					double sum_snatch_preference_multi = 0.0;
					double sum_snatch_contribution_multi = 0.0;
					double sum_abundance_multi = 0.0;
					while (genenames.hasNext())
					{
						String thisgenename = genenames.next();
						Gene gene_t = genes_t.get(thisgenename);
						double snatch_preference_this_gene = gene_t.getSnatchPreference();
						sum_snatch_preference_multi = sum_snatch_preference_multi + snatch_preference_this_gene;
					
						double snatch_contribution_this_gene = gene_t.getSnatchContributionNormalized();
						sum_snatch_contribution_multi = sum_snatch_contribution_multi + snatch_contribution_this_gene;
						
						double abundance_this_gene = gene_t.getAbundance();
						sum_abundance_multi = sum_abundance_multi + abundance_this_gene;
						
					}
					
					genenames = genes.iterator();
					while (genenames.hasNext())
					{
						String thisgenename = genenames.next();
						Gene gene_t = genes_t.get(thisgenename);
						Gene gene_t_plus_1 = genes_t_plus_1.get(thisgenename);
						double multi_read_allocation = multiplicity * (gene_t.getAbundance()/sum_abundance_multi) * (gene_t.getSnatchPreference()/sum_snatch_preference_multi) *(gene_t.getSnatchContributionNormalized()/sum_snatch_contribution_multi);
						gene_t_plus_1.setSnatchReadCounts(gene_t_plus_1.getSnatchReadCounts() + multi_read_allocation);
					}
					
				}
				Iterator<String> it_genes = all_geneNames.iterator();
				it_genes = all_geneNames.iterator();
				double sum_snatch_read_counts_t_plus_1 = 0.0;
				double sum_snatch_contributions_t_plus_1 = 0.0;
				double total_snatch_plus_cap_reads_t_plus_1 = 0.0;
				while (it_genes.hasNext())
				{
					String genename = it_genes.next();
					Gene gene_t_plus_1 = genes_t_plus_1.get(genename);
					sum_snatch_read_counts_t_plus_1 = sum_snatch_read_counts_t_plus_1 + gene_t_plus_1.getSnatchReadCounts();
					sum_snatch_contributions_t_plus_1 = sum_snatch_contributions_t_plus_1 + gene_t_plus_1.getSnatchContributionUnNormalized();
					total_snatch_plus_cap_reads_t_plus_1 = total_snatch_plus_cap_reads_t_plus_1 + gene_t_plus_1.getSnatchReadCounts() + gene_t_plus_1.getCapReadCounts();
				}
				
				it_genes = all_geneNames.iterator();
				while (it_genes.hasNext())
				{
					String genename = it_genes.next();
					Gene gene_t_plus_1 = genes_t_plus_1.get(genename);
					gene_t_plus_1.init(sum_snatch_contributions_t_plus_1, sum_snatch_read_counts_t_plus_1, total_snatch_plus_cap_reads_t_plus_1);
					genes_t_plus_1.put(genename, gene_t_plus_1);
				}
				
				
				
				log_l_t = 0.0;
				it = all_geneNames.iterator();
				while (it.hasNext()) {
					String genename = it.next();
					Gene gene = genes_t.get(genename);
					log_l_t = log_l_t + gene.getSnatchReadCounts()*(Math.log10(gene.getSnatchPreference()) + Math.log10(gene.getSnatchContributionNormalized()) + Math.log10(gene.getAbundance()));
				}
				
				log_l_t_plus_1 = 0.0;
				it = all_geneNames.iterator();
				while (it.hasNext()) {
					Gene gene = genes_t_plus_1.get(it.next());
					log_l_t_plus_1 = log_l_t_plus_1 + gene.getSnatchReadCounts()*(Math.log10(gene.getSnatchPreference()) + Math.log10(gene.getSnatchContributionNormalized()) + Math.log10(gene.getAbundance()));
				}
				
				genes_t = new HashMap<String, Gene>();
				Iterator<String> it_gene2 = all_geneNames.iterator();
				while (it_gene2.hasNext())
				{
					String genename = it_gene2.next();
					Gene newgene = new Gene();
					Gene gene_t_plus_1 = genes_t_plus_1.get(genename);
					newgene.setCapReadCounts(gene_t_plus_1.getCapReadCounts());
					newgene.setSnatchReadCounts(gene_t_plus_1.getSnatchReadCounts());
					newgene.init(sum_snatch_contributions_t_plus_1, sum_snatch_read_counts_t_plus_1, total_snatch_plus_cap_reads_t_plus_1);
					genes_t.put(genename, newgene);
				}
				System.out.println(log_l_t + "***" + log_l_t_plus_1);
			}
			
			
			System.out.println(genes_t.get("A").getSnatchReadCounts()+"***"+genes_t.get("B").getSnatchReadCounts()+ "***"+genes_t.get("C").getSnatchReadCounts());

			
	}
	

}
