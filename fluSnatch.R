readCounts_file <- commandArgs(TRUE)[1]
multimapper_file <- commandArgs(TRUE)[2]

df_readCounts <- as.data.frame(read.table(readCounts_file, header=TRUE))

#genes_t_df <- data.frame(row.names=df_readCounts$gene,snatchReadCounts=df_readCounts$snatchReadCounts,capReadCounts=df_readCounts$capReadCounts)
genes_t_df <- data.frame(row.names=df_readCounts$gene,snatchReadCounts=df_readCounts$snatchReadCounts+1,capReadCounts=NA)
#genes_t_df$snatch_preference <- genes_t_df$snatchReadCounts/sum(genes_t_df$snatchReadCounts)
genes_t_df$snatch_preference <- 1/nrow(genes_t_df)
genes_t_df$snatch_contribution_unnormalized <- genes_t_df$snatchReadCounts/(genes_t_df$snatchReadCounts+genes_t_df$capReadCounts)
genes_t_df$snatch_contribution_normalized <- genes_t_df$snatch_contribution_unnormalized/sum(genes_t_df$snatch_contribution_unnormalized)
genes_t_df$abundance <- (genes_t_df$snatchReadCounts+genes_t_df$capReadCounts)/sum(genes_t_df$snatchReadCounts+genes_t_df$capReadCounts)
genes_t_df$multireadassignment <- 0.0
genes_t_plus_1_df <- data.frame(row.names=rownames(genes_t_df),
                                snatchReadCounts=genes_t_df$snatchReadCounts,
                                capReadCounts=genes_t_df$capReadCounts,
                                snatch_preference=genes_t_df$snatch_preference,
                                snatch_contribution_unnormalized=genes_t_df$snatch_contribution_unnormalized,
                                snatch_contribution_normalized=genes_t_df$snatch_contribution_normalized,
                                abundance=genes_t_df$abundance,
                                multireadassignment=genes_t_df$multireadassignment)
diff <- 0.01
multimap_readCounts <- as.data.frame(read.table(multimapper_file, header=TRUE))
print(genes_t_df)
while (diff > 0.0001)
{
        genes_t_plus_1_df$multireadassignment <- 0
        for (i in 1:nrow(multimap_readCounts))
        {
                genes <- as.character(multimap_readCounts[i,"genenames"])
                multiplicity <- multimap_readCounts[i,"multiplicity"]
                genes_vector <- unlist(strsplit(genes, split=","))
                sum_snatch_preference <- sum(genes_t_df[genes_vector,"snatch_preference"])
                sum_snatch_contribution <- sum(genes_t_df[genes_vector,"snatch_contribution_normalized"])
                for (gene in genes_vector)
                {
                        #genes_t_plus_1_df[gene,"snatchReadCounts"] <- genes_t_df[gene,"snatchReadCounts"] + (multiplicity * (genes_t_df[gene,"snatch_preference"]/sum_snatch_preference) *  (genes_t_df[gene,"snatch_contribution_normalized"]/sum_snatch_contribution))
                        genes_t_plus_1_df[gene,"multireadassignment"] <- genes_t_plus_1_df[gene,"multireadassignment"] + (multiplicity * (genes_t_df[gene,"snatch_preference"]/sum_snatch_preference))
                }
        }
        genes_t_plus_1_df$snatchReadCounts <- genes_t_df$snatchReadCounts + genes_t_plus_1_df$multireadassignment
        genes_t_plus_1_df$snatch_preference <- genes_t_plus_1_df$snatchReadCounts/sum(genes_t_plus_1_df$snatchReadCounts)
        genes_t_plus_1_df$snatch_contribution_unnormalized <- genes_t_plus_1_df$snatchReadCounts/(genes_t_plus_1_df$snatchReadCounts+genes_t_plus_1_df$capReadCounts)
        genes_t_plus_1_df$snatch_contribution_normalized <- genes_t_plus_1_df$snatch_contribution_unnormalized/sum(genes_t_plus_1_df$snatch_contribution_unnormalized)
        genes_t_plus_1_df$abundance <- (genes_t_plus_1_df$snatchReadCounts+genes_t_plus_1_df$capReadCounts)/sum(genes_t_plus_1_df$snatchReadCounts+genes_t_plus_1_df$capReadCounts)
        diff <- sum(abs(genes_t_plus_1_df$snatch_preference - genes_t_df$snatch_preference))
        genes_t_df <- genes_t_plus_1_df

}
print(genes_t_df)
