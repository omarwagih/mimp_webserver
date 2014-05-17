z = sc[,!names(sc) %in% c('perc_bg', 'perc_fg')]
rownames(z) = NULL

names(z) = c('gene', 'sample_count', 'cancer_type', 'kinase', 'psite_pos', 'pubmed', 'seq_wt', 'seq_mt', 'mut_dist', 'wt_score', 'mt_score', 'pwm', 'mut', 'effect', 'log_ratio')

z$log_ratio = log2(z$mt_score/z$wt_score)
z$kinase[z$kinase == ''] = '.'

write.table(z, '~/Development/mimp/public/R/generate_data/tcga_rewiring_events.tab', quote=F, sep='\t', row.names=F)