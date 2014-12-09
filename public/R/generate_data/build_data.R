source('~/Development/phd/reconstruct_kinase/init.R')
require(parallel)
require(seqinr)


#print(load('neg.ps.data.10000.rsav'))
#print(load('pp.data.new.rsav'))

setwd('~/Desktop/MSc project/rewiring/data/new/')

ps.data = readRDS('psites_combined_mapped_08_2014.rds')
mut.data = readRDS('tcga_mut_data_raw.rds')
sc.data = readRDS('../../results/new/sc_kin_new_2014.rsav')


map = readRDS('gene2nm_longest_iso.rds')
map = unlist(split(map$rseq, map$gene))
mut.data$gene_acc = map[mut.data$gene]

sc.data = sc.data[,c('gene', 'gene_acc', 'n_sample', 'cancer_types', 'mut','kinase', 'position', 'npmid', 'mut.dist', 'flank.wt', 'flank.mt', 'wt.score', 'mt.score', 'pwm', 'effect', 'log.ratio')]
names(sc.data) = gsub('\\.', '_', names(sc.data))
neg.list = readRDS('neg_phosphosites_10k_uniq_noterm.rds')


ks = readRDS('ksr_final_08_2014.rds')
ks_fam = readRDS('ksr_family_final_08_2014.rds')
ks_ref = readRDS('ksr_refined_90p.rds')
ks_fam_ref = readRDS('ksr_family_refined_90p.rds')

ks_newman = readRDS('ksr_newman.rds')
ks_newman = lapply(ks_newman, unique)
ks_newman = ks_newman[sapply(ks_newman, length) >= 10]

aucs = readRDS('~/Desktop/MSc project/rewiring/data/new/auc_cv_kinases_k10_with_both.rds')
aucs_fam = readRDS('~/Desktop/MSc project/rewiring/data/new/auc_cv_kinase_families_k10_both.rds')
auc.cutoff = 0.6
keep = rownames(aucs)[aucs$fam.bg >= auc.cutoff]
ks_ref = ks_ref[keep]

keep = rownames(aucs_fam)[aucs_fam$fam.bg >= auc.cutoff]
ks_fam_ref = ks_fam_ref[keep]

#-----------------------------------------------------

setwd('~/Development/mimp_webserver/public/R/generate_data/')
source('build_func.R')
source('~/Development/mimp/R/Rmimp.R')

# Percentile refine all models
if(F){
  ks_all = tapply(pp[[1]], pp[[2]], c)
  ks_all = cleanKS(ks_all, min.seqs=10)
  names(ks_all) = sapply(names(ks_all), function(i){
    gsub(pattern='\\;|\\-|\\s|\\/|\\|', replacement='_', i)
  })
  # Refine things
  myapply = lapply
  myapply = mclapply
  ks_ref = myapply(names(ks_all), function(n){
    writeLines(sprintf('Processing %s', n))
    percentileRefine(ks_all, neg.list, n, min.seqs=10, perc.cutoff=0.9)
  })
  names(ks_ref) = names(ks_all)
  ks = ks_ref
  saveRDS(ks, 'ks_refined_90.rsav')

}

write.ksr <- function(ksr, file.out){
  ksr = lapply(ksr, function(x) strsplit(paste0(x, collapse=''), '')[[1]])
  write.fasta(ksr, names=names(ksr), file.out=file.out, nbchar=15)
}

write.ksr(ks, 'ks_all.fasta')
write.ksr(ks_ref, 'ks_refined_90.fasta')
write.ksr(ks_fam, 'ks_fam_all.fasta')
write.ksr(ks_fam_ref, 'ks_fam_refined_90.fasta')

z = neg.list
names(z) = c(sprintf('%s negative background psites', names(z)))
write.ksr(z, 'negative_psites.fasta')


stop('x')
# ADD EXCEPTIONS

write.table(sc.data, file = 'tcga_rewiring_events.tab', quote = F, sep = '\t', row.names=F, col.names=T)
write.table(ps.data, file = 'phosphorylation_data.tab', quote = F, sep = '\t', row.names=F, col.names=T)
write.table(mut.data[,c('gene', 'gene_acc', 'cancer_type', 'sample_id', 'wt_residue', 'mut_residue', 'position')], 
            file = 'mutation_data.tab', quote = F, sep = '\t', row.names=F, col.names=T)

# Density distributions in json format for webserver
js <- function(z){
  d = density(z)
  x = signif(d$x,5) 
  y = signif(d$y,5) 
  return( paste0('[{', paste0('"x":', x, ', "y":',y, collapse='},{'), '}]') )
}


maj <- function(x){
  x = substr(x,8,8)
  return(names(which.max(table(x))))
}

getPosNeg <- function(ks, neg, name){
  pos = ks[[name]]
  m = maj(pos)
  if(m == 'S' | m == 'T') neg = neg$ST
  if(m == 'Y') neg = neg$Y
  
  return(list(pos=pos, neg=neg))
}

require(RWebLogo)

processKSR <- function(ks, is.fam=F, is.newman=F){
  set.seed(123)
  
  # Generate PWM models
  writeLines(sprintf('Computing %s PWMs', length(ks)))
  pwms = lapply(ks, PWM)
  
  # Compute background and foreground distributions
  writeLines('Computing distributions')
  dist = lapply(names(pwms), function(n){
    pwm = pwms[[n]]
    z = getPosNeg(ks, neg.list, n)
    fg = mss(z$pos, pwm, na.rm=T)
    bg = mss(z$neg, pwm, na.rm=T)
    list(fg=fg, bg=bg)
  })
  names(dist) = names(pwms)
  
  # Percentile rank functions
  writeLines('Computing percentile rank functions')
  perc_fg = lapply(dist, function(n){
    ecdf(n$fg)
  })
  
  perc_bg = lapply(dist, function(n){
    ecdf(sample(n$bg, 500))
  })
  perc = list(fg=perc_fg, bg=perc_bg)
  
  
  # Compute cutoffs
  writeLines('Computing cutoffs')
  sq = seq(0,1,0.01)
  cutoffs = lapply(dist, function(n){
    a_fg = quantile(n$fg, sq)
    a_bg = quantile(n$bg, sq)
    list(fg=a_fg, bg=a_bg)
  })
  
  writeLines('Combining all data into list ...')
  data = list(pwms=pwms, 
              perc=perc, 
              cutoffs=cutoffs)
  
  writeLines('Computing density distributions')
  folder = ifelse(is.fam, 'density_fam', 'density')
  if(is.newman) folder = 'density_newman'
  dir.create(file.path(folder), recursive=T)
  dens = sapply(names(dist), function(name){
    n = dist[[name]]
    x = paste0('[', js(n$bg), ',', js(n$fg), ']')
    writeLines(x, file.path(folder, paste0(name, '.json')))
  })
  
  logo_folder = ifelse(is.fam, 'logos_fam', 'logos')
  if(is.newman) logo_folder = 'logos_newman'
  dir.create(file.path(folder), recursive=T)
  for(kin in names(ks)){
#     pdf.out = sprintf('%s/%s.pdf', logo_folder, kin)
#     png.out = sprintf('%s/%s.png', logo_folder, kin)
#     weblogo(ks[[kin]], open = F, file.out = pdf.out, yaxis = 4, errorbars = F, title = kin, annotate = -7:7, color.scheme='chemistry3')
#     system(sprintf('convert -density 200 %s %s', pdf.out, png.out))
#     unlink(pdf.out)
    
    
    svg.out = sprintf('%s/%s.svg', logo_folder, kin)
    weblogo(ks[[kin]], open = F, file.out = svg.out, format='svg', yaxis = 4, errorbars = F, title = kin, annotate = -7:7, color.scheme='chemistry3')
  }
  
  return(data)
  
}

if(F){
  
  mimp_data = processKSR(ks_ref, is.fam=F, 
                         auc.file='~/Desktop/MSc project/rewiring/data/new/auc_cv_kinases_k10_with_both.rds')
  saveRDS(mimp_data, file='~/Development/mimp/inst/extdata/mimp_data.rds')
  saveRDS(mimp_data, file='mimp_data.rds')
  
  mimp_data_fam = processKSR(ks_fam_ref, is.fam=T,
                             auc.file='~/Desktop/MSc project/rewiring/data/new/auc_cv_kinase_families_k10_both.rds')
  saveRDS(mimp_data_fam, file='~/Development/mimp/inst/extdata/mimp_data_fam.rds')
  saveRDS(mimp_data_fam, file='mimp_data_fam.rds')
  
  mimp_data_newman = processKSR(ks_newman, is.newman=T)
  saveRDS(mimp_data_newman, file='~/Development/mimp/inst/extdata/mimp_data_newman.rds')
  saveRDS(mimp_data_newman, file='mimp_data_newman.rds')
  
  r_base_logos = '~/Development/mimp/inst/extdata/html/images/logos/'
  web_base_logos = '~/Development/mimp_webserver/public/images/logos/'
  
  r_base_logos_fam = '~/Development/mimp/inst/extdata/html/images/logos_fam/'
  web_base_logos_fam = '~/Development/mimp_webserver/public/images/logos_fam/'
  
  r_base_logos_newman = '~/Development/mimp/inst/extdata/html/images/logos_newman/'
  web_base_logos_newman = '~/Development/mimp_webserver/public/images/logos_newman/'
  
  t = c(r_base_logos, r_base_logos_fam, web_base_logos, web_base_logos_fam, r_base_logos_newman, web_base_logos_newman)
  
  for(dd in t) unlink(dd, recursive = T)
  system(sprintf('cp -r logos %s', r_base_logos))
  system(sprintf('cp -r logos %s', web_base_logos))
  
  system(sprintf('cp -r logos_fam %s', r_base_logos_fam))
  system(sprintf('cp -r logos_fam %s', web_base_logos_fam))
  
  system(sprintf('cp -r logos_newman %s', r_base_logos_newman))
  system(sprintf('cp -r logos_newman %s', web_base_logos_newman))
}


# WRITE MULTIPAGE LOGOS FOR SUPPL 
weblogoBatch(ks_ref, titles = names(ks_ref), file.out = 'kinase_logos_refined.pdf', probability = F)
weblogoBatch(ks_fam_ref, titles = names(ks_fam_ref), file.out = 'kinase_family_logos_refined.pdf', probability = F)
weblogoBatch(ks_newman, titles = names(ks_newman), file.out = 'kinase_logos_newman.pdf', probability = F)

writeLines('All done!')
