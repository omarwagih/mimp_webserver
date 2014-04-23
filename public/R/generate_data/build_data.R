require(parallel)
require(seqinr)

VERSION = '1.0'

setwd('~/Development/mimp/public/R/generate_data/')
source('build_func.R')
source('../mimp/R/Rmimp.R')

detachPackage <- function(pkg){
  pkg = sprintf("package:%s", pkg)
  res = tryCatch({
    detach(pkg, unload=TRUE, character.only=T, force=T)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(res)
}

build_package <- function(){
  require(devtools)
  # Compile things
  document('../mimp/')
  system('R CMD BUILD ../mimp/')
  targz = sprintf('MIMP_%s.tar.gz', VERSION)
  # Move up one directory
  newf = file.path('..', targz)
  re = file.rename(targz, newf)
  # Install package
  system(sprintf('R CMD INSTALL %s', newf))
  # Reload in current environment
  detachPackage('MIMP')
  require(MIMP)
}
stop('x')
print(load('neg.ps.data.10000.rsav'))
print(load('pp.data.new.rsav'))



# Percentile refine all models
if(T){
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
  
  # Write fasta files
  ks_all_fa = lapply(ks_all, function(x) strsplit(paste0(x, collapse=''), '')[[1]])
  write.fasta(ks_all_fa, names=names(ks_all), file.out='ks_all.fasta', nbchar=15)
  
  ks_refined_fa = lapply(ks, function(x) strsplit(paste0(x, collapse=''), '')[[1]])
  write.fasta(ks_refined_fa, names=names(ks), file.out='ks_refined_90.fasta', nbchar=15)
  
  
  ks_refined_fa = lapply(ks, function(x) strsplit(paste0(x, collapse=''), '')[[1]])
  write.fasta(ks_refined_fa, names=names(ks), file.out='ks_refined_90.fasta', nbchar=15)
  
  
  z = neg.list
  names(z) = c(sprintf('%s negative background psites', names(z)))
  z = lapply(z, function(x) strsplit(paste0(x, collapse=''), '')[[1]])
  write.fasta(z, names=names(z), file.out='negative_psites.fasta', nbchar=15)
}

# Read in refined sequences
ks = readRDS('ks_refined_90.rsav')

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

saveRDS(data, file='../mimp/inst/extdata/mimp_data.rds')


# Density distributions in json format for webserver
js <- function(z){
  d = density(z)
  x = signif(d$x,5) 
  y = signif(d$y,5) 
  return( paste0('[{', paste0('"x":', x, ', "y":',y, collapse='},{'), '}]') )
}

writeLines('Computing density distributions')
dir.create(file.path('density'), recursive=T)
dens = sapply(names(dist), function(name){
  n = dist[[name]]
  x = paste0('[', js(n$bg), ',', js(n$fg), ']')
  writeLines(x, file.path('density', paste0(name, '.json')))
})


# Build package
build_package()

stop('Done')

bestPercSample <- function(x, n.samples=1000, sample.size=500){
  v = seq(0,1,0.01)
  actual = ecdf(x)(v)
  
  # Generate the samples
  samples = lapply(1:n.samples, function(i){
    sample(x, sample.size)
  })
  
  # Sum squared difference
  sse = sapply(samples, function(s){
    approx = ecdf(s)(v)
    sum( (approx - actual)^2 )
  })
  
  best = samples[[which.min(sse)]]
  worst = samples[[which.max(sse)]]
  
  return(best)
}

# require(RWebLogo)
#plotBatchLogos(ks, '/Users/omarwagih/Desktop/percentile_refine/')