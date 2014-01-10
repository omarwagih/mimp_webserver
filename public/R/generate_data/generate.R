
setwd('~/Development/mimp/public/R/generate_data/')
source('../mimp/R/anRpackage-internal.R')
print(load('neg.ps.data.10000.rsav'))
print(load('pp.data.new.rsav'))
# ks_all = tapply(pp[[1]], pp[[2]], c)
# ks_all = .cleanKS(ks_all, min.seqs=10)
# names(ks_all) = sapply(names(ks_all), function(i){
#   gsub(pattern='\\;|\\-|\\s|\\/|\\|', replacement='_', i)
# })
# # Refine things
# 
# ks_ref = lapply(names(ks_all), function(n){
#   print(n)
#   .percentileRefine(ks_all, neg.list, n, min.seqs=10, perc.cutoff=0.9)
# })
# names(ks_ref) = names(ks_all)
# ks = ks_ref
# save(ks, file='ks_refined_90.rsav')

print(load('ks_refined_90.rsav'))
names(ks) = sapply(names(ks), function(i){
  gsub(pattern='\\;|\\-|\\s|\\/|\\|', replacement='_', i)
})

# Generate PWM models
writeLines(sprintf('Computing %s PWMs', length(ks)))
pwms = lapply(ks, .PWM)

# Compute background and foreground distributions
writeLines('Computing distributions')
dist = lapply(names(pwms), function(n){
  pwm = pwms[[n]]
  z = .getPosNeg(ks, neg.list, n)
  fg = .scorePWM(z$pos, pwm, na.rm=T)
  bg = .scorePWM(z$neg, pwm, na.rm=T)
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


save(perc, file=file.path('perc.rsav'))
save(pwms, file=file.path('pwms.rsav'))
save(cutoffs, file=file.path('cutoffs.rsav'))



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
#plotBatchLogos(ks, '/Users/omarwagih/Desktop/percentile_refine/')