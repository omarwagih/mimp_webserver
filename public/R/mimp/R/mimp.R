mimp <-
function(md.file, sd.file, pd.file=NULL, perc.neg=90, perc.pos=10, display.results=T){
  flank=7
  method='MSS'
  MUT_REGEX = '^[A-Z]\\d+[A-Z]$'
  DIG_REGEX = '^\\d+$'
  
  perc.neg = as.integer(perc.neg)
  perc.pos = as.integer(perc.pos)
  if(!is.numeric(perc.neg) | perc.neg > 100 | perc.neg < 0) stop('perc.neg must be an integer between 0-100')
  if(!is.numeric(perc.pos) | perc.pos > 100 | perc.pos < 0) stop('perc.pos must be an integer between 0-100')
  
  # 1. Read sequence data
  writeLines('Reading fasta data from file ...')
  sd = .read.fasta(sd.file, seqtype='AA', as.string=T)
  
  # 2. Read mutation data
  writeLines('Reading mutation data from file ...')
  md = read.table(md.file, header=F, stringsAsFactors=F)
  names(md) = c('gene', 'mut')
  
  # 2. Validate mutation data
  # Ensure valid regex for mutations
  if(!(all(grepl(MUT_REGEX, md$mut) ))) 
    stop('Mutations must follow the following format X123Y for example: A78R!')
  
  # Ensure gene in mutation data is matched to fasta file, if not ignore these rows
  z = setdiff(md$gene, names(sd))
  if(length(z) > 0){
    warning(sprintf('%s genes found in mutation data could not be matched to headers of the fasta file and will be ignored. Genes: %s', 
                    length(z), paste0(z, collapse=', ')))
    md = md[!md$gene %in% z,]
  }
  
  # 3. Read p-site data
  if(is.null(pd.file)){
    # No psite file, generate psites using all STYs
  	writeLines('No phosphosite data found, generating potential phosphosites using all STYs ...')
    sd_sp = strsplit(unlist(sd), '')
    pd = lapply(names(sd_sp), function(n){
      ind = grep('[STY]', sd_sp[[n]])
      data.frame(gene=n, pos=ind, stringsAsFactors=F)
    })
    pd = do.call(rbind, pd)
  }else{
  	writeLines('Reading phosphosites data from file ...')
    pd = read.table(pd.file, header=F, stringsAsFactors=F)
    names(pd) = c('gene', 'pos')
  }
  
  # 3. Validate p-site data
  # Ensure numerical values in p-site position data
  if(!(all(grepl(DIG_REGEX, pd$pos) ))) 
    stop('All positions of p-sites must be non-negative and non-zero!')
  
  # Ensure gene in p-site data is matched to fasta file, if not ignore these rows
  z = setdiff(pd$gene, names(sd))
  if(length(z) > 0){
    warning(sprintf('%s genes found in p-site data could not be matched to headers of the fasta file and will be ignored. Genes: %s', length(z), paste0(z, collapse=', ')))
    pd = pd[!pd$gene %in% z,]
  }
  
  pd$flank.wt = .flankingSequence(sd[pd$gene], pd$pos, flank)
  
  md.ps = .mutateFlanking(md, pd, sd, flank)
  
  
  writeLines('Loading kinase specificity models ...')
  load(system.file("extdata", "pwms.rsav", package = "MIMP"))
  load(system.file("extdata", "cutoffs.rsav", package = "MIMP"))
  load(system.file("extdata", "perc.rsav", package = "MIMP"))
  
  
  
  writeLines('Predicting kinase rewiring events ...')
  pb <- txtProgressBar(min = 0, max = length(pwms), style = 3)
  
  
  ct = t(sapply(names(pwms), function(name){
    c_bg = cutoffs[[name]]$bg[perc.neg+1]
    c_fg = cutoffs[[name]]$fg[perc.pos+1]
    c(c_bg, c_fg)
  }))
  ct = as.data.frame(ct, stringsAsFactors=F)
  ct$pwm = names(pwms)
  colnames(ct) = c('bg', 'fg', 'pwm')
  
  scored = lapply(1:length(pwms), function(i){
    setTxtProgressBar(pb, i)
    pwm = pwms[[i]]
    name = names(pwms)[i]
    md.ps$wt.score = .scorePWM(md.ps$flank.wt, pwm, method)
    md.ps$mt.score = .scorePWM(md.ps$flank.mt, pwm, method)
    md.ps$pwm = name
    
    md.ps = md.ps[!is.na(md.ps$wt.score) & !is.na(md.ps$mt.score),]
    
    #Filter
    c_bg = cutoffs[[name]]$bg[perc.neg+1]
    c_fg = cutoffs[[name]]$fg[perc.pos+1]
    
    a = md.ps$wt.score > c_fg & md.ps$mt.score < c_bg
    b = md.ps$mt.score > c_fg & md.ps$wt.score < c_bg
    
    md.ps[a|b,]
  })
  close(pb)
  
  # Merge all data frames
  s = do.call('rbind', scored)
  attr(s, 'cutoffs') = ct
  
  s$effect = 'Gain'
  s$effect[s$wt.score > s$mt.score] = 'Loss'
  
  s$log.ratio = -log(s$wt.score/s$mt.score)
  
  # Percentile rank
  pfg = unlist( perc$fg[s$pwm] )
  pbg = unlist( perc$bg[s$pwm] )
  
  
  r = sapply(1:nrow(s), function(i){
    e = s$effect[i]
    fg = pfg[[i]]
    bg = pbg[[i]]
    wt = s$wt.score[i]
    mt = s$mt.score[i]
    if(e == "Loss"){
      c( fg(wt), bg(mt))
    }else{
      c( bg(wt), fg(mt))
    }
  })
  r = as.data.frame(t(r))
  names(r) = c('wt.perc', 'mt.perc')
  s$wt.perc = r$wt.perc
  s$mt.perc = r$mt.perc
  
  if(display.results == TRUE) results2html(s)
  
  
  writeLines('Analysis complete!')
  return(s)
}

dohtml <- function(x, LOGO_DIR){
	
    x$wt.score = signif(x$wt.score, 3)
    x$mt.score = signif(x$mt.score, 3)
    x$log.ratio = signif(x$log.ratio, 3)
    x$wt.perc = signif(x$wt.perc, 3)
    x$mt.perc = signif(x$mt.perc, 3)
    
    x$pwm = gsub('-', '_', x$pwm)
    
    n_gain = sum(x$effect == 'Loss')
    n_loss = sum(x$effect == 'Gain')
    n_mut = nrow(unique(x[,c('gene', 'mut')]))
    lines = sapply(1:nrow(x), function(i){
      r = x[i,]
      d = r$mut.dist
      seq = sprintf('%s<br>%s', .htmlSeq(r$flank.wt, d), .htmlSeq(r$flank.mt, d))
      scr = sprintf('%s<br>%s', r$wt.score, r$mt.score)
      prc = sprintf('%s<br>%s', r$wt.perc, r$mt.perc)
      logo = file.path(LOGO_DIR, paste0(r$pwm, '.png'))
      t = sprintf('<a class="hide name">%s</a><img src="%s" class="logo" alt=%s/>', r$pwm, logo, r$pwm)
      
      if(r$effect == 'Loss'){
        eff = '<a class="loss">Loss</a>'
      }else{
        eff = '<a class="gain">Gain</a>'
      }
      
      gene = sprintf('<a target="_blank" class="gene-link" href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s">%s</a>', r$gene, r$gene)
     
      # sub - _ logos
      sprintf('<tr> <td class="gene-name">%s</td>
                    <td class="psite-pos">%s</td>
                    <td class="mut-abbr">%s</td>
                    <td class="mut-dist">%s</td>
                    <td class="sequence">%s</td>
                    <td class="wt-score">%s</td>
                    <td class="mt-score">%s</td>
                    <td class="wt-perc">%s</td>
                    <td class="mt-perc">%s</td>
                    <td class="log-ratio">%s</td>
                    <td class="effect">%s</td>
                    <td class="seq-logo">%s</td>
              </tr>', 
              r$gene, r$pos, r$mut, r$mut.dist, seq, r$wt.score, r$mt.score, r$wt.perc, r$mt.perc, r$log.ratio, eff, t)
     
    })
    
    tt = '<div id="%s" style="display:none">%s</div>'
    
    lines = append(lines,
                   c(sprintf(fmt=tt, 'n_mut', n_mut),
                     sprintf(fmt=tt, 'n_gain', n_gain),
                     sprintf(fmt=tt, 'n_loss', n_loss)))
	return(lines)
}

results2html <-
  function(x){
  	BASE_DIR = system.file("extdata", "", package = "MIMP")
    LOGO_DIR = file.path(BASE_DIR, 'html', 'images', 'logos')
    lines = dohtml(x, LOGO_DIR)
    save = file.path(BASE_DIR, 'html', 'MIMP_results.html')
    zz <- file(save,"w")
    tt = readLines( file.path(BASE_DIR, 'html', 'index.html'))
    ind = grep('<DATA>', tt)
    writeLines(tt[1: (ind-1)],con=zz,sep="\n")
    writeLines(lines,con=zz,sep="\n")
    writeLines(tt[(ind+1): length(tt)],con=zz,sep="\n")
    close(zz)
    
    browseURL(save)
  }
