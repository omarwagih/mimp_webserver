

.computeData <-
function(ks, neg.list, save.dir){
  
  # Compute PWM models
  writeLines(sprintf('Computing %s PWMs', length(ks)))
  pwms = lapply(ks, .PWM)
  
  # Compute cutoffs
  sq = seq(0,1,0.01)
  cutoffs = lapply(names(pwms), function(n){
    writeLines(sprintf('Computing cutoffs for %s', n))
    pwm = pwms[[n]]
    z = getPosNeg(ks, neg.list, n)
    a_fg = quantile(.scorePWM(z$pos, pwm, na.rm=T), sq)
    a_bg = quantile(.scorePWM(z$neg, pwm, na.rm=T), sq)
    list(fg=a_fg, bg=a_bg)
  })
  names(cutoffs) = names(pwms)
  
  save(pwms, file=file.path(save.dir, 'pwms.rsav'))
  save(cutoffs, file=file.path(save.dir, 'cutoffs.rsav'))
}

.removeMinority <- function(seqs, seqs2){
  cent = substr(seqs, 8,8)
  if(missing(seqs2)){
    t = table(cent)
    t = names(t)[which.max(t)]
  }else{
    t = table(substr(seqs2, 8,8))
    t = names(t[which.max(t)])
  }
  if(t == 'Y'){
    seqs = seqs[cent == 'Y']
  }else{
    seqs = seqs[cent != 'Y']
  }
  return(seqs)
}

.maj <- function(x){
  x = substr(x,8,8)
  return(names(which.max(table(x))))
}

.getPosNeg <- function(ks, neg, name){
  pos = ks[[name]]
  m = .maj(pos)
  if(m == 'S' | m == 'T') neg = neg$ST
  if(m == 'Y') neg = neg$Y
  
  return(list(pos=pos, neg=neg))
}

.cleanKS <- function(ks, min.seqs=10){
  writeLines(sprintf('Reading ksr map of %s kinases', length(ks)))
  ks = lapply(ks, unique)
  
  writeLines('Removing minority')
  ks = lapply(ks, .removeMinority)
  
  ks = ks[sapply(ks, length) >= min.seqs]
  writeLines(sprintf('Removed kinases with less than %s substrates. %s remaining', min.seqs, length(ks)))
  
  return(ks)
}

.percentileRefine <- function(ks, neg.list, name, min.seqs=10, perc.cutoff=.90){
  
  z = .getPosNeg(ks, neg.list, name)
  pos = z$pos
  neg = z$neg
  
  pos2 = pos
  
  iter = 1
  while(TRUE){
    pwm = .PWM(pos2)
    sc.fg = .scorePWM(pos2, pwm, na.rm=T)
    sc.bg = .scorePWM(neg, pwm, na.rm=T)
    
    cut = quantile(sc.bg, perc.cutoff)
    print(cut)
    
    dump = names(sc.fg[sc.fg < cut])
    
    writeLines(sprintf('Iteration: %s\n\tCutoff: %s\n\tDump size: %s', iter, cut, length(dump)))
    
    if(length(dump) == 0)  break
    if(length(pos2) - length(dump) < min.seqs) break
    # remove the dump
    pos2 = setdiff(pos2, dump)
    iter = iter +1
  }
  n = signif((length(pos)-length(pos2))/length(pos)*100, 2)
  writeLines(sprintf('Refining for %s done, %s perc. lost', name,n))
  attr(pos2, 'percentLost') = n
  attr(pos2, 'iterations') = iter
  return(pos2)
}


.flankingSequence <- function(seqs, inds, flank=7, empty.char='-'){
  if(length(seqs) == 1 & length(inds) >= length(seqs)) seqs = rep(seqs, length(inds))
  if(length(seqs) != length(inds)) stop('Length of sequences must be equal to length of positions')
  
  border = paste0(rep(empty.char, flank), collapse="")
  seqs = sapply(seqs, function(s) paste0(border,s,border))
  
  ret = sapply(1:length(seqs), function(i){
    s = seqs[i]
    p = flank + inds[i]
    substr(s, p-flank, p+flank)
  })
  return(ret)
}
.htmlSeq <- function(s, dist){
  s = strsplit(s,'')[[1]]
  s[8] = sprintf('<a class="psite">%s</a>', s[8])
  p = 8 + dist
  s[p] = sprintf('<a class="mut">%s</a>', s[p])
  paste0(s, collapse='')
}

.mutateFlanking <- function(md, pd, sd, flank=7){
  psites = split(pd, pd$gene)
  n = nrow(md)
  
  tt = t( sapply(strsplit(md$mut, ''), function(s){
    c(s[1], s[length(s)], paste0(s[2:(length(s)-1)],collapse=''))
  }) )
  tt = as.data.frame(tt, stringsAsFactors=F)
  names(tt) = c('mut.orig', 'mut.change', 'mut.pos')
  tt$mut.pos = as.numeric(tt$mut.pos)
  md = cbind(md, tt)
  
  mut.ps = lapply(1:n, function(i){
    prot.name = md$gene[i]
    mut.pos = md$mut.pos[i]
    
    prot.psites = psites[[prot.name]]
    
    # Check if the mutation exists in any of the psites 
    t = ( mut.pos <= (prot.psites$pos + flank) 
          & mut.pos >= (prot.psites$pos - flank) )
    if(sum(t) == 0){
      NULL
    }else{
      # Mutate sequence
      mut.seq = .replaceStr(sd[[prot.name]], md[i,]$mut.change, mut.pos)
      # Append mutation information
      prot.psites = prot.psites[t,]
      # Append mutated flanking sequence
      prot.psites$flank.mt = .flankingSequence(mut.seq, prot.psites$pos, flank)
      # Bind everything
      suppressWarnings( ret <- cbind(md[i,], prot.psites) )
      ret$mut.dist = mut.pos - prot.psites$pos
      ret
    }
  })
  
  mut.ps = mut.ps[!sapply(mut.ps, is.null)]
  mut.ps = do.call('rbind', mut.ps)
  mut.ps = mut.ps[, setdiff(names(mut.ps), c(names(tt), 'gene.1'))]
  return(mut.ps)
}

.PWM <- function(seqs, pseudocount=0.001, log.bg=F, relative.freq=T, type='AA', 
                priors=c(A=0.070, R=0.056, N=0.036, D=0.048,C=0.023,Q=0.047,E=0.071,G=0.066,H=0.026,I=0.044,
                         L=0.100,K=0.058,M=0.021,F=0.037,P=0.063,S=0.083,T=0.053,W=0.012,Y=0.027,V=0.060), 
                pwm.prior=numeric(0) ){
  
  # Ensure same length characters 
  seq.len = sapply(seqs, nchar)
  num.pos = seq.len[1]
  if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  
  # Type validity
  if(!type %in% c('AA', 'DNA')) stop('Type must be AA or DNA')
  
  # List of valid amino acids, sorted
  aa = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
  if(type == 'DNA') aa = c('A', 'C', 'G', 'T')
  
  # Validate priors, if is na priors are equal
  if(any(is.na(priors))){
    priors = rep(1/length(aa), length(aa))
    names(priors) = aa
  }else{
    # We have priors, ensure they match to namespace
    if(any(is.na( match(names(priors), aa) ))) stop('Please provide priors for all AA or DNA')
  }
  # Match priors to aa 
  bg.prob = priors[match(aa, names(priors))]
  # Make matrix of letters
  split = unlist( sapply(seqs, function(seq){strsplit(seq, '')}) )
  m = t( matrix(split, seq.len, length(split)/num.pos) )
  
  i = 1
  # Construct PWM
  pwm.matrix = apply(m, 2, function(pos.data){
    # Get frequencies 
    t = table(pos.data)
    # Match to aa
    ind = match(aa, names(t))
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    
    if(length(pwm.prior) > 0){
      bg.prob = pwm.prior[,i]
      #print(i)
      #print(bg.prob)
    }
    
    # Do pseudocounts
    col = col + (pseudocount * bg.prob) #Pseudocounts
    
    # Do relative frequencies
    if(relative.freq){
      col = col / sum(col)
    }
    if(log.bg){
      col = col * log(col  / bg.prob) 
    }
    i = i + 1
    col
  })
  
  attr(pwm.matrix, 'pseudocount') = pseudocount
  attr(pwm.matrix, 'log.bg') = log.bg
  ic = apply(pwm.matrix, 2, function(col){
    #sum(col * log(col))
    sum(col* ( log(20*col)), na.rm=T)
  })
  ic2 = apply(pwm.matrix, 2, function(col){
    #sum(col * log(col))
    sum(col* ( log(col/bg.prob, base=2)), na.rm=T)
  })
  attr(pwm.matrix, 'ic') = ic
  attr(pwm.matrix, 'ic2') = ic2
  
  # Assign AA names to rows/pos col
  rownames(pwm.matrix) = aa
  colnames(pwm.matrix) = 1:num.pos
  return(pwm.matrix)
}
.replaceStr <- function(seq, replacement, index){
  sp = strsplit(seq, '')
  ret = sapply(sp, function(s){
    s[index] = replacement
    paste(s, collapse='')
  })
  return(ret)
}
.scoreArray <- function(seqs, pwm){
  # Split sequence
  sp = strsplit(seqs, '')
  
  seq.lens = sapply(sp, length)
  seq.len = seq.lens[1]
  if(any(seq.lens != seq.len)) stop('Input sequences must be same length')
  
  # Iterate through sequences
  dat = lapply(sp, function(seq){
    # Match sequence to the PWM
    mat = matrix(c(match(seq, rownames(pwm)), 1:seq.len), seq.len, 2)
    prob.vector = pwm[ mat ]
    prob.vector
  })
  names(dat) = seqs
  return(dat)
}

.scoreLOG <- function(seqs, pwm, ignore.central=T, central.res="*"){
  # Get central residue index
  central.ind = ceiling(ncol(pwm)/2)
  
  score.arr = .scoreArray(seqs, pwm)
  scores = sapply(score.arr, function(sa){
    na = is.na(sa)
    na[central.ind] = ignore.central
    
    score.final = sum( sa[!na] , na.rm=T ) 
    score.final
  })
  
  keep = grepl(central.res, substr(seqs, central.ind, central.ind))
  scores[!keep] = NA
  return(scores)
}

.scoreMSS <- function(seqs, pwm, ignore.central=T, central.res="*"){
  
  # Central residue index
  central.ind = ceiling(ncol(pwm)/2)
  # Best/worst sequence match
  oa = .scoreArray(.bestSequence(pwm), pwm)[[1]]
  wa = .scoreArray(.worstSequence(pwm), pwm)[[1]]
  
  # Information content for MSS method
  #   I = apply(pwm, 2, function(col){ 
  #     sum(col* ( log(20*col) ) , na.rm=T)
  #   })
  I = attr(pwm, 'ic2')
  
  score.arr = .scoreArray(seqs, pwm)
  scores = sapply(score.arr, function(sa){
    na = is.na(sa)
    na[central.ind] = ignore.central
    
    # Get information content of non-NA values
    IC = I[!na]
    
    # MSS method
    curr.score  = sum( IC * (sa [!na]), na.rm=T ) 
    opt.score   = sum( IC * (oa [!na]), na.rm=T )
    worst.score = sum( IC * (wa [!na]), na.rm=T )
    score.final = ( (curr.score - worst.score) / (opt.score - worst.score) )
    score.final
  })
  
  keep = grepl(central.res, substr(seqs, central.ind, central.ind))
  scores[!keep] = NA
  return(scores)
}

.scorePWM <- function(seqs, pwm, method='MSS', ignore.central=TRUE, na.rm=F){
  method = toupper(method)
  
  log.bg = attr(pwm, 'log.bg') 
  if(is.null(log.bg)) log.bg = TRUE
  if(method == 'LOG' & !log.bg) warning('LOG method chosen but PWM not constructed based on background frequencies')
  
  # Only score sequences which have a central residue S/T or Y depending on the PWM
  kinase.type = names(which.max(pwm[,ceiling(ncol(pwm)/2)]))
  if(grepl('S|T', kinase.type)){kinase.type = 'S|T'} else{kinase.type='Y'}
  
  # Do scoring
  if(! method %in% c('MSS', 'LOG')) stop(paste(method, ': unknown scoring method'))
  if(method == "MSS") ret = .scoreMSS(seqs, pwm, ignore.central, central.res=kinase.type)
  if(method == "LOG") ret = .scoreLOG(seqs, pwm, ignore.central, central.res=kinase.type)
  
  # Remove NA if requested
  if(na.rm) ret = ret[!is.na(ret)]
  return(ret)
}

.worstSequence <- function(pwm){
  w = rownames(pwm)[apply(pwm, 2, which.min)]
  return(paste(w, collapse=''))
}
.bestSequence <- function(pwm){
  b = rownames(pwm)[apply(pwm, 2, which.max)]
  return(paste(b, collapse=''))
}


.read.fasta = function (file, 
          seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE, 
          set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, 
          strip.desc = FALSE, bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong, 
          endian = .Platform$endian, apply.mask = TRUE) 
{
  seqtype <- match.arg(seqtype)
  if (!bfa) {
    lines <- readLines(file)
    if (legacy.mode) {
      comments <- grep("^;", lines)
      if (length(comments) > 0) 
        lines <- lines[-comments]
    }
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseq <- length(ind)
    if (nseq == 0) {
      stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                         collapse = ""))
    if (seqonly) 
      return(sequences)
    nomseq <- lapply(seq_len(nseq), function(i) {
      firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
      substr(firstword, 2, nchar(firstword))
    })
    if (seqtype == "DNA") {
      if (forceDNAtolower) {
        sequences <- as.list(tolower(sequences))
      }
    }
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        Annot <- lines[ind[i]]
        if (strip.desc) 
          Annot <- substr(Annot, 2L, nchar(Annot))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = switch(seqtype, AA = "SeqFastaAA", 
                                                                         DNA = "SeqFastadna"))
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
  if (bfa) {
    if (seqtype != "DNA") 
      stop("binary fasta file available for DNA sequences only")
    mycon <- file(file, open = "rb")
    r2s <- words(4)
    readOneBFARecord <- function(con, sizeof.longlong, endian, 
                                 apply.mask) {
      len <- readBin(con, n = 1, what = "int", endian = endian)
      if (length(len) == 0) 
        return(NULL)
      name <- readBin(con, n = 1, what = "character", endian = endian)
      ori_len <- readBin(con, n = 1, what = "int", endian = endian)
      len <- readBin(con, n = 1, what = "int", endian = endian)
      seq <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                     size = 1, endian = endian)
      mask <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                      size = 1, endian = endian)
      if (endian == "little") {
        neword <- sizeof.longlong:1 + rep(seq(0, (len - 
                                                    1) * sizeof.longlong, by = sizeof.longlong), 
                                          each = sizeof.longlong)
        seq <- seq[neword]
        mask <- mask[neword]
      }
      seq4 <- c2s(r2s[as.integer(seq) + 1])
      seq4 <- substr(seq4, 1, ori_len)
      if (apply.mask) {
        mask4 <- c2s(r2s[as.integer(mask) + 1])
        mask4 <- substr(mask4, 1, ori_len)
        npos <- gregexpr("a", mask4, fixed = TRUE)[[1]]
        for (i in npos) substr(seq4, i, i + 1) <- "n"
      }
      return(list(seq = seq4, name = name))
    }
    sequences <- vector(mode = "list")
    nomseq <- vector(mode = "list")
    i <- 1
    repeat {
      res <- readOneBFARecord(mycon, sizeof.longlong, endian, 
                              apply.mask)
      if (is.null(res)) 
        break
      sequences[[i]] <- res$seq
      nomseq[[i]] <- res$name
      i <- i + 1
    }
    close(mycon)
    nseq <- length(sequences)
    if (seqonly) 
      return(sequences)
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        if (!strip.desc) 
          Annot <- c2s(c(">", nomseq[[i]]))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = "SeqFastadna")
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
}
