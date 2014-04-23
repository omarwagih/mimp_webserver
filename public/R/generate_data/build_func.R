
#' Find percentile cutoff for a position weight matrix
#' 
#' Compute percentile cutoff using matrix similarity scores of background sequences
#'
#' @param pwm Position weight matrix
#' @param neg.seqs Negative k-mers 
#' @param pos.seqs Positive k-mers 
#' @param perc.bg Percentile of background distribution 
#' @param perc.fg Percentile of foreground distribution (only provide if \code{pos.seqs} not missing)
#' @param is.kinase.pwm True if pwm represents a kinase. Used for central residue exclusion.
#' @keywords cutoff threshold
#' @export
#' @examples
#' # No examples
pwmCutoff <- function(pwm, neg.seqs, pos.seqs, perc.bg=0.9, perc.fg=0.1, is.kinase.pwm=T){
  
  scr.neg = mss(neg.seqs, pwm, is.kinase.pwm, na.rm=T)
  if(missing(pos.seqs)){
    return( list(thresh_bg=quantile(scr.neg, perc.bg)) )
  }
  scr.pos = mss(pos.seqs, pwm, is.kinase.pwm, na.rm=T)
  return( list(thresh_bg = quantile(scr.neg, perc.bg), 
               thresh_fg = quantile(scr.pos, perc.fg)) )
  
}


# Removes minority of ST or Y p-sites using frequency counts from the first argument or second argument (if provided)
#
# Args:
#   seqs: p-site sequences (15mers)
#   seqs2: p-site sequences (15mers)
#
# Returns:
#   Input sequences with minority ST or Y central residue removed
removeMinority <- function(seqs, seqs2){
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

# Majority count of ST or Y
#
# Args:
#   x: p-site sequences (15mers)
#
# Returns:
#   S, T or Y depending on which is most occuring central residue
maj <- function(x){
  x = substr(x,8,8)
  return(names(which.max(table(x))))
}

getPosNeg <- function(ks, neg, name){
  pos = ks[[name]]
  m = .maj(pos)
  if(m == 'S' | m == 'T') neg = neg$ST
  if(m == 'Y') neg = neg$Y
  
  return(list(pos=pos, neg=neg))
}

cleanKS <- function(ks, min.seqs=10){
  writeLines(sprintf('Reading ksr map of %s kinases', length(ks)))
  ks = lapply(ks, unique)
  
  writeLines('Removing minority')
  ks = lapply(ks, removeMinority)
  
  ks = ks[sapply(ks, length) >= min.seqs]
  writeLines(sprintf('Removed kinases with less than %s substrates. %s remaining', min.seqs, length(ks)))
  
  return(ks)
}

percentileRefine <- function(ks, neg.list, name, min.seqs=10, perc.cutoff=.90){
  
  z = getPosNeg(ks, neg.list, name)
  pos = z$pos
  neg = z$neg
  
  pos2 = pos
  
  iter = 1
  while(TRUE){
    pwm = PWM(pos2)
    sc.fg = mss(pos2, pwm, na.rm=T)
    sc.bg = mss(neg, pwm, na.rm=T)
    
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

