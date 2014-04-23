suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MIMP))
options(warn=-1) #suppress warnings


option_list <- list(
  make_option(c("-f", "--fasta"), help="Fasta file" ),
  make_option(c("-m", "--mut"), help="Mutation data", default=NULL),
  make_option(c("-p", "--phos"), help="Phosphorylation data" ),
  make_option(c("-i", "--beta"), help="Lower percentile", type="integer", default=90),
  make_option(c("-j", "--alpha"), help="Upper percentile" , type="integer", default=10),
  make_option(c("-u", "--jobid"), help="Job ID")
)
# Parse the agruments
ARGS=parse_args(OptionParser(option_list = option_list));

# Process without output
sink("/dev/null"); 
data = mimp(muts=ARGS$mut, seqs=ARGS$fasta, psites=ARGS$phos, perc.bg=ARGS$beta, perc.fg=ARGS$alpha, display.results=F)
sink()


if(nrow(data) == 0){
  no_data = file.path("public", "jobs", ARGS$jobid, "no_data");
  writeLines('', no_data)
}else{
  data_path = file.path("public", "jobs", ARGS$jobid, "results.txt");
  cutoffs_path = file.path("public", "jobs", ARGS$jobid, "cutoffs.txt");
  html_path = file.path("public", "jobs", ARGS$jobid, "html.txt");
  write.table(data, data_path, quote=F, sep="\t", row.names=F);
  
  cutoffs = attr(data, "cutoffs")
  js = paste0('"',cutoffs$pwm,'":{','"bg":',cutoffs$bg, ',"fg":', cutoffs$fg, '}', collapse=",")
  js = paste0('cutoffs = {', js, '}');
  writeLines(js, cutoffs_path)
  
  html = dohtml(data, '/assets/images/logos/')
  writeLines(html, html_path)
}


# Output html for java
# writeLines(html);
