suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rmimp))
options(warn=-1) #suppress warnings


option_list <- list(
  make_option(c("-f", "--fasta"), help="Fasta file" ),
  make_option(c("-m", "--mut"), help="Mutation data", default=NULL),
  make_option(c("-p", "--phos"), help="Phosphorylation data" ),
  make_option(c("-i", "--beta"), help="Lower percentile", type="integer", default=90),
  make_option(c("-j", "--alpha"), help="Upper percentile" , type="integer", default=10),
  make_option(c("-x", "--mdata"), help="Use kinase family models or predicted models", default="hconf"),
  make_option(c("-u", "--jobid"), help="Job ID")
)
# Parse the agruments
ARGS=parse_args(OptionParser(option_list = option_list));
fam = ifelse(is.null(ARGS$fam), FALSE, TRUE )
ARGS$fam = fam
print(ARGS)

# Process without output
sink("/dev/null"); 
data = mimp(muts=ARGS$mut, seqs=ARGS$fasta, psites=ARGS$phos, perc.bg=ARGS$beta, perc.fg=ARGS$alpha, 
            display.results=F, include.cent=T, model.data = ARGS$mdata)
sink()

if(is.null(data)) data = data.frame()

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
  
  z = sprintf('/assets/R/generate_data/logos%s/', c('', '_fam', '_newman'))
  names(z) = c('hconf', 'hconf-fam', 'lconf')
  logo_path = z[ARGS$mdata]
  
  hl_path = '/assets/R/generate_data/highlight/'
  html = dohtml(data, LOGO_DIR = logo_path, HL_DIR = hl_path)
  writeLines(html, html_path)
}


# Output html for java
# writeLines(html);
