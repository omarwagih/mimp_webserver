suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rmimp))
options(warn=-1) #suppress warnings


option_list <- list(
  make_option(c("-f", "--fasta"), help="Fasta file" ),
  make_option(c("-m", "--mut"), help="Mutation data", default=NULL),
  make_option(c("-p", "--phos"), help="Phosphorylation data" ),
  make_option(c("-i", "--pthresh"), help="Probability threshold", type="numeric", default=0.5),
  make_option(c("-j", "--logthresh"), help="Log2 threshold" , type="numeric", default=1),
  make_option(c("-c", "--cent"), help="Do central", type="character", default="no"),
  make_option(c("-x", "--mdata"), help="Use kinase family models or predicted models", default="hconf"),
  make_option(c("-u", "--jobid"), help="Job ID")
)
# Parse the agruments
ARGS=parse_args(OptionParser(option_list = option_list));
fam = ifelse(is.null(ARGS$fam), FALSE, TRUE )
ARGS$fam = fam

ARGS$cent = ifelse(ARGS$cent == "no", FALSE, TRUE)

print(ARGS)

# Process without output
sink("/dev/null"); 
data = mimp(muts=ARGS$mut, seqs=ARGS$fasta, psites=ARGS$phos, prob.thresh=ARGS$pthresh, log2.thresh=ARGS$logthresh, 
            display.results=F, include.cent=ARGS$cent, model.data = ARGS$mdata)
sink()

if(is.null(data)) data = data.frame()

if(nrow(data) == 0){
  no_data = file.path("public", "jobs", ARGS$jobid, "no_data");
  writeLines('', no_data)
}else{
  data_path = file.path("public", "jobs", ARGS$jobid, "results.txt");
  html_path = file.path("public", "jobs", ARGS$jobid, "html.txt");

  write.table(format(data, digits=4), data_path, quote=F, sep="\t", row.names=F);
  
  z = sprintf('/assets/R/generate_data/logos%s/', c('', '_fam', '_newman'))
  names(z) = c('hconf', 'hconf-fam', 'lconf')
  logo_path = z[ARGS$mdata]
  
  hl_path = '/assets/R/generate_data/highlight/'
  html = dohtml(data, LOGO_DIR = logo_path, HL_DIR = hl_path)
  writeLines(html, html_path)
}


# Output html for java
# writeLines(html);
