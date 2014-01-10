suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MIMP))
options(warn=-1) #suppress warnings


option_list <- list(
    make_option(c("-f", "--fasta"), help="Fasta file" ),
    make_option(c("-m", "--mut"), help="Mutation data", default=NULL),
    make_option(c("-p", "--phos"), help="Phosphorylation data" ),
    make_option(c("-i", "--amin"), help="Lower percentile", type="integer", default=5),
    make_option(c("-j", "--amax"), help="Upper percentile" , type="integer", default=5),
    make_option(c("-u", "--jobid"), help="Job ID")
)

# Parse the agruments
args=parse_args(OptionParser(option_list = option_list));

# Process without output
sink("/dev/null"); 
data = mimp(args$mut, args$fasta, pd.file=args$phos, perc.neg=args$amin, perc.pos=args$amax, display.results=F)
sink();

data_path = file.path("public", "jobs", args$jobid, "results.txt");
cutoffs_path = file.path("public", "jobs", args$jobid, "cutoffs.txt");
html_path = file.path("public", "jobs", args$jobid, "html.txt");
write.table(data, data_path, quote=F, sep="\t", row.names=F);

cutoffs = attr(data, "cutoffs")
js = paste0('"',cutoffs$pwm,'":{','"bg":',cutoffs$bg, ',"fg":', cutoffs$fg, '}', collapse=",")
js = paste0('cutoffs = {', js, '}');
writeLines(js, cutoffs_path)

html = dohtml(data, '/assets/images/logos/');
writeLines(html, html_path);

# Output html for java
writeLines(html);
