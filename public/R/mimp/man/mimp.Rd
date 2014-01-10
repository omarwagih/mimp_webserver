\name{mimp}
\alias{mimp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assess the impact of single variants on phosphorylation.
}
\description{
This function takes in mutation, sequence and phosphorylation data to predict the impact the mutation has on phosphorylation.
}
\usage{
mimp(md.file, sd.file, pd.file = NULL, perc.neg = 90, perc.pos = 10, display.results = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{md.file}{
Mutation data file: a space delimited text file containing two columns (1) geneand (1) mutation. Example:
  \tabular{ll}{
    TP53 \tab R282W\cr
    CTNNB1 \tab S33C\cr
    CTNNB1 \tab S37F\cr
  }
}
  \item{sd.file}{
Sequence data file: text file containing protein sequences in FASTA format
}
  \item{pd.file}{
Phosphorylation data file (optional): a space delimited text file containing positions of phosphorylation sites. Example:
  \tabular{ll}{
    TP53 \tab 280\cr
    CTNNB1 \tab 29\cr
    CTNNB1 \tab 44\cr
  }
}
  \item{perc.neg}{
Percentile value between 0 - 100. This value is used to compute a threshold, α- from the negative distribution of scores. By default this is the 90th percentile.
}
  \item{perc.pos}{
Percentile value between 0 - 100. This value is used to compute a threshold, α+ from the positive distribution of scores. By default this is the 10th percentile.
}
  \item{display.results}{
Logical value indicating if results should be displayed in the browser. By default this is \code{TRUE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The data is returned in a \code{data.frame} with the following columns:
%%  If it is a LIST, use
\item{gene}{gene with the rewiring event}
\item{mut}{mutation causing the rewiring event}
\item{pos}{position of the central residue of the phosphosite}
\item{flank.wt}{sequence of the wildtype phosphosite (before the mutation) }
\item{flank.mt}{sequence of the mutated phosphosite (after the mutation)}
\item{mut.dist}{distance of the mutation from the central phosphosite}
\item{wt.score}{matrix similarity score of the wildtype phosphosite }
\item{mt.score}{matrix similarity score of the mutated phosphosite }
\item{pwm}{name of the kinase being rewiried}
\item{effect}{effect the mutation has on phosphorylation}
\item{log.ratio}{log of the wildtype/mutant score}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Get the path to example mutation data 
mut.file = system.file("extdata", "mutation_data.txt", package = "MIMP")

# Get the path to example FASTA sequence data 
seq.file = system.file("extdata", "sequence_data.txt", package = "MIMP")

# View the files in a text editor
browseURL(mut.file)
browseURL(seq.file)

# Run rewiring analysis
results = mimp(mut.file, seq.file, display.results=TRUE)

# Show head of results
head(results)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MIMP }
\keyword{ rewiring }
\keyword{ kinase }
\keyword{ phosphorylation }
