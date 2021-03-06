% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_hap_files.R
\name{update_hap_files}
\alias{update_hap_files}
\title{Edit the (upzipped) haplotype files}
\usage{
update_hap_files(con1, con2, repl_lines, line_nos)
}
\arguments{
\item{con1}{name of the original unzipped haplotype file}

\item{con2}{name of the new haplotype file, also unzipped}

\item{repl_lines}{replacement lines - obtained from readLines of the file to replace}

\item{line_nos}{line numbers to be replaced}
}
\value{
R output is text stating what has been accomplished. The new files have been written.
}
\description{
After creating new pruned variants with the *create_new_mac* function,
edit the unzipped haplotypes. A new haplotype file is created.
This needs to be done before variants are removed to ensure the line numbers align correctly.
}
