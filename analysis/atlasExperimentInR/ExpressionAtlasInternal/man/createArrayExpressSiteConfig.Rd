\name{createArrayExpressSiteConfig}
\alias{createArrayExpressSiteConfig}
\title{Parse ArrayExpress YAML site config file.}
\description{
    A function to read in the YAML format site config file which
    contains paths to commonly-used files/directories etc.
}
\usage{
    createArrayExpressSiteConfig()
}
\value{
    A list with elements from the YAML file.
}
\examples{
\dontrun{
    # Create a new ArrayExpressSiteConfig object.
    arrayExpressSiteConfig <- createArrayExpressSiteConfig()

    # Get the path to ArrayExpress load directory.
    aeLoadDir <- arrayExpressSiteConfig$get_AE2_LOAD_DIR
}
}
