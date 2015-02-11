\name{createAtlasSiteConfig}
\alias{createAtlasSiteConfig}
\title{Parse Expression Atlas YAML site config file.}
\description{
    A pacakge to ready in the YAML format site config file which
    contains paths to commonly-used files/directories etc.
}
\usage{
    createAtlasSiteConfig()
}
\value{
    An AtlasSiteConfig object containing a list with elements from the YAML file.
}
\examples{

    # Create a new AtlasSiteConfig object.
    atlasSiteConfig <- createAtlasSiteConfig()

    # Get the list containing the config values.
    siteConfigList <- site_config_list( atlasSiteConfig )

    # Get the path to Ensembl annotation sources.
    annotationSourceDir <- siteConfigList$annotation_source_dir
}
