FROM papaemmelab/docker-facets:v0.1.1

COPY . /facets-ch

RUN \
    # Install Facets-CH
    Rscript -e 'devtools::install_local("/facets-ch", build_vignettes=FALSE, upgrade="never")'

ENTRYPOINT [ "Rscript", "/facets-ch/facets-ch.R" ]


