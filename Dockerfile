FROM quay.io/jupyter/r-notebook
COPY cliquePartitioneR cliquePartitioneR

RUN Rscript -e "devtools::install_github('https://github.com/p100mma/HCCsim')"
RUN Rscript -e "devtools::install_local('cliquePartitioneR')"
