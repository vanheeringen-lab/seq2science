# Frequently Asked Questions (FAQ)

## The pipeline starts creating conda environments and crashes (CreateCondaEnvironmentException)
The first thing the pipeline does is making separate conda environments for the rules that will be run. The pipeline installs these environments in the (hidden) folder `.snakemake`. One thing that causes this error is that there isn't enough space available on the device to install these environments. Sometimes, even when there is enough space on the device, the installation still fails, and we haven't been able to pinpoint exactly what is causing this. What usually seems to work is just to remove the hidden snakemake folder (`rm -r .snakemake`) and try again. 

## I changed the config/samples file but seq2science does not rerun
TODO: we should support rerunning!
Seq2science (actually Snakemake) has a "lazy" policy regarding the generation of files. 

## Failed to call external services.
TODO explain that downloading can fail on server side. Retry later
