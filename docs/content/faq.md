# Frequently Asked Questions (FAQ)

## The pipeline starts creating conda environments and crashes (CreateCondaEnvironmentException)
The first thing the pipeline does is making separate conda environments for the rules that will be run. One thing that causes this error is that there isn't enough space available on the device to install these environments. Sometimes, even when there is enough space on the device, the installation still fails, and we haven't been able to pinpoint exactly what is causing this. What usually seems to work is just to remove the installed environemnts (`seq2science clean`) and try again. 

A different type of CreateCondaEnvironmentException occurs when you have conda configured to strict, and this will give a UnsatisfiableError. We haven't been able to solve all problems with this, and recommend to try again with conda set to flexible instead of strict:

```console
user@comp:~$ conda config --set channel_priority flexible
```

## I changed the config/samples file but seq2science does not rerun
TODO: we should support rerunning!
Seq2science (actually Snakemake) has a "lazy" policy regarding the generation of files, and will normally only rerun jobs if the input is younger than the output. 

To push seq2science to do this anyway, you need to remove one or two downstream files. We suggest deleting the MultiQC file, and a fastqc file of the samples that need to rerun.

## Failed to call external services.
TODO explain that downloading can fail on server side. Retry later
