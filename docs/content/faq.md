# Frequently Asked Questions (FAQ)

## One of the rules failed and the log turns red!
Even though we (the developers) think that seq2science works perfectly in all situations, this, unfortunately, is not always the case. Luckily most rules keep a log, where the error message of that rule is being printed to. When a rule fails seq2science prints a couple things, and one of those thing is where the log is stored. Open that file and check the error message. Sometimes rules fail because for instance the storage of your server is full or the internet connection during the downloading of a genome/sample got interrupted. If you think the rule failed because of a mistake on our side, please check our [issues page](https://github.com/vanheeringen-lab/seq2science/issues) and see if someone already made an issue about this. If not, we invite you to start a new issue and we'll try to help you as soon as possible.

p.s. updating to the newest version (if you did not already) also might solve your problem!

## The pipeline starts creating conda environments and crashes (CreateCondaEnvironmentException)
The first thing the pipeline does is making separate conda environments for the rules that will be run. One thing that causes this error is that there isn't enough space available on the device to install these environments. Sometimes, even when there is enough space on the device, the installation still fails, and we haven't been able to pinpoint exactly what is causing this. What usually seems to work is just to remove the installed environemnts (`seq2science clean`) and try again. 

A different type of CreateCondaEnvironmentException occurs when you have conda configured to strict, and this will give a UnsatisfiableError. We haven't been able to solve all problems with this, and recommend to try again with conda set to flexible instead of strict:

```console
user@comp:~$ conda config --set channel_priority flexible
```

## What if I change the configuration or samples file after running seq2science?
Seq2science starts each run with checking if it was already run before, and if so, if any settings were changed. Seq2science then automatically derives which files need to be changed and only starts the necessary jobs.

This comes in handy when for example you want to add more data later on in the analysis or change for instance the colors in the trackhub.

## Failed to call external services.
When downloading samples / looking up their layout online, seq2science makes use of online resources. Sometimes it can happen that those services are not online at the moment of running seq2science, you do not have internet, or you lose connection with the service. Usually just re-running seq2science solves these issues, either directly or a couple hours later.

## Unknown provider
When downloading a genome assembly, seq2science uses genomepy to look up online which assemblies exist, and whether or not the one you try to download can actually be downloaded. In some rare cases this fails, which gives the error `Unknown provider`. This is ussually resolved by rerunning seq2science, either directly or a couple hours later.

## storage exhausted / No Space Left on Device
Some rules in seq2science make a lot of temporary intermediate files. Generally your operating system stores these temporary files on the location of variable `$TMPDIR` which usually is on `/tmp`. Some servers/computers have a small `$TMPDIR` and seq2science does not handle those cases well! You can manually set the tmpdir to another location with commands:

```
mkdir -p /scratch/${USER}/tmp
export TMPDIR=/scratch/${USER}/tmp
```

These commands make a folder tmp on the scratch disc (assuming this disc exists), and set the `$TMPDIR` to there. 

Another, more simple reason why these errors can happen is that there simply is no more space/storage left on the server. This can simply be resolved by deleting unused/old files you no longer need.

## I want to run the pipeline until a certain rule, but not farther
let's say you want to download and trim reads with seq2science, but you do not want to align and do all the other stuff since you have your own methods for that. You can simply tell seq2science to run until some rules, and not farther:

```
seq2science run alignment --cores 48 --snakemakeOptions until=["trim_galore_PE","trim_galore_SE"]
```

## RNA-seq: this gene should be expressed!
Seq2science supports 3 gene-quantifiers: HTSeq-count, FeatureCounts and Salmon.
All three are valid methods, and a large population of the expressed genes is found by all quantifiers.
There are differences between these methods however, and caution is advised with the genes not found by another method.

Additionally, Seq2science employs strict quality checks by default.
You could try reducing the `min_mapping_quality` (for HTSeq and FeatureCounts) or change the parameters to the quantifier.

The biggest difference in expressed genes found is in switching to or from Salmon.
