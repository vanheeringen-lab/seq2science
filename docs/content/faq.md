# Frequently Asked Questions (FAQ)

## One of the rules failed and the log turns red!
Even though we (the developers) think that seq2science works perfectly in all situations, this is, unfortunately, not always the case. Luckily most rules keep a log, where the error message of that rule is being printed to. When a rule fails seq2science prints a couple of things, and one of those things is where the log is stored. Open that file and check the error message. Sometimes rules fail because for instance the storage of your server is full or the internet connection during the downloading of a genome/sample got interrupted. If you think the rule failed because of a mistake on our side, please check our [issues page](https://github.com/vanheeringen-lab/seq2science/issues) and see if someone already made an issue about this. If not, we invite you to start a new issue and we'll try to help you as soon as possible.

p.s. updating to the newest version (if you did not already) also might solve your problem!

## The pipeline starts creating conda environments and crashes (CreateCondaEnvironmentException)
The first thing the pipeline does is make separate conda environments for the rules that will be run. One thing that causes this error is that there isn't enough space available on the device to install these environments. Sometimes, even when there is enough space on the device, the installation still fails, and we haven't been able to pinpoint exactly what is causing this. What usually seems to work is just to remove the installed environments (`seq2science clean`) and try again. 

A different type of CreateCondaEnvironmentException occurs when you have conda configured to strict, and this will give a UnsatisfiableError. We haven't been able to solve all problems with this, and recommend trying again with conda set to flexible instead of strict:

```console
user@comp:~$ conda config --set channel_priority flexible
```
## Encountered problems while solving (nothing provides ... needed by seq2science) 

This error usually has something to do with not having all the channels properly added. Running this should resolve it:

```console
user@comp:~$ conda config --add channels defaults
user@comp:~$ conda config --add channels bioconda
user@comp:~$ conda config --add channels conda-forge
```

## What if I change the configuration or samples file after running seq2science?
Seq2science starts each run by checking if it was already run before, and if so, if any settings were changed. Seq2science then automatically derives which files need to be changed and only starts the necessary jobs.

This comes in handy when for example you want to add more data later on in the analysis or change for instance the colors in the trackhub.

## Failed to call external services.
When downloading samples / looking up their layout online, seq2science makes use of online resources. Sometimes it can happen that those services are not online at the moment of running seq2science, you do not have internet, or you lose connection with the service. Usually just re-running seq2science solves these issues, either directly or a couple of hours later.

## Unknown provider
When downloading a genome assembly, seq2science uses genomepy to look up online which assemblies exist, and whether or not the one you try to download can actually be downloaded. In some rare cases, this fails, which gives the error `Unknown provider`. This is usually resolved by rerunning seq2science, either directly or a couple of hours later.

## Storage exhausted / No Space Left on Device
Some rules in seq2science make a lot of temporary intermediate files. Generally, your operating system stores these temporary files on the location of variable `$TMPDIR` which usually is on `/tmp`. Some servers/computers have a small `$TMPDIR` and seq2science does not handle those cases well! You can manually set the tmpdir to another location with the commands:

```
mkdir -p /scratch/${USER}/tmp
export TMPDIR=/scratch/${USER}/tmp
```

These commands make a folder tmp on the scratch disc (assuming this disc exists), and set the `$TMPDIR` to there. 

Another, more simple reason why these errors can happen is that there simply is no more space/storage left on the server. This can be resolved by deleting unused/old files you no longer need.

## I want to run the pipeline until a certain rule, but not farther
let's say you want to download and trim reads with seq2science, but you do not want to align and do all the other stuff since you have your own methods for that. You can tell seq2science to run until some rules, and no further:

```
seq2science run alignment --cores 48 --snakemakeOptions until=["trim_galore_PE","trim_galore_SE"]
```

## My run is stuck on 'Select jobs to execute...' What should I do?
Snakemake has a job scheduler to manage the execution of jobs within a workflow, ensuring they are executed in the correct order and maintaining workflow integrity. By default, Snakemake utilizes a non-greedy solver, which prioritizes global optimization over immediate job availability. However, this approach may result in delays, particularly when working with a large number of samples in an analysis.

To address this, you can switch to a greedy solver by using the '--snakemakeOptions' flag when running the seq2science run alignment command. By specifying 'scheduler=greedy' within the options, Snakemake will prioritize the execution of jobs based on their immediate availability. Here's an example of how to use this command:

```
seq2science run alignment --cores 48 --snakemakeOptions scheduler=greedy
```

## ChIP-seq: MACS2: Too few paired peaks so I can not build the model!
This error is due to the fact that MACS2 cannot properly make its internal shifting model, it wants to shift since a peak might not have the same signal on both strands of the genome.
This generally shouldn't happen with good samples with peaks/enrichment.
Did you perhaps add the input/control as a sample?
The input/control shouldn't really have any peaks, so macs2 will fail on those samples here.
You should add the input/control in the `control` column of the samples.tsv, see the ChIP-seq docs.

## RNA-seq: this gene should be expressed!
Seq2science supports 3 gene quantifiers: HTSeq-count, FeatureCounts, and Salmon.
All three are valid methods, and a large population of the expressed genes is found by all quantifiers.
There are differences between these methods, however, and caution is advised with the genes not found by another method.

Additionally, Seq2science employs strict quality checks by default.
You could try reducing the `min_mapping_quality` (for HTSeq and FeatureCounts) or changing the parameters to the quantifier.

The biggest differences in expressed genes can be found by switching gene annotation, followed by switching quantifier (to or from Salmon).

## scRNA-seq: Kallisto|Bustools: \[~warn\] no reads pseudoaligned
This error means that no reads were (pseudo)aligned to the reference genome! 
Oh no! 
This is probably because you did not specify the chemistry correctly. 
Take a look at how to [specify the BUS format](https://vanheeringen-lab.github.io/seq2science/content/workflows/scrna_seq.html#bus-barcode-umi-set-format) for kallisto|bustools.

## scRNA-seq: Kallisto|Bustools: what():  std::bad_alloc
bad alloc generally means that kallisto|bustools tried to reserve some memory, but there was none left to be reserved...
Are there other programs running that take up a lot of memory?
For us (the developers) it also seems as if this rule sometimes happens, seemingly at random.
Simply restarting seq2science might just solve the problem!

## UCSC trackhub displays no data
Sometimes this happens because even though the genome assembly used is supported by UCSC, somehow the chromosome ids are different. 
For instance chr1 vs Chr1 vs 1, etc.
To solve this one can try to change assembly, with the downside that all the alignment, etc. has to be re-run again.
Or you can add `force_assembly_hub: true` in the config so that only the trackhub gets remade, but this time as a so-called assembly hub.
An assembly hub does not use a genome assembly that's provided by UCSC so that no chromosome id mismatches can occur.

## ATAC-seq: bwa-mem2: Runs are aborted due to a segmentation fault. 
If you get an error while running the rules bwa-mem2/samtools_presort, yet the log bwa-mem2_align and samtools_presort log files show no errors, this may be due to segmentation faults that sometimes occur when running bwa-mem2. This can be circumvented by using the aligner bwa-mem instead.
