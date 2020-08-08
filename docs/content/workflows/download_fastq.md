## Downloading fastqs
Downloading public data in bulk from the NCBI, ENA, and DDBJ databases has never been easier!

### Pipeline steps
<p align="center">
  <img src="../../_static/download_fastq.png">
</p>

#### Download SRA file
The three largest databases that store sequencing data are National Center for Biotechnology Information (NCBI), the European Nucleotide Archive (ENA) and the DNA Data Bank of Japan (DDBJ). Only the ENA stores the actual fastq files, but all three of them store the raw data (as a sra file) from which a fastq can be derived. For this reason for each sample will first be checked if it can be downloaded from ENA. Otherwise we will download the samples in its raw format. To convert this data to a fastq it has to be *dumped*. 

### Filling out the samples.tsv
Before running a workflow you will have to specify which samples you want to run the workflow on. Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs. As an example, the `samples.tsv` could look something like this:
```
sample
GSM123
GSM321
GSMabc
GSMxzy
GSM890
```

#### Sample column
This column is necessary for all workflows, not just the atac-seq workflow. If you use the pipeline on public data this should be the name of the accession (e.g. GSM2837484), if you use the pipeline on local data this should be the *basename* of the file without the *extension*. For instance `/home/user/myfastqs/sample1.fastq.gz` would be `sample1`.

#### Final notes
- Make sure that the samples.tsv is a tab separated values file when running the pipeline.
- Feel free to add columns to your liking (these will be ignored).

### Best practices
#### Downloading with ascp
ascp is a downloading protocol that allows for (sometimes much) faster speeds. If you have ascp installed you can make the workflow download through this protocol. You can specify the path to the binary in the config.yaml by key `ascp_path` and the ascp key in the config by key `ascp_key`:

See this [gist](https://gist.github.com/mfansler/71f09c8b6c9a95ec4e759a8ffc488be3) for an easy installation of ascp.

#### Filetype extensions and paired-end suffix
The pipeline will save fastq files in the *fastq_dir* directory, which is located in the *result_dir* directory by default.

People and tools have different preferences for storing their data. One of these differences is how to name the fastq.gz files. Some people/tools prefer `fastq.gz`, while others prefer `fq.gz`. By default the pipeline names files with the fastq.gz extension. However if you prefer to change this you can set the variable *fqsuffix* accordingly.
The same goes for paired-end suffix conventions, most people prefer `sample_R1` which is our default. However you can use your preferred suffix by setting variables *fqext1* and *fqext2*.

In the`config.yaml:
```
fastq_dir: ./my_first_fastq_dir

fqsuffix: fastq
fqext1: R1
fqext2: R2
```

#### Downloading many samples
Before the pipeline starts downloading it will first determine whether the samples are paired-end or single-end. This is done by doing two requests per sample with the ncbi tools. Unregistered accounts can only make 3 requests per second. This means that initialization takes +/- 0.66s per sample. If you are an extremely impatient person, like us, you can up these numbers a with a [ncbi api key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) (their default one allows 10 request p/s). You can set the variables *ncbi_key* and *ncbi_requests* (the nr of requests per second) accordingly in the config.yaml.

The results of these lookups get cached, so when re-running parts of the workflow with the same samples no lookup online is required.

#### Fastq-dump options
We made a selection of [reasonable parameters](https://edwards.sdsu.edu/research/fastq-dump/) for the dumping of SRAs to fastq. You can change these settings by changing the variable *splot* for single-end dumping, or the variable *split* for paired-end dumping, in the config.yaml, or on the command line. 
