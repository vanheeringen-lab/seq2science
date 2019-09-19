def get_strandedness(wildcards):
    sample = f"{wildcards.sample}"
    strandedness = samples["strandedness"].loc[sample]
    return strandedness

rule bam_stranded_bigwig:
    """
    Convert a bam file into two bigwig files, one for each strand    
    """
    input:
        bam=expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config),
        bai=expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bai", **config)
    output:
        forward=expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.fwd.bw", **config),
        reverse=expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.rev.bw", **config),
    params:
        flags=config['bam_bigwig']['deeptools'],
        strandedness=get_strandedness
    wildcard_constraints:
        sample=any_sample(),
        sorting=config['bam_sort_order']
    log:
        expand("{log_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 20
    resources:
        deeptools_limit=1
    shell:
        """       
        direction1=forward
        direction2=reverse
        if [ {params.strandedness} == 'reverse' ]; then
            direction1=reverse
            direction2=forward
        fi
                    
        bamCoverage --bam {input.bam} --outFileName {output.forward} --filterRNAstrand $direction1 --numberOfProcessors {threads} {params.flags} --verbose >> {log} 2>&1 &&        
        bamCoverage --bam {input.bam} --outFileName {output.reverse} --filterRNAstrand $direction2 --numberOfProcessors {threads} {params.flags} --verbose >> {log} 2>&1
        """

rule bam_bigwig:
    """
    Convert a bam file into a bigwig file
    """
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bw", **config)
    params:
        config['bam_bigwig']['deeptools']
    wildcard_constraints:
        sample=any_sample(),
        sorting=config['bam_sort_order']
    log:
        expand("{log_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 20
    resources:
        deeptools_limit=1
    shell:
        """
        bamCoverage --bam {input} --outFileName {output} --numberOfProcessors {threads} {params} --verbose >> {log} 2>&1
        """


# rule narrowpeak_bignarrowpeak:
#     """
#     Convert a narrowpeak file into a bignarrowpeak file.
#     """
#     input:
#         narrowpeak= expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_peaks.narrowPeak", **config),
#         genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
#     output:
#         out=     expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}.bigNarrowPeak", **config),
#         tmp=temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}.tmp.narrowPeak", **config))
#     log:
#         expand("{log_dir}/narrowpeak_bignarrowpeak/{{peak_caller}}/{{assembly}}-{{sample}}.log", **config)
#     benchmark:
#         expand("{benchmark_dir}/bedgraphish_to_bedgraph/{{assembly}}-{{sample}}-{{peak_caller}}.log", **config)[0]
#     conda:
#         "../envs/ucsc.yaml"
#     shell:
#         """
#         # keep first 10 columns, idr adds extra columns we do not need for our bigpeak
#         cut -d$'\t' -f 1-10 {input.narrowpeak} |
#         bedSort /dev/stdin {output.tmp} > {log} 2>&1;
#         bedToBigBed -type=bed4+6 -as=../../schemas/bigNarrowPeak.as {output.tmp} {input.genome_size} {output.out} > {log} 2>&1
#         """


# def get_bigfiles(wildcards):
#     bigfiles = {}
#     bigfiles['bigwigs'] = []; bigfiles['bigpeaks'] = []
#
#     # get all the peak files for all replicates or for the replicates combined
#     if 'condition' in samples:
#         for condition in set(samples['condition']):
#             for assembly in set(samples[samples['condition'] == condition]['assembly']):
#                 bigfiles['bigpeaks'].extend(expand(f"{{result_dir}}/{{peak_caller}}/{condition}-{assembly}.bigNarrowPeak", **config))
#     else:
#         for sample in samples.index:
#             bigfiles['bigpeaks'].extend(expand(f"{{result_dir}}/{{peak_caller}}/{sample}-{samples.loc[sample, 'assembly']}.bigNarrowPeak", **config))
#
#     # get all the bigwigs
#     if config.get('combine_replicates', '') == 'merge' and 'condition' in samples:
#         for condition in set(samples['condition']):
#             for assembly in set(samples[samples['condition'] == condition]['assembly']):
#                 bigfiles['bigwigs'].extend(expand(f"{{result_dir}}/{{peak_caller}}/{condition}-{assembly}.bw", **config))
#     else:
#         for sample in samples.index:
#             bigfiles['bigwigs'].extend(expand(f"{{result_dir}}/{{peak_caller}}/{sample}-{samples.loc[sample, 'assembly']}.bw", **config))
#
#     return bigfiles
#
#
# rule trackhub:
#     """
#     Generate a trackhub which has to be hosted on your own machine, but can then be viewed through the UCSC genome
#     browser.
#     """
#     input:
#         unpack(get_bigfiles)
#     output:
#         directory(expand("{result_dir}/trackhub/", **config))
#     log:
#         "log/trackhub.log"
#     benchmark:
#         expand("{benchmark_dir}/trackhub.log", **config)[0]
#     run:
#         import os
#         import re
#         import trackhub
#         from contextlib import redirect_stdout
#
#         with open(str(log), 'w') as f:
#             with redirect_stdout(f):
#                 # start a shared hub
#                 hub = trackhub.Hub(hub=f"{os.path.basename(os.getcwd())} trackhub",
#                                    short_label=f"{os.path.basename(os.getcwd())} trackhub",
#                                    long_label="Automated trackhub generated by the snakemake-workflows tool: \n"
#                                               "https://github.com/vanheeringen-lab/snakemake-workflows",
#                                    email=config.get('email', 'none@provided.com'))
#
#                 # link the genomes file to the hub
#                 genomes_file = trackhub.genomes_file.GenomesFile()
#                 hub.add_genomes_file(genomes_file)
#
#                 for assembly in set(samples['assembly']):
#                     # TODO: add assembly hub support
#                     # now add each assembly to the genomes_file
#                     genome = trackhub.Genome(assembly)
#                     genomes_file.add_genome(genome)
#
#                     # each trackdb is added to the genome
#                     trackdb = trackhub.trackdb.TrackDb()
#                     genome.add_trackdb(trackdb)
#                     priority = 1
#
#                     for peak_caller in config['peak_caller']:
#                         conditions = set()
#                         for sample in samples[samples['assembly'] == assembly].index:
#                             if 'condition' in samples:
#                                 if samples.loc[sample, 'condition'] not in conditions:
#                                     bigpeak = f"{config['result_dir']}/{peak_caller}/{samples.loc[sample, 'condition']}-{assembly}.bigNarrowPeak"
#                                 else:
#                                     bigpeak = False
#                                 conditions.add(samples.loc[sample, 'condition'])
#                                 sample_name = f"{samples.loc[sample, 'condition']}{peak_caller}PEAK"
#                             else:
#                                 bigpeak = f"{config['result_dir']}/{peak_caller}/{sample}-{assembly}.bigNarrowPeak"
#                                 sample_name = f"{sample}{peak_caller}PEAK"
#
#                             if bigpeak:
#                                 track = trackhub.Track(
#                                     name=sample_name,           # track names can't have any spaces or special chars.
#                                     source=bigpeak,             # filename to build this track from
#                                     visibility='dense',         # shows the full signal
#                                     tracktype='bigNarrowPeak',  # required when making a track
#                                     priority=priority
#                                 )
#                                 priority += 1
#                                 trackdb.add_tracks(track)
#
#                             bigwig = f"{config['result_dir']}/{peak_caller}/{sample}-{assembly}.bw"
#                             sample_name = f"{sample}{peak_caller}BW"
#
#                             track = trackhub.Track(
#                                 name=sample_name,    # track names can't have any spaces or special chars.
#                                 source=bigwig,       # filename to build this track from
#                                 visibility='full',   # shows the full signal
#                                 color='0,0,0',       # black
#                                 autoScale='on',      # allow the track to autoscale
#                                 tracktype='bigWig',  # required when making a track
#                                 priority = priority
#                             )
#
#                             # each track is added to the trackdb
#                             trackdb.add_tracks(track)
#                             priority += 1
#
#                 # now finish by storing the result
#                 trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=output[0])
