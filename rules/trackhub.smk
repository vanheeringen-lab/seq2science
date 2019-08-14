rule bedgraphish_to_bedgraph:
    input:
        expand("{result_dir}/genrich/{{sample}}-{{assembly}}.bdgish", **config)
    output:
        bedgraph=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.bedgraph", **config)
    shell:
        """
        splits=$(grep -Pno "([^\/]*)(?=\.bam)" {input})
        splits=($splits)
        lncnt=$(wc -l {input} | echo "$(grep -Po ".*(?=\ )")":)
        splits+=($lncnt)

        counter=1
        for split in "${{splits[@]::${{#splits[@]}}-1}}";
        do
            filename=$(grep -Po "(?<=:).*" <<< $split);
            if [[ $filename =~ {wildcards.sample} ]]; then
                startnr=$(grep -Po ".*(?=\:)" <<< $split);
                endnr=$(grep -Po ".*(?=\:)" <<< ${{splits[counter]}});

                lines="NR>=startnr && NR<=endnr {{ print \$1, \$2, \$3, \$4 }}"
                lines=${{lines/startnr/$((startnr + 2))}}
                lines=${{lines/endnr/$((endnr - 1))}}

                awk "$lines" {input} > {output}
            fi
            ((counter++))
        done
        """

def find_bedgraph(wildcards):
    if wildcards.peak_caller == 'genrich':
        suffix = '.bedgraph'
    elif wildcards.peak_caller == 'macs2':
        suffix = '_treat_pileup.bdg'
    else:
        raise NotImplementedError

    return f"{config['result_dir']}/{{peak_caller}}/{wildcards.sample}-{wildcards.assembly}{suffix}"


rule bedgraph_bigwig:
    input:
        bedgraph=find_bedgraph,
        genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        out=expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}.bw", **config),
        tmp=temp(expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}.bedgraphtmp", **config))
    log:
        "logs/bedgraph_bigwig/{peak_caller}/{sample}-{assembly}.log"
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        awk -v OFS='\\t' '{{print $1, $2, $3, $4}}' {input.bedgraph} | sed '/experimental/d' |
        bedSort /dev/stdin {output.tmp} 2>&1;
        bedGraphToBigWig {output.tmp} {input.genome_size} {output.out}  2>&1
        """

def find_narrowpeak_to_big(wildcards):
    if wildcards.peak_caller == 'genrich':
        suffix = '.narrowPeak'
    elif wildcards.peak_caller == 'macs2':
        suffix = '_peaks.narrowPeak'

    return f"{config['result_dir']}/{{peak_caller}}/{{sample}}-{{assembly}}{suffix}"

rule narrowpeak_bignarrowpeak:
    input:
        narrowpeak= expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}_peaks.narrowPeak", **config),
        genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        out=     expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}.bigNarrowPeak", **config),
        tmp=temp(expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}.tmp.narrowPeak", **config))
    conda:
        "../envs/ucsc.yaml"
    shell:
         "LC_COLLATE=C sort -k1,1 -k2,2n {input.narrowpeak} > {output.tmp}; "
         "bedToBigBed -type=bed4+6 -as=../../bigNarrowPeak.as {output.tmp} {input.genome_size} {output.out} 2>&1"

def get_bigfiles(wildcards):
    bigfiles = {}
    bigfiles['bigwigs'] = []; bigfiles['bigpeaks'] = []

    if 'condition' in samples:
        for condition in set(samples['condition']):
            for assembly in samples[samples['condition'] == condition]['assembly']:
                bigfiles['bigpeaks'].append(f"{config['result_dir']}/genrich/{condition}-{assembly}.bigNarrowPeak")
    else:
        for sample in samples.index:
            bigfiles['bigpeaks'].append(f"{config['result_dir']}/genrich/{sample}-{samples.loc[sample, 'assembly']}.bigNarrowPeak")

    for sample in samples.index:
        bigfiles['bigwigs'].append(f"{config['result_dir']}/genrich/{sample}-{samples.loc[sample, 'assembly']}.bw")

    return bigfiles


rule trackhub:
    input:
        unpack(get_bigfiles)
    output:
        directory(expand("{result_dir}/trackhub/", **config))
    log:
        "log/trackhub.log"
    run:
        import os
        import re
        import trackhub
        from contextlib import redirect_stdout

        with open(str(log), 'w') as f:
            with redirect_stdout(f):
                # start a shared hub
                hub = trackhub.Hub(hub=f"test pipeline",
                                   short_label=f"test pipeline",
                                   long_label="todo: long description",
                                   email=config.get('email', 'none@provided.com'))

                # link the genomes file to the hub
                genomes_file = trackhub.genomes_file.GenomesFile()
                hub.add_genomes_file(genomes_file)

                for assembly in set(samples['assembly']):
                    # TODO: add assembly hub support
                    # now add each assembly to the genomes_file
                    genome = trackhub.Genome(assembly)
                    genomes_file.add_genome(genome)

                    # each trackdb is added to the genome
                    trackdb = trackhub.trackdb.TrackDb()
                    genome.add_trackdb(trackdb)
                    # TODO: Keep the ordering of the samples.tsv
                    for priority, bigpeak in enumerate([f for f in input.bigpeaks if f"{assembly}" in f]):
                        name = trackhub.helpers.sanitize(os.path.basename(bigpeak))
                        # name = os.path.basename(bigpeak)[:os.path.basename(bigpeak).find(assembly)]

                        print(name, bigpeak)
                        track = trackhub.Track(
                            name=name,  # track names can't have any spaces or special chars.
                            source=bigpeak,  # filename to build this track from
                            visibility='dense',  # shows the full signal
                            color='0,0,0',  # black
                            tracktype='bigNarrowPeak',  # required when making a track
                            priority=priority
                        )
                        trackdb.add_tracks(track)

                        samcon, project = re.split('-|/', bigpeak)[-3:-1]
                        bigwigs = [bw for bw in input.bigwigs if project in bw and samcon in bw]

                        for bigwig in bigwigs:
                            print('adding ', bigwig)

                            short_label = name # TODO
                            name = trackhub.helpers.sanitize(os.path.basename(bigwig))

                            track = trackhub.Track(
                                # name=f"{name}bigwig",          # track names can't have any spaces or special chars.
                                name=name,          # track names can't have any spaces or special chars.
                                short_label=short_label, #
                                source=bigwig,      # filename to build this track from
                                visibility='full',  # shows the full signal
                                color='0,0,0',      # black
                                autoScale='on',     # allow the track to autoscale
                                tracktype='bigWig',  # required when making a track
                                priority = priority
                            )

                            # each track is added to the trackdb
                            trackdb.add_tracks(track)

                # now finish by storing the result
                trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=output[0])
