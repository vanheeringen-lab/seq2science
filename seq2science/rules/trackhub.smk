import os.path
import trackhub
from Bio import SeqIO
from multiprocessing import Pool


rule twobit:
    """
    Generate a 2bit file for each assembly
    """
    input:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.2bit", **config),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.2bit.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.2bit.benchmark.txt", **config)[0]
    conda:
        "../envs/ucsc.yaml"
    shell:
        "faToTwoBit {input} {output} >> {log} 2>&1"


rule gcPercent:
    """
    Generate a gc content track

    source: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/gc5Base/
    """
    input:
        twobit=expand("{genome_dir}/{{assembly}}/{{assembly}}.2bit", **config),
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        gcpcvar=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.gc5Base.wigVarStep.txt.gz", **config)),
        gcpc=expand("{genome_dir}/{{assembly}}/{{assembly}}.gc5Base.bw", **config),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.gc5Base.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.gc5Base.benchmark.txt", **config)[0]
    resources:
        mem_gb=5,
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 {wildcards.assembly} {input.twobit} \
            | gzip > {output.gcpcvar}

        wigToBigWig {output.gcpcvar} {input.sizes} {output.gcpc} >> {log} 2>&1
        """


rule cytoband:
    """
    Generate a cytoband track for each assembly

    source: http://genomewiki.ucsc.edu/index.php/Assembly_Hubs#Cytoband_Track
    """
    input:
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        cytoband_bb=expand("{genome_dir}/{{assembly}}/cytoBandIdeo.bb", **config),
        cytoband_bd=temp(expand("{genome_dir}/{{assembly}}/cytoBandIdeo.bed", **config)),
    params:
        schema=f"{config['rule_dir']}/../schemas/cytoBand.as",
    log:
        expand("{log_dir}/trackhub/{{assembly}}.cytoband.log", **config),
    benchmark:
        expand("{log_dir}/trackhub/{{assembly}}.cytoband.benchmark.txt", **config)[0]
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        cat {input.sizes} | 
        bedSort /dev/stdin /dev/stdout | 
        awk '{{print $1,0,$2,$1,"gneg"}}' > {output.cytoband_bd}

        bedToBigBed -type=bed4 {output.cytoband_bd} -as={params.schema} \
        {input.sizes} {output.cytoband_bb} >> {log} 2>&1
        """


def get_masked_regions(contig):
    masked_regions = ""

    masked_seq = str(contig.seq)
    inMasked = False
    mmEnd = 0

    length = len(masked_seq) - 1
    for x in range(0, length + 1):
        # mark the starting position of a softmasked region
        if masked_seq[x].islower() and inMasked == False:
            mmStart = x + 1
            inMasked = True

        # mark end position of softmasked region (can be end of contig)
        elif (not (masked_seq[x].islower()) or x == length) and inMasked == True:
            mmEnd = x
            inMasked = False

            # store softmasked region in a bed3 (chr, start, end) file
            masked_regions += contig.id + "\t" + str(mmStart) + "\t" + str(mmEnd) + "\n"

    return masked_regions


rule softmask_track_1:
    """
    Generate a track of all softmasked regions

    source: https://github.com/Gaius-Augustus/MakeHub/blob/master/make_hub.py
    """
    input:
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    output:
        mask_unsorted=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}_softmasking_unsorted.bed", **config)),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.softmask1.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.softmask1.benchmark.txt", **config)[0]
    threads: 4
    resources:
        mem_gb=2,
    run:
        with open(str(input.genome), "r") as genome_handle, open(str(output.mask_unsorted), "w+") as bed_handle:
            p = Pool(threads)
            # seqIO.parse returns contig data.
            # Each contig is scanned by get_masked_regions (in parallel by imap_unordered).
            # As soon as a contig is scanned, the output is yielded and written to file.
            for softmasked_regions_per_contig in p.imap_unordered(get_masked_regions, SeqIO.parse(genome_handle, "fasta")):
                bed_handle.write(softmasked_regions_per_contig)



rule softmask_track_2:
    """
    Generate a track of all softmasked regions

    source: https://github.com/Gaius-Augustus/MakeHub/blob/master/make_hub.py
    """
    input:
        mask_unsorted=expand("{genome_dir}/{{assembly}}/{{assembly}}_softmasking_unsorted.bed", **config),
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        maskbed=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}_softmasking.bed", **config)),
        mask=expand("{genome_dir}/{{assembly}}/{{assembly}}_softmasking.bb", **config),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.softmask2.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.softmask2.benchmark.txt", **config)[0]
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        bedSort {input.mask_unsorted} {output.maskbed} >> {log} 2>&1

        bedToBigBed -type=bed3 {output.maskbed} {input.sizes} {output.mask} >> {log} 2>&1
        """


rule trackhub_index:
    """
    Generate a searchable annotation & index for each assembly

    source: https://genome.ucsc.edu/goldenPath/help/hubQuickStartSearch.html
    """
    input:
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config), # TODO: add gtf back to input once checkpoints are fixed
    params:
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
    output:
        genePred=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.gp", **config)),
        genePredbed=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.gp.bed", **config)),
        genePredbigbed=expand("{genome_dir}/{{assembly}}/{{assembly}}.bb", **config),
        info=temp(expand("{genome_dir}/{{assembly}}/info.txt", **config)),
        indexinfo=temp(expand("{genome_dir}/{{assembly}}/indexinfo.txt", **config)),
        ix=expand("{genome_dir}/{{assembly}}/{{assembly}}.ix", **config),
        ixx=expand("{genome_dir}/{{assembly}}/{{assembly}}.ixx", **config),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.index.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.index.benchmark.txt", **config)[0]
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        # generate annotation files
        gtfToGenePred -allErrors -geneNameAsName2 -genePredExt {params.gtf} {output.genePred} -infoOut={output.info} >> {log} 2>&1

        genePredToBed {output.genePred} {output.genePredbed} >> {log} 2>&1

        bedSort {output.genePredbed} {output.genePredbed} >> {log} 2>&1

        bedToBigBed -extraIndex=name {output.genePredbed} {input.sizes} {output.genePredbigbed} >> {log} 2>&1

        # generate searchable indexes (by 2: geneId, 8: proteinID, 9: geneName, 10: transcriptName and 1: transcriptID)
        grep -v "^#" {output.info} | awk '{{print $1, $2, $8, $9, $10, $1}}' > {output.indexinfo}

        ixIxx {output.indexinfo} {output.ix} {output.ixx} >> {log} 2>&1
        """


def get_ucsc_name(assembly):
    """
    Returns as first value (bool) whether or not the assembly was found to be in the
    ucsc genome browser, and as second value the name of the assembly according to ucsc
    "convention".
    """
    # strip custom prefix, if present
    assembly = ori_assembly(assembly)

    # patches are not relevant for which assembly it belongs to
    # (at least not human and mouse)
    assembly_np = [split for split in re.split(r"(.+)(?=\.p\d)", assembly) if split != ""][0].lower()

    # check if the assembly matches a ucsc assembly name
    if assembly_np in ucsc_assemblies:
        return True, ucsc_assemblies[assembly_np][0]

    # else check if it is part of the description
    for ucsc_assembly, desc in ucsc_assemblies.values():
        assemblies = desc[desc.find("(") + 1 : desc.find(")")].split("/")
        assemblies = [val.lower() for val in assemblies]
        if assembly_np in assemblies:
            return True, ucsc_assembly

    # if not found, return the original name
    return False, assembly


def get_defaultPos(sizefile):
    # extract a default position spanning the first scaffold/chromosome in the sizefile.
    with open(sizefile, "r") as file:
        dflt = file.readline().strip("\n").split("\t")
    return dflt[0] + ":0-" + str(min(int(dflt[1]), 100000))


def trackhub_data(wildcards):
    """
    generate a workflow specific dictionary with
    all metadata and files to control the trackhub.

    each samples metadata can contain arguments as can be found here:
    https://daler.github.io/trackhub/autodocs/trackhub.BaseTrack.html#trackhub.BaseTrack


    ATAC-/ChIP-seq dict:

    track_data
    ├── assembly_1
    |   ├── peak_caller_1
    |   |   ├── biological_replicate_1
    |   |   |   ├── biological_replicate_1
    |   |   |   |   └── {filepath, name, visibility, etc.}
    |   |   |   ├── technical_replicate_1a
    |   |   |   |   └── {filepath, name, visibility, etc.}
    |   |   |   └── technical_replicate_1b
    |   |   |       └── {filepath, name, visibility, etc.}
    |   |   |
    |   |   └── biological_replicate_2
    |   |
    |   ├── peak_caller_2
    |   |
    |   └── hubfiles
    |       └── {genome, annotation, etc.}
    |
    └── assembly_2


    Alignment/RNA-seq dict:

    track_data
    ├── assembly_1
    |   ├── aligner
    |   |   ├── technical_replicate_1
    |   |   |   ├── unstranded
    |   |   |   |   └── {filepath, name, visibility, etc.}
    |   |   |   ├── forward
    |   |   |   |   └── {filepath, name, visibility, etc.}
    |   |   |   └── reverse
    |   |   |       └── {filepath, name, visibility, etc.}
    |   |   |
    |   |   └── technical_replicate_2
    |   |
    |   └── hubfiles
    |       └── {genome, annotation, etc.}
    |
    └── assembly_2

    """
    track_data = {}
    for assembly in all_assemblies:
        asmbly = ori_assembly(assembly)  # no custom suffix, if present
        track_data[assembly] = {}
        priority = 10.0

        # check if the trackhub exists on UCSC, or if we need to make an assembly hub
        assembly_hub = not get_ucsc_name(assembly)[0]
        if assembly_hub:
            class RequiredFiles:
                """
                Store filepath similar to Track class objects,
                so the rule trackhub inputFunction can extract files easily.
                """
                def __init__(self, source):
                    self.source = source

            track_data[assembly]["hubfiles"] = {}

            sizes_file = f"{config['genome_dir']}/{assembly}/{assembly}.fa.sizes"
            dflt = get_defaultPos(sizes_file) if os.path.exists(sizes_file) else "chr1:0-10000"  # placeholder
            track_data[assembly]["hubfiles"]["genome"] = trackhub.Assembly(
                genome=asmbly,
                twobit_file=f"{config['genome_dir']}/{assembly}/{assembly}.2bit",  # only file input not named "source"
                organism=asmbly,
                defaultPos=dflt,
                scientificName=asmbly,
                description=asmbly,
            )
            # the sizes file is required above
            track_data[assembly]["hubfiles"]["sizes_file"] = RequiredFiles(sizes_file)

            track_data[assembly]["hubfiles"]["cytobands"] = trackhub.Track(
                name="cytoBandIdeo",
                source=f"{config['genome_dir']}/{assembly}/cytoBandIdeo.bb",
                tracktype="bigBed",
                visibility="dense",
                color="0,0,0",  # black
                priority=10.1,
            )

            # add gtf-dependent file(s) only if the gtf has been found
            if has_annotation(assembly):
                track_data[assembly]["hubfiles"]["annotation"] = trackhub.Track(
                    name="annotation",
                    source=f"{config['genome_dir']}/{assembly}/{assembly}.bb",
                    tracktype="bigBed 12",
                    visibility="pack",
                    color="140,43,69",  # bourgundy
                    priority=10.2,
                    searchIndex="name",
                    searchTrix=f"{assembly}.ix",
                )
                # the ix and ixx files are required above
                track_data[assembly]["hubfiles"]["ix"] = RequiredFiles(
                    f"{config['genome_dir']}/{assembly}/{assembly}.ix"
                )
                track_data[assembly]["hubfiles"]["ixx"] = RequiredFiles(
                    f"{config['genome_dir']}/{assembly}/{assembly}.ixx"
                )

            track_data[assembly]["hubfiles"]["gcPercent"] = trackhub.Track(
                name="gcPercent",
                source=f"{config['genome_dir']}/{assembly}/{assembly}.gc5Base.bw",
                tracktype="bigWig",
                visibility="dense",
                color="59,189,191",  # cyan
                priority=10.3,
            )

            track_data[assembly]["hubfiles"]["RMsoft"] = trackhub.Track(
                name="softmasked",
                source=f"{config['genome_dir']}/{assembly}/{assembly}_softmasking.bb",
                tracktype="bigBed",
                visibility="dense",
                color="128,128,128",  # grey
                priority=10.4,
            )

        # workflow specific data
        if get_workflow() in ["atac_seq", "chip_seq"]:
            for peak_caller in config["peak_caller"]:
                track_data[assembly][peak_caller] = {}

                ftype = get_ftype(peak_caller)
                peak_caller_suffix = "" if len(config["peak_caller"]) == 1 else f"_{peak_caller}"
                for brep in set(breps[breps["assembly"] == asmbly].index):
                    track_data[assembly][peak_caller][brep] = {brep: {}}

                    # the biological replicate
                    priority += 1.0
                    track_data[assembly][peak_caller][brep][brep] = trackhub.Track(
                        name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(brep, brep=True)}{peak_caller_suffix}_pk"),
                        tracktype       = "bigNarrowPeak" if ftype == "narrowPeak" else "bigBed",
                        short_label     = None,  # 17 characters max
                        long_label      = None,
                        subgroups       = {},
                        source          = f"{config['result_dir']}/{peak_caller}/{assembly}-{brep}.big{ftype}",  # filename to build this track from
                        visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
                        color           = "0,0,0",  # black
                        priority        = priority  # change the order this track will appear in
                    )
                    if track_data[assembly][peak_caller][brep][brep].tracktype != "bigNarrowPeak":
                        track_data[assembly][peak_caller][brep][brep].autoScale       = "on",  # allow the track to autoscale
                        track_data[assembly][peak_caller][brep][brep].maxHeightPixels = "100:32:8"

                    # the technical replicate(s) that comprise this biological replicate
                    for trep in treps_from_brep[(brep, asmbly)]:
                        priority += 1.0
                        track_data[assembly][peak_caller][brep][trep] = trackhub.Track(
                            name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(trep)}{peak_caller_suffix}_bw"),
                            tracktype       = "bigWig",  # required when making a track
                            short_label     = None,  # 17 characters max
                            long_label      = None,
                            subgroups       = {},
                            source          = f"{config['result_dir']}/{peak_caller}/{assembly}-{trep}.bw",  # filename to build this track from
                            visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
                            color           = "0,0,0",  # black
                            autoScale       = "on",  # allow the track to autoscale
                            maxHeightPixels = "100:32:8",
                            priority        = priority  # change the order this track will appear in
                        )

        elif get_workflow() in ["alignment", "rna_seq"]:
            # add aligner to mirror the other dict structure
            track_data[assembly][config["aligner"]] = {}
            for trep in treps[treps["assembly"] == asmbly].index:
                track_data[assembly][config["aligner"]][trep] = {}
                for bw in strandedness_to_trackhub(trep):
                    folder = "unstranded" if bw == "" else ("forward" if bw == ".fwd" else "reverse")
                    priority += 1.0
                    track_data[assembly][config["aligner"]][trep][folder] = trackhub.Track(
                        name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(trep)}{bw}"),
                        tracktype       = "bigWig",  # required when making a track
                        short_label     = None,  # 17 characters max
                        long_label      = None,
                        subgroups       = {},
                        source          = f"{config['bigwig_dir']}/{assembly}-{sample}.{config['bam_sorter']}-{config['bam_sort_order']}{bw}.bw",  # filename to build this track from
                        visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
                        color           = "0,0,0",  # black
                        autoScale       = "on",  # allow the track to autoscale
                        maxHeightPixels = "100:32:8",
                        priority        = priority  # change the order this track will appear in
                    )

    return track_data


def get_trackhub_files(wildcards):
    """
    extract all files from the trackhub_data dict
    """
    input_files = []
    track_data = trackhub_data(wildcards)
    for assembly in all_assemblies:
        for div in track_data[assembly]:

            # assembly hub files
            if div == "hubfiles":
                for file in track_data[assembly]["hubfiles"].values():
                    hubfile = file.source if hasattr(file, "source") else file.twobit.source
                    input_files.append(hubfile)
                continue

            # sample & replicate files
            for rep in track_data[assembly][div].values():  # rep   = brep/trep
                for track in rep.values():                  # track = peaks & reads/reads per strand
                        input_files.append(track.source)

    return input_files


rule trackhub:
    """
    Generate a UCSC track hub/assembly hub. 
    
    To view the hub, the output directory must be hosted on an web accessible location, 
    and the location of the hub.txt file given the UCSC genome browser at 
    My Data > Track Hubs > My Hubs
    """
    input:
        get_trackhub_files
    output:
        directory(f"{config['result_dir']}/trackhub"),
    params:
        trackhub_data
    message: explain_rule("trackhub")
    log:
        expand("{log_dir}/trackhub/trackhub.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/trackhub.benchmark.txt", **config)[0]
    run:
        import re
        import sys
        import trackhub

        with open(log[0], "w") as f:
            sys.stderr = sys.stdout = f
            track_data = params[0]

            # start a shared hub
            hub = trackhub.Hub(
                hub=config.get("hubname", "trackhub"),
                short_label=config.get("shortlabel", "trackhub"),  # 17 characters max
                long_label=config.get(
                    "longlabel",
                    "Automated trackhub generated by seq2science: \nhttps://github.com/vanheeringen-lab/seq2scsience",
                ),
                email=config.get("email", "none@provided.com"),
            )

            # link a genomes file to the hub
            genomes_file = trackhub.genomes_file.GenomesFile()
            hub.add_genomes_file(genomes_file)

            for assembly in all_assemblies:
                trackdb = trackhub.trackdb.TrackDb()
                priority = 11.0

                # trackhub/assembly hub
                if "hubfiles" in track_data[assembly]:
                    # create an assembly hub for this genome
                    hubfiles = track_data[assembly]["hubfiles"]
                    genome = hubfiles["genome"]

                    # add the assembly hub tracks to the trackdb
                    trackdb.add_tracks(hubfiles["cytobands"])
                    trackdb.add_tracks(hubfiles["gcPercent"])
                    trackdb.add_tracks(hubfiles["RMsoft"])
                    if "annotation" in hubfiles:
                        trackdb.add_tracks(hubfiles["annotation"])

                        # copy the trix files (not supported by the Trackhub package)
                        dir = os.path.join(str(output), assembly)
                        shell(f"mkdir -p {dir}")
                        basename = f"{config['genome_dir']}/{assembly}/{assembly}"
                        for ext in [".ix", ".ixx"]:
                            file_loc = basename + ext
                            link_loc = os.path.join(dir, assembly + ext)
                            shell(f"ln {file_loc} {link_loc}")
                else:
                    # link this data to an existing trackhub
                    assembly_uscs = get_ucsc_name(assembly)[1]
                    genome = trackhub.Genome(assembly_uscs)
                genomes_file.add_genome(genome)
                genome.add_trackdb(trackdb)

                # add the workflow specific files
                for div in [d for d in track_data[assembly] if d != "hubfiles"]:  # div  = peak_caller/aligner
                    for rep in track_data[assembly][div]:                         # rep  = brep/trep
                        for track in track_data[assembly][div][rep]:              # track= peaks & reads/reads per strand
                            sample_track = track_data[assembly][div][rep][track]
                            trackdb.add_tracks(sample_track)

            # now finish by storing the result
            trackhub.upload.upload_hub(hub=hub, host="localhost", remote_dir=output[0])
