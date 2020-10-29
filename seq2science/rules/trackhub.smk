import os.path
import trackhub
from Bio import SeqIO
from multiprocessing import Pool
import colorsys
import numpy as np
from seq2science.util import color_picker, color_gradient, hsv_to_ucsc, unique, shorten
from math import ceil


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
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
    output:
        genePred=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.gp", **config)),
        genePredbed=temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.gp.bed", **config)),
        genePredbigbed=expand("{genome_dir}/{{assembly}}/annotation.bedBed", **config),
        info=temp(expand("{genome_dir}/{{assembly}}/info.txt", **config)),
        indexinfo=temp(expand("{genome_dir}/{{assembly}}/indexinfo.txt", **config)),
        ix=expand("{genome_dir}/{{assembly}}/annotation.ix", **config),
        ixx=expand("{genome_dir}/{{assembly}}/annotation.ixx", **config),
    log:
        expand("{log_dir}/trackhub/{{assembly}}.index.log", **config),
    benchmark:
        expand("{benchmark_dir}/trackhub/{{assembly}}.index.benchmark.txt", **config)[0]
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        # generate annotation files
        gtfToGenePred -allErrors -geneNameAsName2 -genePredExt {input.gtf} {output.genePred} -infoOut={output.info} >> {log} 2>&1

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


def get_defaultpos(sizefile):
    # extract a default position spanning the first scaffold/chromosome in the sizefile.
    with open(sizefile, "r") as file:
        dflt = file.readline().strip("\n").split("\t")
    return dflt[0] + ":0-" + str(min(int(dflt[1]), 100000))


def get_colors(asmbly):
    """
    Return a dictionary with colors for each track.

    First picks a color for each main track (biological replicates for ATAC-/ChIP-seq or
    forward stranded technical replicate for alignment/RNA-seq).
    Main track colors can be specified in the samples.tsv, in column 'colors'.

    Then add paler colors for each sub track.
    """
    palletes = {}

    # pick colors for each main track
    main_tracks = unique(breps[breps["assembly"] == asmbly].index)
    if "colors" in breps:
        mc = breps[breps.index.isin(main_tracks)].colors\
            .reset_index().drop_duplicates(breps.index.name).set_index(cols[0])\
            ['colors'].to_list()
    else:
        mc = color_picker(len(main_tracks))

    # create a gradient for each main track color
    for n, rep in enumerate(main_tracks):
        tracks_per_main_track = (1 + len(treps_from_brep[(rep, asmbly)])) # ATAC-/ChIP-seq: 1 for breps + 1 per trep. RNA-seq: 1 for forward/reverse strand
        palletes[rep] = color_gradient(mc[n], tracks_per_main_track)

    return palletes


# def trackhub_data(wildcards):
#     """
#     generate a workflow specific dictionary with
#     all metadata and files to control the trackhub.
#
#     each samples metadata can contain arguments as can be found here:
#     https://daler.github.io/trackhub/autodocs/trackhub.BaseTrack.html#trackhub.BaseTrack
#
#
#     ATAC-/ChIP-seq dict:
#
#     track_data
#     ├── assembly_1
#     |   ├── peak_caller_1
#     |   |   ├── biological_replicate_1
#     |   |   |   ├── biological_replicate_1
#     |   |   |   |   └── {filepath, name, visibility, etc.}
#     |   |   |   ├── technical_replicate_1a
#     |   |   |   |   └── {filepath, name, visibility, etc.}
#     |   |   |   └── technical_replicate_1b
#     |   |   |       └── {filepath, name, visibility, etc.}
#     |   |   |
#     |   |   └── biological_replicate_2
#     |   |
#     |   ├── peak_caller_2
#     |   |
#     |   └── hubfiles
#     |       └── {genome, annotation, etc.}
#     |
#     └── assembly_2
#
#
#     Alignment/RNA-seq dict:
#
#     track_data
#     ├── assembly_1
#     |   ├── aligner
#     |   |   ├── technical_replicate_1
#     |   |   |   ├── unstranded
#     |   |   |   |   └── {filepath, name, visibility, etc.}
#     |   |   |   ├── forward
#     |   |   |   |   └── {filepath, name, visibility, etc.}
#     |   |   |   └── reverse
#     |   |   |       └── {filepath, name, visibility, etc.}
#     |   |   |
#     |   |   └── technical_replicate_2
#     |   |
#     |   └── hubfiles
#     |       └── {genome, annotation, etc.}
#     |
#     └── assembly_2
#
#     """
#     track_data = {}
#     for assembly in all_assemblies:
#         asmbly = ori_assembly(assembly)  # no custom suffix, if present
#         track_data[assembly] = {}
#         palletes = get_colors(asmbly)
#         priority = 10.0
#
#         # check if the trackhub exists on UCSC, or if we need to make an assembly hub
#         assembly_hub = not get_ucsc_name(assembly)[0]
#         if assembly_hub:
#             class RequiredFiles:
#                 """
#                 Store filepath similar to Track class objects,
#                 so the rule trackhub inputFunction can extract files easily.
#                 """
#                 def __init__(self, source):
#                     self.source = source
#
#             track_data[assembly]["hubfiles"] = {}
#
#             sizes_file = f"{config['genome_dir']}/{assembly}/{assembly}.fa.sizes"
#             dflt = get_defaultpos(sizes_file) if os.path.exists(sizes_file) else "chr1:0-10000"  # placeholder
#             track_data[assembly]["hubfiles"]["genome"] = trackhub.Assembly(
#                 genome=asmbly,
#                 twobit_file=f"{config['genome_dir']}/{assembly}/{assembly}.2bit",  # only file input not named "source"
#                 organism=asmbly,
#                 defaultPos=dflt,
#                 scientificName=asmbly,
#                 description=asmbly,
#             )
#             # the sizes file is required above
#             track_data[assembly]["hubfiles"]["sizes_file"] = RequiredFiles(sizes_file)
#
#             track_data[assembly]["hubfiles"]["cytobands"] = trackhub.Track(
#                 name="cytoBandIdeo",
#                 source=f"{config['genome_dir']}/{assembly}/cytoBandIdeo.bb",
#                 tracktype="bigBed",
#                 visibility="dense",
#                 color="0,0,0",  # black
#                 priority=10.1,
#             )
#
#             track_data[assembly]["hubfiles"]["gcPercent"] = trackhub.Track(
#                 name="gcPercent",
#                 source=f"{config['genome_dir']}/{assembly}/{assembly}.gc5Base.bw",
#                 tracktype="bigWig",
#                 visibility="dense",
#                 color="59,189,191",  # cyan
#                 priority=10.2,
#             )
#
#             track_data[assembly]["hubfiles"]["RMsoft"] = trackhub.Track(
#                 name="softmasked",
#                 source=f"{config['genome_dir']}/{assembly}/{assembly}_softmasking.bb",
#                 tracktype="bigBed",
#                 visibility="dense",
#                 color="128,128,128",  # grey
#                 priority=10.3,
#             )
#
#             # add gtf-dependent file(s) only if the gtf has been found
#             if has_annotation(assembly):
#                 track_data[assembly]["hubfiles"]["annotation"] = trackhub.Track(
#                     name="annotation",
#                     source=f"{config['genome_dir']}/{assembly}/annotation.bedBed",
#                     tracktype="bigBed 12",
#                     visibility="pack",
#                     color="140,43,69",  # bourgundy
#                     priority=10.4,
#                     searchIndex="name",
#                     searchTrix="annotation.ix",
#                 )
#                 # the ix and ixx files are required above
#                 track_data[assembly]["hubfiles"]["ix"] = RequiredFiles(
#                     f"{config['genome_dir']}/{assembly}/annotation.ix"
#                 )
#                 track_data[assembly]["hubfiles"]["ixx"] = RequiredFiles(
#                     f"{config['genome_dir']}/{assembly}/annotation.ixx"
#                 )
#
#         # workflow specific data
#         if get_workflow() in ["atac_seq", "chip_seq"]:
#             for peak_caller in config["peak_caller"]:
#                 track_data[assembly][peak_caller] = {}
#
#                 ftype = get_ftype(peak_caller)
#                 peak_caller_suffix = "" if len(config["peak_caller"]) == 1 else f" {peak_caller}"
#                 pcs = shorten(peak_caller_suffix, 5)
#                 for brep in unique(breps[breps["assembly"] == asmbly].index):
#                     track_data[assembly][peak_caller][brep] = {"files": []}
#
#                     descriptive = rep_to_descriptive(brep, brep=True)
#                     label = f"{descriptive}{peak_caller_suffix}"
#                     composite = trackhub.CompositeTrack(
#                         name=trackhub.helpers.sanitize(label),
#                         short_label = f"{shorten(descriptive, 17-len(pcs), ['signs', 'vowels', 'center'])}{pcs}",
#                         long_label  = label,
#                         dimensions=f"dimX=peaks dimY=signal",
#                         # dimensions=f"dimX={safe_brep}pk dimY={safe_brep}bw  dimA={peak_caller}",
#                         # filterComposite="dimA",
#
#                         # The availalbe options here are the `name` attributes of each subgroup.
#                         # sortOrder='num=+ kind=-',
#                         #track_type_override="compositeTrack",
#                         tracktype="bigWig",
#                         # tracktype="compositeTrack",
#                         visibility="dense",
#                     )
#
#                     view = trackhub.ViewTrack(
#                         name=trackhub.helpers.sanitize(f"{label} peak files"),
#                         short_label=f"{shorten(descriptive, 13-len(pcs), ['signs', 'vowels', 'center'])}{pcs} pks",  # <= 17 characters suggested
#                         long_label=f"{label} peak files",
#                         view="peaks",
#                         visibility="full",
#                         tracktype="bigNarrowPeak" if ftype == "narrowPeak" else "bigBed",
#                     )
#                     composite.add_view(view)
#
#                     file = f"{config['result_dir']}/{peak_caller}/{assembly}-{brep}.big{ftype}"
#                     track_data[assembly][peak_caller][brep]["files"].append(file)
#                     priority += 1.0
#                     track = trackhub.Track(
#                         name            = trackhub.helpers.sanitize(f"{label} pk"),
#                         tracktype       = "bigNarrowPeak" if ftype == "narrowPeak" else "bigBed",
#                         short_label     = f"{shorten(descriptive, 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} pk",  # <= 17 characters suggested
#                         long_label      = f"{label} peaks",
#                         subgroups       = {},
#                         source          = file,  # filename to build this track from
#                         visibility      = "full",  # full/squish/pack/dense/hide visibility of the track
#                         color           = hsv_to_ucsc(palletes[brep][0]),
#                         priority        = priority  # change the order this track will appear in
#                     )
#                     if track.tracktype != "bigNarrowPeak":
#                         track.autoScale       = "on",  # allow the track to autoscale
#                         track.maxHeightPixels = "100:32:8"
#                     view.add_tracks(track)
#
#
#
#                     signal = trackhub.ViewTrack(
#                         name=trackhub.helpers.sanitize(f"{label} signals"),
#                         short_label=f"{shorten(descriptive, 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} sg",  # <= 17 characters suggested
#                         long_label=f"{label} signals",
#                         view="signal",
#                         visibility="full",
#                         tracktype="bigWig",
#                     )
#                     composite.add_view(signal)
#
#                     # the technical replicate(s) that comprise this biological replicate
#                     for n, trep in enumerate(treps_from_brep[(brep, asmbly)]):
#                         # folder = f"{trep}_bw"  # prevents overwriting if brep == trep
#                         file = f"{config['result_dir']}/{peak_caller}/{assembly}-{trep}.bw"
#                         track_data[assembly][peak_caller][brep]["files"].append(file)
#                         priority += 1.0
#                         track = trackhub.Track(
#                             name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(trep)}{peak_caller_suffix} bw"),
#                             tracktype       = "bigWig",  # required when making a track
#                             short_label     = f"{shorten(rep_to_descriptive(trep), 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} bw",  # <= 17 characters suggested
#                             long_label      = f"{rep_to_descriptive(trep)}{peak_caller_suffix} bigWig",
#                             subgroups       = {},
#                             source          = file,  # filename to build this track from
#                             visibility      = "full",  # full/squish/pack/dense/hide visibility of the track
#                             color           = hsv_to_ucsc(palletes[brep][n+1]),
#                             autoScale       = "on",  # allow the track to autoscale
#                             maxHeightPixels = "100:32:8",
#                             priority        = priority  # change the order this track will appear in
#                         )
#                         signal.add_tracks(track)
#
#                     # print(composite)
#                     # print()
#                     # print(dir(composite))
#                     # # print(composite)
#                     # # print(composite.source)
#                     # print(composite.children)
#                     # # print(composite.filename)
#                     # print(composite.leaves)
#                     # # print(composite.subgroups)
#                     # # print(composite.subtracks)
#                     # print(composite.views)
#                     # # print(composite.validate())
#                     # for child in composite.children:
#                     #     print(dir(child.children))
#                     #     print(child.children)
#                     # exit(0)
#                     track_data[assembly][peak_caller][brep]["composite"] = composite
#
#
#
#
#                     # # the biological replicate
#                     # priority += 1.0
#                     # track_data[assembly][peak_caller][brep][brep] = trackhub.Track(
#                     #     name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(brep, brep=True)}{peak_caller_suffix} pk"),
#                     #     tracktype       = "bigNarrowPeak" if ftype == "narrowPeak" else "bigBed",
#                     #     short_label     = f"{shorten(rep_to_descriptive(brep, brep=True), 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} pk",  # <= 17 characters suggested
#                     #     long_label      = f"{rep_to_descriptive(brep, brep=True)}{peak_caller_suffix} peaks",
#                     #     subgroups       = {},
#                     #     source          = f"{config['result_dir']}/{peak_caller}/{assembly}-{brep}.big{ftype}",  # filename to build this track from
#                     #     visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
#                     #     color           = hsv_to_ucsc(palletes[brep][0]),
#                     #     priority        = priority  # change the order this track will appear in
#                     # )
#                     # if track_data[assembly][peak_caller][brep][brep].tracktype != "bigNarrowPeak":
#                     #     track_data[assembly][peak_caller][brep][brep].autoScale       = "on",  # allow the track to autoscale
#                     #     track_data[assembly][peak_caller][brep][brep].maxHeightPixels = "100:32:8"
#
#                     # # the technical replicate(s) that comprise this biological replicate
#                     # for n, trep in enumerate(treps_from_brep[(brep, asmbly)]):
#                     #     folder = f"{trep}_bw"  # prevents overwriting if brep == trep
#                     #     priority += 1.0
#                     #     track_data[assembly][peak_caller][brep][folder] = trackhub.Track(
#                     #         name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(trep)}{peak_caller_suffix} bw"),
#                     #         tracktype       = "bigWig",  # required when making a track
#                     #         short_label     = f"{shorten(rep_to_descriptive(trep), 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} bw",  # <= 17 characters suggested
#                     #         long_label      = f"{rep_to_descriptive(trep)}{peak_caller_suffix} bigWig",
#                     #         subgroups       = {},
#                     #         source          = f"{config['result_dir']}/{peak_caller}/{assembly}-{trep}.bw",  # filename to build this track from
#                     #         visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
#                     #         color           = hsv_to_ucsc(palletes[brep][n+1]),
#                     #         autoScale       = "on",  # allow the track to autoscale
#                     #         maxHeightPixels = "100:32:8",
#                     #         priority        = priority  # change the order this track will appear in
#                     #     )
#
#         elif get_workflow() in ["alignment", "rna_seq"]:
#             # add aligner to mirror the other dict structure
#             track_data[assembly][config["aligner"]] = {}
#             for trep in treps[treps["assembly"] == asmbly].index:
#                 track_data[assembly][config["aligner"]][trep] = {}
#
#                 for n, bw in enumerate(strandedness_to_trackhub(trep)):
#                     folder = "unstranded" if bw == "" else ("forward" if bw == ".fwd" else "reverse")
#                     priority += 1.0
#                     track_data[assembly][config["aligner"]][trep][folder] = trackhub.Track(
#                         name            = trackhub.helpers.sanitize(f"{rep_to_descriptive(trep)}{bw}"),
#                         tracktype       = "bigWig",  # required when making a track
#                         short_label     = shorten(rep_to_descriptive(trep), 14 if bw else 17, ['signs', 'vowels', 'center']) +
#                                           ("" if bw == "" else (" fw" if bw == ".fwd" else " rv")),  # <= 17 characters suggested
#                         long_label      = rep_to_descriptive(trep) + ("" if bw == "" else (" forward" if bw == ".fwd" else " reverse")),
#                         subgroups       = {},
#                         source          = f"{config['bigwig_dir']}/{assembly}-{sample}.{config['bam_sorter']}-{config['bam_sort_order']}{bw}.bw",  # filename to build this track from
#                         visibility      = "dense",  # full/squish/pack/dense/hide visibility of the track
#                         color           = hsv_to_ucsc(palletes[trep][n]),
#                         autoScale       = "on",  # allow the track to autoscale
#                         maxHeightPixels = "100:32:8",
#                         priority        = priority  # change the order this track will appear in
#                     )
#
#     return track_data
#
#
# def get_trackhub_files(wildcards):
#     """
#     extract all files from the trackhub_data dict
#     """
#     input_files = []
#     track_data = trackhub_data(wildcards)
#     for assembly in all_assemblies:
#         for div in track_data[assembly]:
#
#             # assembly hub files
#             if div == "hubfiles":
#                 for file in track_data[assembly]["hubfiles"].values():
#                     hubfile = file.source if hasattr(file, "source") else file.twobit.source
#                     input_files.append(hubfile)
#                 continue
#
#             # sample & replicate files
#             for rep in track_data[assembly][div].values():  # rep   = brep/trep
#                 input_files.extend(rep["files"])
#                 # for track in rep.values():                  # track = peaks & reads/reads per strand
#                 #         input_files.append(track.source)
#
#     return input_files


def trackhub_data(wildcards):
    """
    Create a UCSC trackhub.
    Returns the hub, and a list of required files.

    Birds eye view of the trackhub system:
    hub
    └──genomes_file
       └──genome (1 per assembly)
          ├──trackdb
          |  ├── peaks, bigwigs
          |  └── annotation files (assembly hub only)
          └──groups_file (assembly hub only)

    trackdb
    ├──(annotation group annotation assembly hub only)
    ├──(gcPercent  group annotation assembly hub only)
    ├──(softmasked group annotation assembly hub only)
    |
    └──composite (1 per peak caller/aligner)
       ├──view: peaks (ChIP- & ATAC-seq only)
       ├──view: signal
       ├──----brep1       peak   (group trackhub - assembly hub only)
       ├──----brep1 trep1 signal (group trackhub - assembly hub only)
       └──----brep1 trep2 signal (group trackhub - assembly hub only)
    """
    out = {
        "hub": None,  # the complete trackhub
        "files": []   # the files required (for the snakemake input)
    }

    # universal hub file
    wf = get_workflow().replace('_', '-')
    hub = trackhub.Hub(
        hub=config.get("hubname", f"{wf}_trackhub"),
        short_label=config.get("shortlabel", f"seq2science {wf} hub"),  # title of the control box in the genome browser (for trackhubs)
        long_label=config.get(
            "longlabel",
            f"Automated {wf} trackhub generated by seq2science: \nhttps://github.com/vanheeringen-lab/seq2scsience",
        ),
        email=config.get("email", "none@provided.com"),
    )
    out["hub"] = hub

    # universal genomes file
    genomes_file = trackhub.genomes_file.GenomesFile()
    hub.add_genomes_file(genomes_file)

    for assembly in all_assemblies:
        asmbly = ori_assembly(assembly)  # no custom suffix, if present
        hub_type = "trackhub" if get_ucsc_name(assembly)[0] else "assembly_hub"
        trackdb = trackhub.trackdb.TrackDb()
        priority = 10.0

        if hub_type == "trackhub":
            # link this data to the existing trackhub
            assembly_uscs = get_ucsc_name(assembly)[1]
            genome = trackhub.Genome(assembly_uscs)

        elif hub_type == "assembly_hub":
            # add the files for an assembly hub

            # the genome
            file = f"{config['genome_dir']}/{assembly}/{assembly}.2bit"
            sizes_file = f"{config['genome_dir']}/{assembly}/{assembly}.fa.sizes"
            dflt = get_defaultpos(sizes_file) if os.path.exists(sizes_file) else "chr1:0-10000"  # placeholder
            genome = trackhub.Assembly(
                genome=asmbly,
                twobit_file=file,
                organism=asmbly,
                defaultPos=dflt,
                scientificName=asmbly,
                description=asmbly,
            )
            out["files"].append(file)
            out["files"].append(sizes_file)

            # groups_file
            annotation_group = trackhub.groups.GroupDefinition(
                name="annotation",
                label="Annotation tracks",
                priority=10,
                default_is_closed=False
            )
            hub_group = trackhub.groups.GroupDefinition(
                name="trackhub",
                label=config.get("shortlabel", f"seq2science {wf} hub"),
                priority=11,
                default_is_closed=False
            )
            groups_file = trackhub.groups.GroupsFile([annotation_group, hub_group])
            genome.add_groups(groups_file)

            # annotation tracks
            file = f"{config['genome_dir']}/{assembly}/cytoBandIdeo.bb"
            track = trackhub.Track(
                name="cytoBandIdeo",
                tracktype="bigBed",
                short_label="cytoBandIdeo",
                long_label="Chromosome ideogram with cytogenetic bands",
                source=file,
            )
            trackdb.add_tracks(track)
            out["files"].append(file)

            file = f"{config['genome_dir']}/{assembly}/{assembly}.gc5Base.bw"
            track = trackhub.Track(
                name="gcPercent",
                source=file,
                tracktype="bigWig",
                visibility="dense",
                color="59,189,191",  # cyan
                group="annotation",
                priority=10.2,
            )
            trackdb.add_tracks(track)
            out["files"].append(file)

            file = f"{config['genome_dir']}/{assembly}/{assembly}_softmasking.bb"
            track = trackhub.Track(
                name="softmasked",
                source=file,
                tracktype="bigBed",
                visibility="dense",
                color="128,128,128",  # grey
                group="annotation",
                priority=10.3,
            )
            trackdb.add_tracks(track)
            out["files"].append(file)

            # add gtf-dependent file(s) only if the gtf has been found
            if has_annotation(asmbly):
                file = f"{config['genome_dir']}/{assembly}/annotation.bedBed"
                track = trackhub.Track(
                    name="annotation",
                    source=f"{config['genome_dir']}/{assembly}/annotation.bedBed",
                    tracktype="bigBed 12",
                    visibility="pack",
                    color="140,43,69",  # bourgundy
                    group="annotation",
                    priority=10.4,
                    searchIndex="name",
                    searchTrix="annotation.ix",
                )
                trackdb.add_tracks(track)
                out["files"].append(file)
                out["files"].append(f"{config['genome_dir']}/{assembly}/annotation.ix")
                out["files"].append(f"{config['genome_dir']}/{assembly}/annotation.ixx")

        genomes_file.add_genome(genome)
        genome.add_trackdb(trackdb)

        # workflow specific data
        if get_workflow() in ["atac_seq", "chip_seq"]:
            for peak_caller in config["peak_caller"]:
                ftype = get_ftype(peak_caller)
                ttype = "bigNarrowPeak" if ftype == "narrowPeak" else "bigBed"
                peak_caller_suffix = "" if len(config["peak_caller"]) == 1 else f" {peak_caller}"
                pcs = shorten(peak_caller_suffix, 5)
                peak_caller_prefix = "" if len(config["peak_caller"]) == 1 else f"{peak_caller} "
                pcp = shorten(peak_caller_prefix, 5)

                # one composite track to rule them all
                name = f"{wf}{peak_caller_suffix} samples"
                safename = trackhub.helpers.sanitize(name)
                composite = trackhub.CompositeTrack(
                    name=safename,
                    short_label=f"{wf}{pcs}",
                    long_label=name,
                    dimensions=f"dimX=view dimY=conditions",
                    tracktype="bigWig",
                    visibility="dense",
                )
                if hub_type == "assembly_hub":
                    composite.add_params(group="trackhub")
                trackdb.add_tracks(composite)

                # two views to guide them
                peaks_view_name = f"{peak_caller_prefix}peaks"
                peaks_view = trackhub.ViewTrack(
                    name=trackhub.helpers.sanitize(peaks_view_name),
                    short_label=f"{pcp}peaks",  # <= 17 characters suggested
                    long_label=peaks_view_name,
                    view="peaks",               # default
                    visibility="dense",         # peaks aren't nice beyond dense
                    tracktype=ttype,
                )
                composite.add_view(peaks_view)

                signal_view_name = f"{peak_caller_prefix}signal"
                signal_view = trackhub.ViewTrack(
                        name=trackhub.helpers.sanitize(signal_view_name),
                        short_label=f"{pcp}signal",    # <= 17 characters suggested
                        long_label=signal_view_name,
                        view="signal",
                        visibility="full",             # default
                        maxHeightPixels = "100:50:8",  # with 50 the y-axis is visible
                        # TODO: vertical limits
                        #autoScale="on",  # on/off
                        tracktype="bigWig",
                    )
                composite.add_view(signal_view)

                # subgroup containing all breps
                subgroup = trackhub.SubGroupDefinition(
                    name='conditions',
                    label='Conditions',
                    mapping={}
                )
                composite.add_subgroups([subgroup])

                for brep in unique(breps[breps["assembly"] == asmbly].index):
                    descriptive = rep_to_descriptive(brep, brep=True)
                    safedescr = trackhub.helpers.sanitize(descriptive)
                    subgroup.mapping[safedescr] = descriptive

                    file = f"{config['result_dir']}/{peak_caller}/{assembly}-{brep}.big{ftype}"
                    priority += 1.0
                    track = trackhub.Track(
                        name            = trackhub.helpers.sanitize(f"{descriptive}{peak_caller_suffix} pk"),
                        tracktype       = ttype,
                        short_label     = f"{shorten(descriptive, 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} pk",  # <= 17 characters suggested
                        long_label      = f"{descriptive}{peak_caller_suffix} peaks",
                        subgroups       = {"conditions": safedescr},
                        source          = file,
                        visibility      = "full",  # full/squish/pack/dense/hide visibility of the track
                        # color           = hsv_to_ucsc(palletes[brep][0]),
                        priority        = priority  # change the order this track will appear in
                    )
                    if ttype != "bigNarrowPeak":
                        track.autoScale       = "on",  # allow the track to autoscale
                        track.maxHeightPixels = "100:50:8"  # with 50 the y-axis is shown
                    out["files"].append(file)
                    peaks_view.add_tracks(track)

                    # the technical replicate(s) that comprise this biological replicate
                    for n, trep in enumerate(treps_from_brep[(brep, asmbly)]):
                        file = f"{config['result_dir']}/{peak_caller}/{assembly}-{trep}.bw"
                        priority += 1.0
                        track = trackhub.Track(
                            name            = trackhub.helpers.sanitize(f"{descriptive}{peak_caller_suffix} bw"),
                            tracktype       = "bigWig",
                            short_label     = f"{shorten(rep_to_descriptive(trep), 14-len(pcs), ['signs', 'vowels', 'center'])}{pcs} bw",  # <= 17 characters suggested
                            long_label      = f"{rep_to_descriptive(trep)}{peak_caller_suffix} bigWig",
                            subgroups       = {"conditions": safedescr},
                            source          = file,
                            visibility      = "full",  # full/squish/pack/dense/hide visibility of the track
                            # color           = hsv_to_ucsc(palletes[brep][n+1]),
                            # autoScale       = "on",
                            # maxHeightPixels = "100:50:8",
                            priority        = priority  # change the order this track will appear in
                        )
                        out["files"].append(file)
                        signal_view.add_tracks(track)

                        # TODO: remove after testing
                        if priority >= 15:
                            return out



    return out


# test the code
# TODO: this needs to be in the actual rule
hub = trackhub_data("")["hub"]
output = [f"{config['result_dir']}/trackhub"]
assembly = asmbly = "XENTR_10.0"

# upload the trix files (not supported by the Trackhub package)
shell(f"mkdir -p {output[0]}/{asmbly}")
for ext in ["ix", "ixx"]:
    src = f"{config['genome_dir']}/{assembly}/annotation.{ext}"
    dst = f"{output[0]}/{asmbly}/annotation.{ext}"
    shell(f"rsync {src} {dst}")

# upload the hub
trackhub.upload.upload_hub(hub=hub, host="localhost", remote_dir=output[0])
exit(0)


def get_trackhub_files(wildcards):
    """
    extract all files from the trackhub_data dict
    """
    return trackhub_data(wildcards)["files"]


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
        print(get_trackhub_files)
        exit(0)

        # import re
        # import sys
        # # import trackhub
        #
        # with open(log[0], "w") as f:
        #     #sys.stderr = sys.stdout = f
        #     track_data = params[0]
        #
        #     # start a shared hub
        #     hub = trackhub.Hub(
        #         hub=config.get("hubname", f"{get_workflow()} trackhub"),
        #         short_label=config.get("shortlabel", f"seq2science {get_workflow()}"),  # title of the control box in the genome browser (for trackhubs)
        #         long_label=config.get(
        #             "longlabel",
        #             f"Automated {get_workflow()} trackhub generated by seq2science: \nhttps://github.com/vanheeringen-lab/seq2scsience",
        #         ),
        #         email=config.get("email", "none@provided.com"),
        #     )
        #
        #     # link a genomes file to the hub
        #     genomes_file = trackhub.genomes_file.GenomesFile()
        #     hub.add_genomes_file(genomes_file)
        #
        #     for assembly in all_assemblies:
        #         trackdb = trackhub.trackdb.TrackDb()
        #         if "hubfiles" in track_data[assembly]:
        #             # create an assembly hub for this genome
        #             hubfiles = track_data[assembly]["hubfiles"]
        #             genome = hubfiles["genome"]
        #
        #             # add the assembly hub tracks to the trackdb
        #             trackdb.add_tracks(hubfiles["cytobands"])
        #             trackdb.add_tracks(hubfiles["gcPercent"])
        #             trackdb.add_tracks(hubfiles["RMsoft"])
        #             if "annotation" in hubfiles:
        #                 trackdb.add_tracks(hubfiles["annotation"])
        #
        #                 # copy the trix files (not supported by the Trackhub package)
        #                 dir = os.path.join(str(output), assembly)
        #                 shell(f"mkdir -p {dir}")
        #                 basename = f"{config['genome_dir']}/{assembly}/{assembly}"
        #                 for ext in ["ix", "ixx"]:
        #                     file_loc = f"{basename}.{ext}"
        #                     link_loc = os.path.join(dir, f"annotation.{ext}")
        #                     shell(f"cp {file_loc} {link_loc}")
        #         else:
        #             # link this data to an existing trackhub
        #             assembly_uscs = get_ucsc_name(assembly)[1]
        #             genome = trackhub.Genome(assembly_uscs)
        #         genomes_file.add_genome(genome)
        #         genome.add_trackdb(trackdb)
        #
        #         # add the workflow specific files
        #         for div in [d for d in track_data[assembly] if d != "hubfiles"]:  # div  = peak_caller/aligner
        #             for rep in track_data[assembly][div].values():                         # rep  = brep/trep
        #                 trackdb.add_tracks(rep["composite"])
        #                 #print(rep["composite"])
        #                 #exit(0)
        #             # for rep in track_data[assembly][div]:                         # rep  = brep/trep
        #             #     for track in track_data[assembly][div][rep]:              # track= peaks & reads/reads per strand
        #             #         sample_track = track_data[assembly][div][rep][track]
        #             #         trackdb.add_tracks(sample_track)
        #
        #     # now finish by storing the result
        #     trackhub.upload.upload_hub(hub=hub, host="localhost", remote_dir=output[0])
