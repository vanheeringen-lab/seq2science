rule all:
    input:
        [
         "output/aligned_and_sort/1.txt",
         "output/aligned_and_sort/2.txt",
        ]


checkpoint trimming:
    output:
        "output/trimmed/{sample}.txt"
    shell:
        "touch {output}; sleep 1"


rule align:
    input:
        "output/trimmed/{sample}.txt"
    output:
        pipe("output/aligned/{sample}.txt")
    shell:
        "touch {output}; sleep 1"


rule sort:
    input:
        "output/aligned/{sample}.txt"
    output:
        "output/aligned_and_sort/{sample}.txt"
    shell:
        "touch {output}; sleep 1"
