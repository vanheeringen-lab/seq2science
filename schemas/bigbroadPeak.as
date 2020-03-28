table bigBroadPeak
"BED6+3 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
<<<<<<< HEAD
    string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
=======
    string name;	     "Name given to a region (preferably unique). Use . if no name is assigned"
>>>>>>> c97e22a6d10cc8c441a266e40bd750997debb646
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
)
