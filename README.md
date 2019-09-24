# Human Removal benchmark
Comparing BBMap to Bowtie2 for removal of host contamination in human
microbiome shotgun sequencing data.

## Results

CPU/mem

    bbmap.fecal_200__V300019092_L01_1.txt
    s           h:m:s    max_rss   max_vms    max_uss   max_pss   io_in     io_out    mean_load
    11741.9862  3:15:41  19936.88  107738.37  19926.11  19930.88  25452.21  32106.87  0.00
    bbmap.vag_200__V300019092_L01_2.txt
    s           h:m:s    max_rss   max_vms    max_uss   max_pss   io_in     io_out    mean_load
    29994.1121  8:19:54  21317.20  112711.51  21316.31  21316.37  15782.80  21128.40  0.00
    bowtie2.fecal_200__V300019092_L01_1.txt
    s           h:m:s    max_rss  max_vms  max_uss  max_pss  io_in    io_out     mean_load
    17687.8932  4:54:47  3438.63  5017.18  3427.23  3428.95  9035.55  110415.20  0.00
    bowtie2.vag_200__V300019092_L01_2.txt
    s           h:m:s     max_rss  max_vms  max_uss  max_pss  io_in   io_out    mean_load
    36815.6363  10:13:35  3518.61  5017.18  3510.13  3510.42  603.71  89250.04  0.00

Performance

    kraken2/bbmap.fecal_200__V300019092_L01_1.human.kreport
    94.74  214697  0  F  9604  Hominidae
    kraken2/bbmap.fecal_200__V300019092_L01_1.kreport
    0.11  134871  0  F  9604  Hominidae
    kraken2/bbmap.vag_200__V300019092_L01_2.human.kreport
    99.47  153230411  0  F  9604  Hominidae
    kraken2/bbmap.vag_200__V300019092_L01_2.kreport
    60.96  3461241  0  F  9604  Hominidae
    kraken2/bowtie2.fecal_200__V300019092_L01_1.human.kreport
    98.80  104459  0  F  9604  Hominidae
    kraken2/bowtie2.fecal_200__V300019092_L01_1.kreport
    0.12  143242  0  F  9604  Hominidae
    kraken2/bowtie2.vag_200__V300019092_L01_2.human.kreport
    99.79  75052857  0  F  9604  Hominidae
    kraken2/bowtie2.vag_200__V300019092_L01_2.kreport
    71.48  5356010  0  F  9604  Hominidae

