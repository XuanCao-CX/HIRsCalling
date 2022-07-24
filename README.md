HIRsCalling: identifing the high-interacting regions for malaria genomes.
=====
To systematically identify HIRs at whole genome, we first called high level intra- and inter-chromosomal interacting regions separately. High level intra-chromosomal interacting regions calling consisted of four steps: (1) For each chromosome, removing contacts along the diagonal regions that were located 20% length of each chromosome around diagonal. (2) Keeping the bins in which contact frequency ranked on the top 3% of each contact matrix. (3) Keeping the continuous bins and removing the separated bins on contact matrix. (4) Mapping bins on contact matrix to linear genome and getting high level interacting regions. High level Inter-chromosomal interacting regions calling contained three steps (2), (3) and (4) that described above in high level intra-chromosomal calling. We then considered overlap regions between high level intra-chromosomal interacting regions as HIRs.

HIRs for intra-chromosomal

![image](https://user-images.githubusercontent.com/57889560/113144974-71e0b680-9260-11eb-9071-0f467de6563c.png)

HIRs for inter-chromosomal

![image](https://user-images.githubusercontent.com/57889560/113145112-9e94ce00-9260-11eb-8300-be619ef191d2.png)

-----------------------

HIRsCalling Usage
-------
Usage: HIRsCalling.R: Call Pf High-interacting Regions(HIRs) from merged hic concact matrix.


Options:
        -i CHARACTER, --input=CHARACTER
                Input merged square matrix file. File has header and rownames.

        -o CHARACTER, --out_dir=CHARACTER
                Output dir contain 4 files and 2 hetmap dir.
Files: HIRs.intra.txt, HIRs.inter.txt, Regions.intra.bdg, Regions.inter.bdg.
Heatmap dir: heatmap_intra and heatmap_inter. Each heatmap_dir contain contact heatmap for each chr

        -r CHARACTER, --regions=CHARACTER
                The genome regions corresponding to each row of matirx.Format: first 3 columns: chr,star,end.

        --rm_diag_ratio=CHARACTER
                For each intra matrix, it's the ratio of each chr bins number to remove the diagonal bins.
                Not used for inter chrs. Default: 0.3

        --quan_min_intra=CHARACTER
                For each intra matrix, the minimum quantile ratio of each matrix to define the bins
                as signal bins to next step or noise bins. Default:0.9

        --quan_min_inter=CHARACTER
                For each intra matrix, the minimum quantile ratio of each matrix to define the bins 
                as signal bins to next step or noise bins. Default:0.9

        --continous_bins=CHARACTER
                For each matrix after removing the noise bins, we kept the continous bins that 
                number >= continous_bins

        --inter_bin_apper=CHARACTER
                For call inter HIRs of each chrs, we kept the continous bins that appeared 
                number >= inter_bin_apper

        -h, --help
                Show this help message and exit
--------------

Run HIRsCalling
-----------
Rscript HIRsCalling.R -i data/3D7_Ring.matrix -o results/ -r Pf3D7_genome_regions_10k.bed --rm_diag_ratio 0.2 --quan_min_intra 0.97 --quan_min_inter 0.92 --continous_bins 3 --inter_bin_apper 3
