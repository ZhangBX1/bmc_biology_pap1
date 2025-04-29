"/home/bingxu/bgiscsoftware/dnbc4tools2.1.3/dnbc4tools" rna mkref \
--fasta "/home/bingxu/pap1_scrna/all/index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa" \
--ingtf "/home/bingxu/pap1_scrna/all/index/Arabidopsis_thaliana.TAIR10.60.gtf" \
--species Arabidopsis \
--genomeDir "/home/bingxu/pap1_scrna/all/index/index/" \
--chrM Mt \
--threads 20

nohup "/home/bingxu/bgiscsoftware/dnbc4tools2.1.3/dnbc4tools" rna run \
	--name mutant \
	--cDNAfastq1 "/home/bingxu/pap1_scrna/all/mutant_cDNA_1.fq.gz" \
	--cDNAfastq2 "/home/bingxu/pap1_scrna/all/mutant_cDNA_2.fq.gz" \
	--oligofastq1 "/home/bingxu/pap1_scrna/all/mutant_oligo_1.fq.gz" \
	--oligofastq2 "/home/bingxu/pap1_scrna/all/mutant_oligo_2.fq.gz" \
	--genomeDir "/home/bingxu/bgiscsoftware/index/new/" \
	--outdir "/home/bingxu/scrna_paper/" \
	--threads 20 &


nohup "/home/bingxu/bgiscsoftware/dnbc4tools2.1.3/dnbc4tools" rna run \
--name wt \
--cDNAfastq1 "/home/bingxu/pap1_scrna/all/wt_cDNA_1.fq.gz" \
--cDNAfastq2 "/home/bingxu/pap1_scrna/all/wt_cDNA_2.fq.gz" \
--oligofastq1 "/home/bingxu/pap1_scrna/all/control_oligo_1.fq.gz" \
--oligofastq2 "/home/bingxu/pap1_scrna/all/control_oligo_2.fq.gz" \
--genomeDir "/home/bingxu/bgiscsoftware/index/new/" \
--outdir "/home/bingxu/scrna_paper/" \
--threads 20 &

