# qtar
**Quickly Target for smRNA**
## Introduction
qtar is a software for small rna target prediction against to  
UTR and/or coding sequence. It consists of two parts: anchor  
sets and target. The first part is the acceleration step, which  
reduce the prediction counts by only targeting those pairs  
with common anchor between small rna and targets. The  
second part is core step, which contains alignment by using  
Smith-Waterman and filteration by MFE, alignment score  
and bubble size.
## Getting started
	git clone https://github.com/zhuqianhua/qtar.git   
	make   
	cd example   
	sh run.sh   
## Options
Options | Type| Descriptions
--- | --- | ---
-s, --smrna | * str | small rna sequence in fasta format
-t, --target | * str | target sequence in fasta format
-o, --output | str | output, default ./qtar_output.xls
-m, --mode | str | mode of target, a for animal and p for plant
  |   | a same as: -l 6 -b 6 -a 15 -e -12  
  |   | p same as: -S -l 7 -b 2 -a 15 -e -9
-S, --strict | bool | only involve pairing of the miRNA seed
-l, --length | int | length of anchor, default 6
-p, --process | int | number of threads, default 6
-b, --bubble | int | length of bubbles, default 6
-a, --aln | float | cutoff of align score, default 10.0
-e, --mfe | float| cutoff of minimal free energy, default -8.0
-q, --quiet | bool | print nothing except serious errors
-h, --help | bool | print this information
## Output


## Notes
It relies on multithreaded execution and your system needs  
to support -lpthread. If you have questions about qtar, you  
may send the questions to zhuqianhua@bgi.com.

## Reference
1)	Nurul-Syakima Ab Mutalib, Siti Aishah Sulaiman, Rahman Jamal. COMPUTATIONAL TOOLS FOR MICRORNA TARGET PREDICTION. Computational Epigenetics and Diseases. 2019
2)	http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138
3)	 Helwak A, Kudla G, Dudnakova T, Tollervey D. Mapping the Human miRNA Interactome by CLASH Reveals Frequent Noncanonical Binding. Cell. 2013;153:654¨C65. pmid:23622248
4)	miRBase: annotating high confidence microRNAs using deep sequencing data.
Kozomara A, Griffiths-Jones S. NAR 2014 42:D68-D73
5)	Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4:e05005
6)	John B, Enright AJ, Aravin A, Tuschl T, Sander C, Marks DS., PLoS Biol. 2005 Jul;3(7)
7)	Rehmsmeier, Marc and Krueger, Jan RNAhybrid: microRNA target prediction easy, fast and flexible, Nucleic Acids Res.
8)	Xinbin Dai, Zhaohong Zhuang and Patrick X. Zhao (2018). psRNATarget: a plant small RNA target analysis server (2017 release). Nucleic Acids Research. doi: 10.1093/nar/gky316.
9)	Wu HJ, Ma YK, Chen T, Wang M, Wang XJ. (2012) PsRobot: a web-based plant small RNA meta-analysis toolbox. Nucleic Acids Res. DOI:10.1093/nar/gks554
10)	Fahlgren N, Carrington JC (2010) miRNA Target Prediction in Plants. Methods in molecular biology (Clifton, NJ) 592
11)	Freier et al., Improved free-energy parameters for predictions of RNA duplex stability. Proc. Natl. Acad. Sci. USA 83: 9373-9377 (1986))