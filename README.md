# gfungi
Assembly annotation process for fungi
## Requirements
### Software
* [Python](https://www.python.org/)
* * Three-party python package
  * [matplotlib](https://matplotlib.org/)
  * [numpy](https://numpy.org/doc/stable/index.html)
* [STAR](https://github.com/alexdobin/STAR)
* [stringtie](https://ccb.jhu.edu/software/stringtie/)
* [gffread](https://github.com/gpertea/gffread)
* [seqclean](https://sourceforge.net/projects/seqclean/files/)
* [PASA](https://github.com/PASApipeline/PASApipeline/)
* [gmst](https://github.com/ConesaLab/SQANTI3)
* [BRAKER2](https://github.com/Gaius-Augustus/BRAKER)
* [Augustus](https://github.com/Gaius-Augustus/Augustus)
* [blast+](https://ftp.ncbi.nih.gov/blast/executables/LATEST/)
* [GeMoMa](http://www.jstacs.de/index.php/GeMoMa)
* [EVM](http://evidencemodeler.github.io/)
### Database
* [UniVec](https://ftp.ncbi.nih.gov/pub/UniVec/)



## Summary of Fungus Annotation Questions
### Long intron length

Common fungal introns are 60-80bp in length. If the fungal intron length reaches 100bp or even more than 200 pb, then you need to look at the intron length of its homologous species. If the intron length of the homologous species is normal, but the intron length of the target sample is abnormal, then it is likely that Denovo has a problem with the prediction. It can be judged by counting the length of introns predicted by Denovo.
Denovo predicts that the intron length of the gene is too long. Solution:

(1)、Identify the species through 18S, go to Augustus to find if there is a trained model of the same species or the same genus.Use the trained model to re-predict.
(2)、Retrain the model. Use gff_filter.py to filter the files of the training model. Retrain the model after filtering.
example:
gff_filter.py trainset.aug.gtf --cds 2 --min_cds 500 --max_cds 4000 --cds_exon 0.5  >trainset.aug_new.gtf
