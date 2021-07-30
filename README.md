# gfungi
Assembly annotation process for fungi
## Summary of Fungus Annotation Questions
### Long intron length

Common fungal introns are 60-80bp in length. If the fungal intron length reaches 100bp or even more than 200 pb, then you need to look at the intron length of its homologous species. If the intron length of the homologous species is normal, but the intron length of the target sample is abnormal, then it is likely that Denovo has a problem with the prediction. It can be judged by counting the length of introns predicted by Denovo.
Denovo predicts that the intron length of the gene is too long. Solution:

(1)、Identify the species through 18S, go to Augustus to find if there is a trained model of the same species or the same genus.Use the trained model to re-predict.
(2)、Retrain the model.
