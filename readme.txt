README file for GENSENG software distribution



1. Compilation
==============

You must have GUN g++ available on your computer.

$ cd src/
$ make

A file called GENSENG will be created after the compilation is completed

2. Running algorithms
=========================

Running involves first setting several parameters outlined below. 

./GENSENG readcountdata round mc1 mc2 map autor mixture tran init human postprocessing.

readcountdata: the window count data generted from the bam file. One example is the reviseddata4hmm_300bpSlides_chr1.txt. We will soon provide the pipeline generating the window count data from the bam file.
round: the em inference rounds [1-]
mc1: the mixture component parameter for the normal state [0-1]
mc2: the mixture component parameter for the CNV state [0-1]
map: whether corrects the mapability bias (1), or not (0)
autor: whether uses the auto regressive component (1), or not (0)
mixture: whether uses the mixture component (1), or not (0)
tran: whether re-estimate the transition probability of HMM (1), or not (0)
init: whether re-estimate the intial probablity of HMM (1), or not (0)
human: the input data is human (1), or mouse (0)
postprocessing: whether the postprocessing step will be executed (1) , or not (0).

For the provided example data, the suggested parameters will be

./GENSENG reviseddata4hmm_300bpSlides_chr1.txt 10 0.01 0.01 1 1 1 0 0 1 1

To run, please also download the file called transition_init.dat, and put it to the same folder of the executable file.

3. Results
========================

There will be two files been generated.

Jingerbread_chr1reviseddata4hmm_300bpSlides_chr1.txtreviseddata4hmm_300bpSlides_chr1.txt_segment.dat (referred as file1) and
Jingerbread_chr1reviseddata4hmm_300bpSlides_chr1.txtreviseddata4hmm_300bpSlides_chr1.txt_SNP.dat (referred as file2)

File1 contains the information for each CNV call. The format is 
chr	start	end	state	cn	sample	score	n	mscore	ave.mprop	ave.gprop	winct.ratio	exp.winct

chr, start, end give the boundary of the CNV call.
state and cn gives the copy number information of the CNV call
sample gives the sample name
score is the confidence score of the CNV
n is number of windows of the CNV has. 
mscore is the average confidence score per window of the CNV
ave.mprop, ave.grpop give the average mapability, the average gc content of the CNV
winct.ratio gives the read depth accessible ratio of the CNV. See our paper for a detail explaination. In our paper, CNV was viewed as detectable using NGS data with winct.ratio < 0.75 for deletion and winct.ratio > 1.25.
exp.winct gives the theoretical read depth of the CNV estimated from the underline copy number, the mappability and the gc content of the CNV region. 

File2 contains the information for each window of the window count data. The format is

chr	str	end	mprop	gprop	wincount	state	stateP	CN 

chr str end give the position of the window
mprop, gprop and wincount give the mapability, the gc content, and the wincount of the particular window
state, CN gives the inferred copy number for this window
stateP gives the confidence score associated with the inferred copy number, for this window 


Please send your comments and suggestions to wangwb09@cs.unc.edu
