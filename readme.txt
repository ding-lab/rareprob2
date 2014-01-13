
------- Title and Author -----------

Title: rareprob2
Author: wenke Wang

-------- Version and Contact information ------------

Software Version Number: 2.0
Release Date: 2013.12.30
Contact Information(email): wwk881211@gmail.com

---------- Requirements/Compatibility -------------

1. R language is necessary
2. Java1.7 or higher is necessary
3. Unix-like system is required.

---------- Source Code Description --------------

Rareprob2 contains the following file
1. A probabilistic method for identifying rare variants: 
		C++ codes: probrare.cpp global.cpp global.h
		R codes: r/estimate.r r/initialX.r r/function.r r/P_value.r
		
2. Implementation of hidden Markov models (HMM):
		C codes: backward.c baum.c forward.c hmmrand.c hmmutil.c nrutil.c sequence.c verbi.c nrutil.h hmm.h

3. Filtering, interval and calculation of initial R vector 
		java codes: intervalR.jar
4. shell script file gathering the intervalR.jar and rareprob
		linux shell: rareprob2.sh

----------- Data Format Requirement ---------------

The RareProb2 can process two kinds genotype file format:

1. File with reference causal rare varients information

	1.1 genotype format ( necessary ):
				-----------------       
				0000101001100010
				0000100000100000A
				0000001000100000A
				0000100001000000A
				0000000000000000C
				0000000000000000C
				0000000000100100C
				-----------------
	The first line of the file represents reference causal rare varients information. And character '1' represents this site is a causal rare varient site, while character '0' represents this site is a noncausal rare varient site.
	Rest of lines represents genotype of each individual. And the 'A' and 'C' of the last character of every line represents the individual is a case or control.
	
2. File without reference causal rare varients information

	2.1 genotype format ( necessary ):
				-----------------       
				0000100000100000A
				0000001000100000A
				0000100001000000A
				0000000000000000C
				0000000000000000C
				0000000000100100C
				-----------------
	each of lines represents genotype of each individual. And the 'A' and 'C' of the last character of every line represents the individual is a case or control.
	
3. the site information 
				---------
				ZMYM1
				ZNF642
				ZNF643
				ZP4
				ZYG11B
				AHCTF1
				BSND
				CLCNKA	
				DENND2D
				DNAJC16
				GBP5
				----------
	each line shoule be correponding to the site name in the genotype format. and the length of row shouble also be equal to the number of site. Here we use the gene name to represent the site name; This file is not the necessary. if you have this information, you can put it into the diectory "./data/input_siteInfo" with the same filename to your input filename.

---------- Usage and Options --------------
1.1 Usage:	
	$ bash rareprob2.sh -i genotypeFileDirectory -x 1/0
	-i: directory of the genotype file 
	-r: directory of the reference panel file
	-x: the causal variant vector X; 0 means no causal vector X is contained in genotype file. 1 means causal vector X is contained in genotype file.

1.2 Selective Options
	-u [mutation_rate]: initial mutation rate (default: 7.78e-4)
	-d [dt]: genetic distance (default: 1.0)
	-n [Iter_num]: the upper bound of iteration times for the estimation of parameter sita,P1 and miu.( default: 20)
	-s [transition_rate]: the initial transition probability for sita(default: 5.0)
	-h [pro_threshold]: the posterior probability threshold for removing the information-poor individual data(defalut:1.0e-6)
	-l [reference_panel_number]. The number of the reference panel. The large this value, the more time need to be considered.
	-o [site_info]. The site name corresponding to the site   

1.3 Example

Note: before run the rareprob2.sh file, you first need to run "make" in the root directory of Rareprob2 to compile the program. 

e.g.
	$ bash rareprob2.sh -i ./data/genotypeFileDir -x 0

	after you run this command, there will be information contained in three directories: 

	(1) ./data/referencePanel which contains the reference panel from the cases and control for the calculation of initial R. Each reference panel corresponds to the input genotype file.

	(2) ./data/result which contains only one file "statistic" recording the result of p-value, number of causal variant in the reference X, number of the selected causal variants; As default, we will output the result directly on the shell screen, you can re-direct the output to the "./data/result/statistic" when necessary.

	(3) ./RFile which contains the result vector R file named of inputfile name. Each file corresponds to the inputfile.   
Note: the parameter following -i is a directory containing only your target file, you need to put your file into a directory. the larger the reference panel, the more time are required to training the parameter to get the result. So do the number of sites.
	


	
	

	
	
	
