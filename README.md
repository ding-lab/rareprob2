RareProb2
===========
RareProb2 is a novel approach focusing on identifying deleterious variants and regions using the germline variants and somatic mutational events from cancer genome sequencing data. Comprehensive testing indicates RareProb2 achieves superior performance than existing approaches.


Usage
-----
        Version 2.0
        Usage:  bash rareprob2.sh [options]

[options]:

	-i: directory of the genotype file 
	-r: directory of the reference panel file
	-x: the causal variant vector X; 0 means no causal vector X is contained in genotype file. 1 means causal vector X is contained in genotype file.

[Selective Options]:

	-u [mutation_rate]: initial mutation rate for imputation processdure(default: 7.78e-4)
	-d [dt]: genetic distance (default: 1.0)
	-n [Iter_num]: the upper bound of iteration times for the estimation of parameter sita,P1 and miu.( default: 20)
	-s [transition_rate]: the initial transition probability for the imputation(default: 5.0)
	-h [pro_threshold]: the posterior probability threshold for removing the information-poor individual data(defalut:1.0e-6)

Note: before run the rareprob2.sh file, you first need to run "make" in the root directory of Rareprob2 to compile the program. After this, you do not need to 



Install
-------
Clone the RareProb2 repo, and build the `rareprob2` binary:

    git clone https://github.com/ding-lab/rareprob2.git
    cd rareprob2
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv rareprob2 /usr/local/bin/


Requirements/Compatibility
------

        1. R language is necessary
        2. Java1.7 or higher is necessary
        3. Unix-like system is required


Data Format Requirement
------

The RareProb2 can process two kinds genotype file format:

1. File with reference causal rare varients information and reference panel information

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
	
1.2 reference panel format ( necessary ) :
        ------------------       
        1000010000100101
        0010100010000000
        0000001000100000
        0100100001000000
        ------------------
	
        the number of site is required to be equal to the number of variant. The reference panel Each element of this matrix represent the genotype of the reference panel. This matrix is necessary since it is needed in the 
	imputation procedure. In Rareprob2 you do not need to prepare the refererence panel file, since rareprob2 can bootstrap them from cases and controls.				

2. File with reference causal rare varients information and reference panel information

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
	
2.2 reference panel format ( necessary ) :
        ----------------
        1000010000100101
        0010100010000000
        0000001000100000
        0100100001000000
        ----------------
        
        the number of site is required to be equal to the number of variant. The reference panel Each element of this matrix represent the genotype of the reference panel. This matrix is necessary since it is needed in the 
        imputation procedure. In Rareprob2 you do not need to prepare the refererence panel file, since rareprob2 can bootstrap them from cases and controls.				


Contact
-------
Please contact Wenke Wang by wwk881211@gmail.com and Jiayin Wang by jwang@genome.wustl.edu if you have any questions.

