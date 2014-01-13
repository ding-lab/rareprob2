#!bin/bash

while getopts "i:s:d:n:t:m:x:h:l:" arg
do
	case $arg in
		#r) reference=$OPTARG;;
		i) input=$OPTARG;;
		s) sita=$OPTARG;;
		d) distance=$OPTARG;;
		n) iterationNum=$OPTARG;;
		t) delta=$OPTARG;;
		m) miu=$OPTARG;;
		x) isXExist=$OPTARG;;
		h) proThreshold=$OPTARG;;
		l) referencePanelNum=$OPTARG;;
		?) echo "unknown parameter";;
	esac
done

if [  -z $input ]
then 	
	echo ""
	echo "input directory is necessary!  "
	echo ""
	echo "e.g. -i ./data/input "
	echo ""
	exit;
fi

if [  ! -d $input  ]
then 	
	echo $input "is not a dir"
	exit;
fi

if [ -z $isXExist  ]
then 
	echo ""
	echo "-x is necessary! only 0 and 1 is allowed  "
	echo ""
	echo "e.g. -x 0 "
	echo "e.g. -x 1"
	echo ""
	exit;
fi




dirname=$input
referencePrefix="./data/referencePanel";

for file in `ls "$dirname" `

	do 
		echo ""
		echo "rareprob2: $file is processing......."
		echo ""

#		echo "1. filter begin"
		inputFileName=$file
		input="$dirname""/""$inputFileName"
		reference="$referencePrefix""/""$file"
		java -jar intervalR.jar -r $reference -i $input `if [ $sita ];then echo " -s $sita";fi` `if [ $distance];then echo " -d $distance";fi` `if [ $iterationNum];then echo " -n $iterationNum";fi` `if [ $delta ]; then echo " -t $delta"; fi ` `if [ $miu ]; then echo " -m $miu"; fi` `if [ $isXExist ]; then echo " -x $isXExist"; fi ` `if [ $proThreshold ]; then echo " -h $proThreshold"; fi ` `if [ $referencePanelNum ]; then echo " -l $referencePanelNum"; fi `
		filteredFile="./data/filtered/""${inputFileName}"
		siteInfo="./data/input_siteInfo/"$inputFileName".siteList"
		
#		echo "2. initialR over"
		if [ "$isXExist" == 0 ]
			then 
			./rareprob $filteredFile $siteInfo 
		elif [ "$isXExist" == 1 ]
			then
				./rareprob -x $filteredFile $siteInfo 
		else			
			echo "error usage rareprob2"
			echo ""
		fi
	
	done
echo "finished!"
echo ""

