usage(){
cat << EOF
INPUT: Folder for incoming nanopore reads and output folder for outgoing nanopore results. Threads.
OUTPUT: Aforementioned nanopore results; now processed with pangia!


OPTIONS: PLOT COLUMN:
    -i  Input directory
    -t  THREADS
    -p  PanGIA code directory
    -d PanGIA database directory
    -x Project directory
    -o  Output directory
    -h  Help
EOF
}
while getopts "i:t:p:d:x:o:h" OPTION
do
    case $OPTION in
    i) IN=$OPTARG
           ;;
    t) THREADS=$OPTARG
           ;;
	  p) TOOLLOC=$OPTARG
		       ;;
		d) DBLOC=$OPTARG
		       ;;
    x) PROJDIR=$OPTARG
           ;;
    o) OUT=$OPTARG
           ;;
    *) usage
    exit
           ;;
    esac
done

if [ -d "$OUT" ]; then rm -Rf $OUT; fi
mkdir -p $OUT
mkdir -p $OUT/processing
mkdir -p $OUT/postprocessing
start=0
#write in counter for maximum number of iterations without a pangia run
while true #continuation loop
do #do 1
  if [ $(ls $IN/ 2> /dev/null | wc -l) != 1 ] #if 1
  then
	#if there is a file in this folder #note different operating systems return different values for an empty folder -- 0 or 1
	mv $IN/*fastq $OUT/processing/ #real time should be employed for nanopore data, which is not paired, so no risk of capturing R1 and not R2
  #	grep -vf $OUT/tmp/done.txt $OUT/tmp/nanopore.txt > $OUT/tmp/run.txt #create run file for nanopore reads to run        while read file #execution loop
	for file in $OUT/processing/*fastq #execution loop #this for behavior may change in different OSs
	do #do 2
		prefix=$(basename $file .fastq) #get prefix
		echo $prefix
		prefix=$OUT/tmp/$prefix #improve to include temporary output directory
		echo $prefix
		awk 'BEGIN { cntr = 0 } /^@/ { cntr++ ; print "@read"cntr } !/^@/ { print $0 }' $file > $prefix.renamed.fastq #rename reads just to be sure we don't hit a read name error
    #run minimap2 on each DB
		for db in $DBLOC/*mmi;
		do #do 3
			db_name=$(basename $db .mmi)
			minimap2 -a -x map-ont -t $THREADS $db $prefix.renamed.fastq -o $prefix.$db_name.sam #removed parameters to better foster alignments
			samtools view -hb -T $DBLOC/$db_name $prefix.$db_name.sam > $prefix.$db_name.bam #testing inclusion of fasta headers using samtools view -T db.fa
			if [ $start -eq 0 ] #if 2
			#if it is the first one run, rename to master.bam because we need to have one o' them fer the merger
				then
				cp $prefix.$db_name.bam $OUT/master.bam
				start=1 #so that nothing else gets named master inappropriately
			fi #if 2
			#samtools merge files to master.sam so we can constantly update and rerun pangia on latest sam file.
			#Note - it is important to use the entire samfile every time because of pangia's cutoffs
			samtools merge --threads $THREADS -uf $OUT/tmp.bam $OUT/master.bam $prefix.$db_name.bam #TEST that we can overwrite file we're reading from when using the -f flag. Otherwise scoot with a tmp file
            mv $OUT/tmp.bam $OUT/master.bam #temp file needed because the -f breaks wheen using the same file
		done #done 3

    #pangia on master.sam
		samtools view $OUT/master.bam > $OUT/master.sam #convert back to samfile so pangia can run
        gawk -F\\\\t '!/^@/ { print }' $OUT/master.sam | gawk -F\\\\t '!and($2,256) && !and($2,2048) { print } END { print NR }' | gawk -F\\\\t '!and($2,4) { print }'|grep -v "|Hc" | awk 'NF>=2' > $OUT/master.postawk.sam #format samfile for pangia
        #gawk commands taken out of pangia for samfile processing. grep -v Hc for removing host chromosomes that trip up pangia. awk to remove anything that doesn't have a second field - some artifact of samtools merge
		python $TOOLLOC/pangia.py -st nanopore --kt --database $DBLOC/database/*mmi -s $OUT/master.postawk.sam -o $OUT/master.results -t $THREADS  #consider using "-s -" to pipe the samfile in from stdout of samtools view above. Saves one more disk write/read

		mv $file $OUT/postprocessing/$prefix.fastq #push to completed directory - note $file contains the filepath, tooo

		cp $OUT/master.results/*report.tsv $PROJDIR
	done  #done 2
  fi #if 1
	sleep 5 #sleep to allow files to be made #Not strictly necessary
done #done 1
