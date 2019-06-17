# Creating 'dictionary' file from 'FASTA reference'. It is needed in GATK-tools
java -Xmx2g -jar /home/sebasfranch/TOOLS/PICARD/CreateSequenceDictionary.jar REFERENCE=/home/sebasfranch/TOOLS/ReferenceGenome_hg19/hsapiens.hs37d5.fasta OUTPUT=/home/sebasfranch/TOOLS/ReferenceGenome_hg19/hsapiens.hs37d5.dict


mkdir Targets
mkdir merged_BEDs


# tumorals SPS "AA3542 AA3543 AA3544 AA3545 AA3546 AA3547 AA3548 AA3549 AA3575 AA3576 AA3577 AA3578 AA3579 AA3580"

###################################################
# Compute coverage (including by bases) for each samples


for i in AA3582 AA3583 AA3584 AA3585 AA3586 AA3587 AA3588 AA3589 AA3596 AA3597 AA3598 AA3599 AA3600
do
	echo $i

	java -Xmx7g -jar /home/sebasfranch/TOOLS/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
 	   -R /home/sebasfranch/TOOLS/ReferenceGenome_hg19/hsapiens.hs37d5.fasta \
 	   -o ./Targets/$i.target \
  	  -I /bkp/sebasfranch_backup/DATA/bam_data/FAMCOLON_10.11__bwa/$i.bam \
  	  -L ./HUMAN_hsapiens.hs37d5_SureSelectHumanAllExon.v5_ELID_S04380110.bed \
  	  --countType COUNT_FRAGMENTS \
   	 --outputFormat table \
   	 --omitLocusTable \
 	   --omitIntervalStatistics \
 	   --includeRefNSites \
  	  --minMappingQuality 20 \
  	  --read_filter BadCigar \
  	  --nBins 65533 --start 1 --stop 65534 \
  	  --logging_level ERROR \
	    -ct 10 -ct 15 -ct 20 -ct 30 

done

for i in AA3550 AA3551 AA3552 U726 U729
do
	echo $i

	java -Xmx7g -jar /home/sebasfranch/TOOLS/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
 	   -R /home/sebasfranch/TOOLS/ReferenceGenome_hg19/hsapiens.hs37d5.fasta \
 	   -o ./Targets/$i.target \
  	  -I /bkp/sebasfranch_backup/DATA/bam_data/$i.merged.bqsr.bam \
  	  -L ./HUMAN_hsapiens.hs37d5_SureSelectHumanAllExon.v5_ELID_S04380110.bed \
  	  --countType COUNT_FRAGMENTS \
   	 --outputFormat table \
   	 --omitLocusTable \
 	   --omitIntervalStatistics \
 	   --includeRefNSites \
  	  --minMappingQuality 20 \
  	  --read_filter BadCigar \
  	  --nBins 65533 --start 1 --stop 65534 \
  	  --logging_level ERROR \
	    -ct 10 -ct 15 -ct 20 -ct 30 

done




for i in AA3582 AA3583 AA3584 AA3585 AA3586 AA3587 AA3588 AA3589 AA3596 AA3597 AA3598 AA3599 AA3600
do
	echo $i


	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 5){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C5.bed

	sortBed -i ./merged_BEDs/$i.C5.bed | mergeBed -i - > ./merged_BEDs/$i.C5.merge.bed


	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 10){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C10.bed

	sortBed -i ./merged_BEDs/$i.C10.bed | mergeBed -i - > ./merged_BEDs/$i.C10.merge.bed

	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 20){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C20.bed

	sortBed -i ./merged_BEDs/$i.C20.bed | mergeBed -i - > ./merged_BEDs/$i.C20.merge.bed


done


for i in AA3550 AA3551 AA3552 U726 U729
do
	echo $i


	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 5){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C5.bed

	sortBed -i ./merged_BEDs/$i.C5.bed | mergeBed -i - > ./merged_BEDs/$i.C5.merge.bed


	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 10){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C10.bed

	sortBed -i ./merged_BEDs/$i.C10.bed | mergeBed -i - > ./merged_BEDs/$i.C10.merge.bed

	cat ./Targets/$i.target | awk '{ if($1 != "Locus"){ if($4 >= 20){ split($1,pos,":"); print pos[1] "\t" pos[2]-1 "\t" pos[2] } } }' > ./merged_BEDs/$i.C20.bed

	sortBed -i ./merged_BEDs/$i.C20.bed | mergeBed -i - > ./merged_BEDs/$i.C20.merge.bed


done

###################################################
# Get bases with coverage >= 5/10/20 and convert to bed format for each samples
#


###################################################

mkdir multiIntersect

multiIntersectBed -i ./merged_BEDs/*C5.merge.bed | sortBed -i > ./multiIntersect/multiIntersect.C5_file.bed

multiIntersectBed -i ./merged_BEDs/*C10.merge.bed | sortBed -i > ./multiIntersect/multiIntersect.C10_file.bed

multiIntersectBed -i ./merged_BEDs/*C20.merge.bed | sortBed -i > ./multiIntersect/multiIntersect.C20_file.bed


mkdir merged_2.0

awk '{ if($4 == 14){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5} }' ./multiIntersect/multiIntersect.C5_file.bed > ./merged_2.0/ALL_SAMPLES_C5_merged.bed​ #agafo les subregions que apareix les 14 mostres

awk '{ if($4 == 14){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5} }' ./multiIntersect/multiIntersect.C10_file.bed > ./merged_2.0/ALL_SAMPLES_C10_merged.bed​ #agafo les subregions que apareix les 14 mostres

awk '{ if($4 == 14){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5} }' ./multiIntersect/multiIntersect.C20_file.bed > ./merged_2.0/ALL_SAMPLES_C20_merged.bed​ #agafo les subregions que apareix les 14 mostres









######################################################################################################
###################################                            #######################################
###################################      -- In R program --    #######################################
###################################                            #######################################
######################################################################################################


different_coverages <- c("C5", "C10", "C20")

for (i in 1:length(different_coverages)){

	samples <- 18
	coverage <- different_coverages[i]
	print(paste("Plotting multiIntersect ", coverage, " comparision", sep=""))
	file_input <- paste("./multiIntersect/multiIntersect.", coverage, "_file.bed", sep="")
	BED <- read.table(file_input, header=F, sep="\t")
	inds_project <- c("AA3550","AA3551","AA3552","AA3582","AA3583","AA3584","AA3585","AA3586","AA3587","AA3588","AA3589","AA3596","AA3597","AA3598","AA3599","AA3600","U726","U729")

	colnames(BED)[6:(6+samples-1)] <- inds_project

	total_per_sample <- apply(BED[,6:(6+samples-1)], 2, sum)
	max_sample <- inds_project[which(total_per_sample==max(total_per_sample))]



	percentage_per_sample <- (apply(BED[,6:(6+samples-1)], 2, sum)/dim(BED)[1])*100

	pdf_name <- paste("./Presence_in_", coverage, "_subregions_ALL_samples.pdf", sep="")
	pdf(pdf_name)
	data = percentage_per_sample
	#data = (c(total_per_sample, dim(merged)[1], dim(intersect)[1])/dim(BED)[1])*100
	names = c(colnames(BED)[6:(6+samples-1)])
	main_name <- paste("Presence in ", coverage, " subregions", sep="")
	barplot(data, main=main_name, names=names, ylab="%", ylim=c(0,100),col="white", cex.names=0.7, las=2)    # intervals closed on the left 
	dev.off()

	sample1 <- max_sample
	list_comparisions <- list()
	for (j in 1:length(inds_project)){
		name_to_compare <- inds_project[j]
		list_comparisions[[paste(sample1, name_to_compare, sep=".")]] <- length(which(BED[,sample1]==1 & BED[,name_to_compare]==1))
	}

	comparisions <- do.call(cbind, list_comparisions)
	cell_control <- which(colnames(comparisions)==paste(sample1,sample1,sep="."))
	data <- (comparisions/comparisions[1,cell_control])*100


	pdf_name <- paste("./Sharing_regions_", coverage, "_ALL_samples.pdf", sep="")
	pdf(pdf_name)
	main_name <- paste("Sharing_regions_", coverage, "_ALL_samples", sep="")
	barplot(data, main="Sharing regions in c.10 subregions_FIVE_samples_NORMALmerged", names=colnames(comparisions), ylab="%", ylim=c(0,100),col="grey", cex.names=0.6, las=2)    # intervals closed on the left 
	dev.off()
}









