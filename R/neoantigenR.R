
	##--------------------------------
	## load required packages
	##--------------------------------


	## check if a package is installed or not
	#' check if a package is installed or not
	#' @param mypkg the name of a R package
	#' @examples
	#' is.installed('seqinr')
	#' @return whether the package is installed, if not, then install it
	is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

	if(is.installed('seqinr')==FALSE){install.packages('seqinr')}
	if(is.installed('msa')==FALSE){install.packages('msa')}
	if(is.installed('GenomicRanges')==FALSE){install.packages('GenomicRanges')}
	if(is.installed('Gviz')==FALSE){install.packages('Gviz')}
	if(is.installed('Biostrings')==FALSE){install.packages('Biostrings')}
	if(is.installed('BSgenome')==FALSE){install.packages('BSgenome')}
	if(is.installed('BSgenome.Hsapiens.UCSC.hg19')==FALSE){install.packages('BSgenome.Hsapiens.UCSC.hg19')}


	library(seqinr)
	library(msa)
	library(GenomicRanges)
	library(Gviz)
	library(Biostrings)
	library(BSgenome)
	library(BSgenome.Hsapiens.UCSC.hg19)


	##--------------------------------
	## required input file names
	##--------------------------------

	#' the path for reference protein database file (Swiss-uniprot)
	#' @examples
	#' protein.database.file.name="swissuniprots.fasta"
	#' @return the path for reference protein database file (Swiss-uniprot)
	protein.database.file.name		   =	"D:\\coding\\bioconductor\\neoantigenR\\neoantigenR\\data\\swissuniprots.fasta"
	#' the path for reference annotated gene database file
	#' @examples
	#' reference.gff.file="gencode.v19.annotation.gff3"
	#' @return the path for reference annotated gene database file
	reference.gff.file					     =	"D:\\coding\\bioconductor\\neoantigenR\\neoantigenR\\data\\gencode.v19.annotation.gff3"
	#' the predicted gff file
	#' @examples
	#' pacbio.gff="cufflinks.gff"
	#' @return the predicted gff file
	pacbio.gff							         =	"D:\\coding\\bioconductor\\neoantigenR\\neoantigenR\\data\\model.gff"
	#' the predicted gff file intersecting with reference gene annotation file by BEDTools
	#' @examples
	#' pacbio.gencode.overlapping.file	 = "bedtool.intersect.overlaps.txt"
	#' @return the predicted gff file intersecting with reference gene annotation file by BEDTools
	pacbio.gencode.overlapping.file	 =	"D:\\coding\\bioconductor\\neoantigenR\\neoantigenR\\data\\bedtool.intersect.overlaps.txt"
	#' the folder of the output files
	#' @examples
	#' output.folder = "analysis"
	#' @return the folder of the output files
	output.folder					         	 =	"D:\\coding\\bioconductor\\neoantigenR\\analysis\\"

	#' the organism used to search protein database
	#' @examples
	#' org="hg19"
	#' @return the organism used to search protein database
	org="hg19"

	##-----------------------------------
	## setup parameters
	##-----------------------------------

	#' indicate whether we will use reference gene annotation or now
	#' @examples
	#' use.reference.annotation.gff=TRUE
	#' @return indicate whether we will use reference gene annotation or now
	use.reference.annotation.gff      = TRUE
	#' indicate whether we will write dna sequence for protein translation
	#' @examples
	#' write.dna.seq.by.reference.genome	=	TRUE
	#' @return indicate whether we will write dna sequence for protein translation
	write.dna.seq.by.reference.genome	=	TRUE 	# if true, name Hsapiens instance
	#' the alternative isoform and reference protein sequence must be 80 percentage or more overlap
	#' @examples
	#' percOverlap	= 0.8
	#' @return the alternative isoform and reference protein sequence must be 80 percentage or more overlap
	percOverlap							=	0.8 	        # used in determing the overlap in protein MSA analysis, disabled !!
	#' mimimum anchor size in the indel region
	#' @examples
	#' minimum.flanking.region.size = 10
	#' @return mimimum anchor size in the indel region
	minimum.flanking.region.size		=	10
	#' at least 3 amino acids must be different from the annotation
	#' @examples
	#' min.mismatch.size=3
	#' @return at least 3 amino acids must be different from the annotation
	min.mismatch.size					=	3  		      # need at least a region >= 3 AAs
	#' at most 1000 amino is different from the annotation
	#' @examples
	#' max.mismatch.size=1000
	#' @return at most 1000 amino is different from the annotation
	max.mismatch.size					=	1000
	#' the alternative isoform and reference protein sequence must be 50 percentage or more overlap
	#' @examples
	#' minimum.sequence.similarity     = 0.5
	#' @return the alternative isoform and reference protein sequence must be 50 percentage or more overlap
	minimum.sequence.similarity     = 0.5   # at least 50 percentage sequence overlap
	#' look for the genomic location of the indel sequence in the predicted gff gene models
	#' @examples
	#' compute.query.exon.details=TRUE
	#' @return look for the genomic location of the indel sequence in the predicted gff gene models
	compute.query.exon.details			=	TRUE
	#' look for the genomic location of the indel sequence in the annotated gff gene models
	#' @examples
	#' compute.database.exon.details=TRUE
	#' @return look for the genomic location of the indel sequence in the annotated gff gene models
	compute.database.exon.details		=	TRUE
	#' whether to validate indels in the annotated protein database
	#' @examples
	#' validate.in.database.boolean=TRUE
	#' @return whether to validate indels in the annotated protein database
	validate.in.database.boolean		=	TRUE

	#' the annotated gene models for chosen genes
	#' @examples
	#' reference.gff.info.chosen=data.frame()
	#' @return the annotated gene models for chosen genes
	reference.gff.info.chosen			=	data.frame()

	#' coverage file to calculate sequence coverage for indel regions with BAM files
	#' @examples
	#' coverage.file= "alignment.bam"
	#' @return coverage file to calculate sequence coverage for indel regions with BAM files
	coverage.file						= ""			# might need to be removed
	#' the predicted gene model's gene ids (prediction names)
	#' @examples
	#' refseq.pacbio.name.ids=c()
	#' @return the predicted gene model's gene ids (prediction names)
	refseq.pacbio.name.ids				= c()  		# might need to be removed


	##-----------------------------------
	## define global variables
	##-----------------------------------

	#' the gff file in the dat frame mode
	#' @examples
	#' gff.file=data.frame()
	#' @return the gff file in the dat frame mode
	gff.file=data.frame()
	#' the gff file's gene names
	#' @examples
	#' gff.names="prediction.gff"
	#' @return the gff file's gene names
	gff.names=""
	#' the directly where code and example data are stored
	#' @examples
	#' mainDir="data"
	#' @return the directly where code and example data are stored
	mainDir=""
	#' the directly where data is stored
	#' @examples
	#' dataDir=""
	#' @return the directly where data is stored
	dataDir=""
	#' the human protein names from uniprot
	#' @examples
	#' uniprot.protein.names.new=""
	#' @return the human protein names from uniprot
	uniprot.protein.names.new=""
	#' the human protein names from uniprot swissport (only the rows or indexes out of the entire uniprot file)
	#' @examples
	#' uniprot.protein.swissprot.index=""
	#' @return the human protein names from uniprot swissport (only the rows or indexes out of the entire uniprot file)
	uniprot.protein.swissprot.index=""
	#' the annotated gene in data frame mode
	#' @examples
	#' reference.gff=data.frame()
	#' @return the annotated gene in data frame mode
	reference.gff=data.frame()
	#' the annotated gene in data frame mode with attributes splitted
	#' @examples
	#' reference.gff.info=data.frame()
	#' @return the annotated gene in data frame mode with attributes splitted
	reference.gff.info=data.frame()
	#' the chosen column of gff file (to remove unused columns from analysis)
	#' @examples
	#' chosen.column.indicies=c()
	#' @return the chosen column of gff file (to remove unused columns from analysis)
	chosen.column.indicies=c()
	#' the uniprot human dataset in data frame
	#' @examples
	#' uniprot.human.database=data.frame()
	#' @return the uniprot human dataset in data frame
	uniprot.human.database=data.frame()
	#' the gene list finally chosen for neoantigen analysis
	#' @examples
	#' target.gene.list=c("MET")
	#' @return the gene list finally chosen for neoantigen analysis
	target.gene.list=c()
	#' the mapping between predicted gene names and reference annotated gene names
	#' @examples
	#' sample.pacbio.gencode.id.mapping=data.frame()
	#' @return the mapping between predicted gene names and reference annotated gene names
	sample.pacbio.gencode.id.mapping=data.frame()
	#' the data frame of predicted gene models
	#' @examples
	#' gene.df.pacbio=data.frame()
	#' @return the data frame of predicted gene models
	gene.df.pacbio=data.frame()
	#' the isoform or gene names from predicted gff model
	#' @examples
	#' pacbio.isof.ids=c()
	#' @return the isoform or gene names from predicted gff model
	pacbio.isof.ids=c()
	#' the attribute information of the gff file
	#' @examples
	#' att.info=c()
	#' @return the attribute information of the gff file
	att.info=c()

	##--------------------------------
	## utility functions
	##--------------------------------


	## use R to split the file into pieces
	#' read the gff gene prediction files
	#' @param gffFile the name of input file
	#' @param nrows the num of rows to skip
	#' @return the data frame of gff input file
	gffRead <- function(gffFile, nrows = -1) {
		 cat("Reading ", gffFile, ": ", sep="")
		 gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
		 header=FALSE, comment.char="#", nrows = nrows,
		 colClasses=c("character", "character", "character", "integer",  "integer", "character", "character", "character", "character"))
		 colnames(gff) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
		 cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, class), collapse=", "), "\n")
		 stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
		 return(gff)
	}


	## use R to map the gene annotation name and prediction name
	#' read the gff gene prediction data frame
	#' @param gff.file the name of a gene prediction data frame
	#' @return a data frame with annotated namme to gene prediction mapping
	retrieve.geneName.pacbioID.mapping.v2<-function(gff.file){
		geneName.info=unlist(lapply(as.vector(gff.file[,18]), function(n){strsplit(n, ";")[[1]][7]}))
		gencode.name=unlist(lapply(geneName.info, function(p){substr(p, 11, nchar(p))}))
		sample.pacbio.gencode.id.mapping=cbind(gff.file$attributes, gencode.name)
		sample.pacbio.gencode.id.mapping=sample.pacbio.gencode.id.mapping[!duplicated(sample.pacbio.gencode.id.mapping),]
		return(sample.pacbio.gencode.id.mapping)
	}

	## write the dna sequence of a predicted gene based on reference genome object
	#' write the dna sequence of a predicted gene based on reference genome object
	#' @param Hsapiens the name of input file
	#' @param gene.df a data frame (gff) of predicted gene's coordinates
	#' @param file.name the output file name to be written
	#' @param gene.name the gene of interest for dna sequence extraction
	#' @return a dna sequence file will be written to output file
	write.dna.sequence.reference<-function(Hsapiens, gene.df, file.name, gene.name){
		sink(file.name, append = TRUE)
		for(transcript in unique(gene.df$pacBio)){
			isof.seq=""
			for(index in which(gene.df$pacBio==transcript)){
				this.seq=getSeq(Hsapiens, gene.df[index, 1], gene.df[index, 2], gene.df[index, 3])
				isof.seq=paste(isof.seq, this.seq, sep="")
			}
			for(strand in c("+", "-")){
				strand.name=if(strand=="+"){"F"}else{"R"}
				for(j in c(0,1,2)){
					if(nchar(isof.seq)>6){
						seqv=seqinr::translate(s2c(isof.seq),sens=strand.name, frame=j)
						stop.count=length(which(seqv=="*"))
						seqv[which(seqv=="*")]="-"  # 10/19/2016, tried to disable it, or change it to "+" doesn't work
						if(stop.count<1000){
							pseq=paste(seqv, collapse="")
							cat(paste(">", gene.name, "_", strand.name, "_", j, sep=""))
							cat("\n")
							cat(pseq)
							cat("\n")
						}
					}
				}
			}
		}
		sink()
	}

	## find the place (exon) of the indel sequence in annotation and predicted genes
	#' find the place (exon) of the indel sequence in annotation and predicted genes
	#' @param chosen.gene the gene of interest
	#' @param chose.pacbio.transcript.id chosen gene's prediction name
	#' @param indel.seq the chosen gene's alternative sequences
	#' @param pacbio.protein.sequences the protein sequence of the predicted gene
	#' @param gene.df.pacbio the data frame of the predicted genes (Gff)
	#' @param reference.gff.info.chosen the data frame of the annnotated genes
	#' @param pacbio.isof.ids the predicted gene's prediction id
	#' @param compute.query.exon.details whether to find the neoantigen's position in predictions
	#' @param compute.database.exon.details whether to find the neoantigen's position in annotations
	#' @return a vector of the predicted gene's genomic locations and locations in annotations
	find.indels.location.by.location<-function(chosen.gene, chose.pacbio.transcript.id, indel.seq,
		  pacbio.protein.sequences, gene.df.pacbio, reference.gff.info.chosen, pacbio.isof.ids,
		  compute.query.exon.details=FALSE,
		  compute.database.exon.details=FALSE){

		dna.target.pos=gregexpr(indel.seq, pacbio.protein.sequences)[[1]][1]*3-2
		this.gene.gff=gene.df.pacbio[which(gene.df.pacbio$feature=="exon" & pacbio.isof.ids==chose.pacbio.transcript.id),]
		reference.gene.gff=reference.gff.info.chosen[which(reference.gff.info.chosen$pacBio == chosen.gene),]
		strand=unique(this.gene.gff[,'strand'])[1]
		target.pos=0
		pacbio.info=c()
		refseq.info=c()

		if(dna.target.pos>0){
			if(compute.query.exon.details==TRUE & nrow(this.gene.gff)>0){
				prev.exon.length.sum=0
				if(strand=="+"){
  				for(i in 1:nrow(this.gene.gff)){
  					sstart=this.gene.gff[i, 'start']
  					send=this.gene.gff[i, 'end']
  					srange=send - sstart + 1 # should add 1 ?
  					if(prev.exon.length.sum<=dna.target.pos & dna.target.pos<=(prev.exon.length.sum+srange)){
  						target.pos=sstart+(dna.target.pos-prev.exon.length.sum)-1
  						pacbio.info=as.character(unlist(this.gene.gff[i,]))
  						break
  					}
  					prev.exon.length.sum=prev.exon.length.sum+srange
  				}
				}else{
				  for(i in nrow(this.gene.gff):1){
				    sstart=this.gene.gff[i, 'start']
				    send=this.gene.gff[i, 'end']
				    srange=send - sstart + 1 # should add 1 ?
				    if(prev.exon.length.sum<=dna.target.pos & dna.target.pos<=(prev.exon.length.sum+srange)){
				      target.pos=send-(dna.target.pos-prev.exon.length.sum)+1
				      pacbio.info=as.character(unlist(this.gene.gff[i,]))
				      break
				    }
				    prev.exon.length.sum=prev.exon.length.sum+srange
				  }
	  		}
			}

			if(compute.database.exon.details==TRUE & nrow(reference.gene.gff)){
				for(i in 1:nrow(reference.gene.gff)){
					sstart=reference.gene.gff[i, 'start']
					send=reference.gene.gff[i, 'end']
					srange=send - sstart +  1 # should I add 1 ?
					if(sstart<=target.pos & target.pos<=send){
						refseq.info=as.character(unlist(reference.gene.gff[i,]))
						break
					}
				}
			}
		}
		return(c(target.pos, pacbio.info, refseq.info))
	}




	##-----------------------------------
	## data preprossing
	##-----------------------------------

	## initialize the parameters and functions for neoantigen prediction
	#' initialize the parameters and functions for neoantigen prediction
	#' @return global variables will be set after initialization
	neoantigenR.initialize<-function(){

    	mainDir		<<-	output.folder
    	#subDir		<<-	"data"
    	dataDir		<<-	file.path(mainDir) #, subDir)
    	#if (!file.exists(subDir)){dir.create(dataDir)}
    	if (!file.exists(dataDir)){dir.create(dataDir)}

    	gene.df.pacbio					=	gffRead(pacbio.gff)
    	pacbio.isof.ids					=	as.vector(gene.df.pacbio[, 'attributes'])
    	pacbio.isof.ids					=	unlist(lapply(pacbio.isof.ids, function(id){id=strsplit(id,  " ")[[1]][4]; substr(id, 2, nchar(id)-2)}))
    	gene.df.pacbio[, 'attributes']	=	pacbio.isof.ids
    	pacbio.isof.ids <<- pacbio.isof.ids
    	gene.df.pacbio <<- gene.df.pacbio

    	gff.file = gffRead(pacbio.gencode.overlapping.file)
    	gff.file = gff.file[which(gff.file[,12]=="exon"),]
    	gff.attributes = as.vector(gff.file[, 18])
    	gff.names = unlist(lapply(gff.attributes, function(g){v=strsplit(g, ";")[[1]][7]; substr(v, 11, nchar(v))}))
    	pacbio.transcript.ids = unlist(lapply(as.vector(gff.file$attributes), function(g){v=strsplit(g, ";")[[1]]; if(length(v)>1){v=v[2]}else{v=v[1]}; v2=strsplit(v, '"')[[1]];if(length(v2)>1){v2=v2[2]}else{v2=v2[1]}; v2}))
    	gff.file$attributes = pacbio.transcript.ids # change "gene_id \"PB.1\"; transcript_id \"PB.1.1\";" ==> "PB.1.1"
    	pacbio.transcript.ids <<- pacbio.transcript.ids
    	gff.file <<- gff.file
    	gff.names <<- gff.names

    	sample.pacbio.gencode.id.mapping <<- retrieve.geneName.pacbioID.mapping.v2(gff.file)
    	target.gene.list<<-unique(gff.names)

    	reference.gff					<<-	gffRead(reference.gff.file)
    	att.info						<<-	data.frame(do.call('rbind', strsplit(as.character(reference.gff[,'attributes']),';',fixed=TRUE)))
    	reference.gff.info				<<-	data.frame(reference.gff, att.info)
    	reference.gff.info				<<-	reference.gff.info[which(reference.gff.info$feature=="exon"),]
    	chosen.column.indicies			<<-	c("seqname", "start", "end", "X7", "score", "strand", "feature", "X3", "X1", "X4", "X14")
    	reference.gff.info.chosen		<<-	reference.gff.info[, chosen.column.indicies]
    	reference.gff.info.chosen[,4]	=	unlist(lapply(as.vector(reference.gff.info.chosen[,4]), function(r){strsplit(r, "=")[[1]][2]}))
    	colnames(reference.gff.info.chosen)	=	c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
    	reference.gff.info.chosen		<<-	data.frame(reference.gff.info.chosen, stringsAsFactors=FALSE)
    	reference.gff.info.chosen[,8]	=	as.character(reference.gff.info.chosen[,8])
    	reference.gff.info.chosen[,9]	=	as.character(reference.gff.info.chosen[,9])
    	reference.gff.info.chosen[,10]	=	as.character(reference.gff.info.chosen[,10])

      uniprot.human.database			<<-	read.fasta(protein.database.file.name, seqtype = "AA" ,  as.string = TRUE)
    	uniprot.protein.names.new		<<-	unlist(lapply(1:length(uniprot.human.database), function(s){name=attributes(uniprot.human.database[[s]])$Annot
    										name.part.pos=regexpr("GN=", name)[1]; pref.name=strsplit(substr(name, name.part.pos, nchar(name)), " ")[[1]][1]; substr(pref.name, 4, nchar(pref.name))}))
    	uniprot.protein.swissprot.index	<<-	unlist(lapply(1:length(uniprot.human.database), function(s){name=attributes(uniprot.human.database[[s]])$name
    										if(substr(name, 1, 2)=="sp"){1}else{0}}))
	}


	##-----------------------------------
	## main function module 1 : isoforms
	##-----------------------------------

	## retrieve the predicted gene's exon coordinates and prepare them for visualization
	#' retrieve the predicted gene's exon coordinates and prepare them for visualization
  #' @return a list of gene models will be produced from the GFF file
	neoantigenR.get.Model<-function(){
  	for(suthee.gene in target.gene.list){

  		eef1a1.positive.genes=gff.file[which(gff.names==suthee.gene),]
  		if(nrow(eef1a1.positive.genes)>2){
  			colnames(eef1a1.positive.genes)=paste("V", 1:dim(eef1a1.positive.genes)[2], sep="")
  			eef1a1.positive.genes=eef1a1.positive.genes[which(eef1a1.positive.genes[,12]=="exon" & eef1a1.positive.genes[,3]=="exon"),]
  			vcf.info=data.frame(do.call('rbind', strsplit(as.character(eef1a1.positive.genes[,18]),';',fixed=TRUE)))
  			eef1a1.positive.genes=data.frame(eef1a1.positive.genes, vcf.info)
  			if(length(unique(eef1a1.positive.genes$V1))>1){
  				eef1a1.positive.genes=eef1a1.positive.genes[which(eef1a1.positive.genes$V1==unique(eef1a1.positive.genes$V1)[1]),]
  			}
  			grtrack.list=list()
  			m=1
  			chosen.column.indicies=c('V1','V4','V5',  'V9', 'V19', 'V7','V3', 'X3','X1','X4','X15')

  			if(dim(eef1a1.positive.genes)[1]>2 & length(intersect(chosen.column.indicies, colnames(eef1a1.positive.genes)))==11){

  				gene.df=eef1a1.positive.genes[,chosen.column.indicies]
  				colnames(gene.df)=c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
  				gene.df=gene.df[!duplicated(gene.df),]
  				gene.df=gene.df[which(gene.df$feature=="exon"),]
  				gene.df$width=0

  				gene.df=gene.df[order(gene.df$pacBio, gene.df$start),]
  				chr.start.end.list=unlist(lapply(1:dim(gene.df)[1], function(n){paste(gene.df[n,1:3], collapse="-")}))
  				gene.df=gene.df[which(duplicated(chr.start.end.list)==FALSE),]
  				if(length(unique(gene.df$pacBio))>10){
  					gene.df=gene.df[which(gene.df$pacBio%in%unique(gene.df$pacBio)[1:10]),]
  					duplication.index=duplicated(gene.df[,c("chromosome","start","end", "pacBio")])
  					gene.df=gene.df[which(duplication.index==FALSE),]
  					gene.df=data.frame(gene.df, stringsAsFactors=FALSE)
  					gene.df$transcript=as.vector(gene.df$transcript)
  				}

  				if(use.reference.annotation.gff==TRUE){
  					reference.gff.info.chosen.gene=reference.gff.info.chosen[which(reference.gff.info.chosen$pacBio==suthee.gene),]
  					gene.df.gff=reference.gff.info.chosen.gene[order(reference.gff.info.chosen.gene$transcript, reference.gff.info.chosen.gene$start),]
  					gene.df=rbind(gene.df, gene.df.gff)
  				}
  				granges.df=makeGRangesFromDataFrame(gene.df)
  				granges.df@seqinfo@genome=org
  				my.gen <- genome(granges.df)
  				my.chr <- as.character(unique(seqnames(granges.df)))
  				my.itrack <- IdeogramTrack(genome = my.gen, chromosome = my.chr)
  				my.gtrack <- GenomeAxisTrack()

    			grtrack.list[[m]]=my.itrack; m=m+1
  				grtrack.list[[m]]=my.gtrack; m=m+1

  				isoform.names=paste(gene.df$pacBio, gene.df$transcript, sep=":")
  				for(isoform in unique(isoform.names)){
  					this.gene.df=gene.df[which(isoform.names==isoform),]	# gene.df[which(gene.df$pacBio==isoform),]
  					granges.df=makeGRangesFromDataFrame(this.gene.df)
  					granges.df@seqinfo@genome=org
  					if(write.dna.seq.by.reference.genome==TRUE){
  						file.name=file.path(dataDir, paste(suthee.gene, ".pacbio.dnaseq.protein.sequence.txt", sep=""))
  						write.dna.sequence.reference(Hsapiens, this.gene.df[, 1:4], file.name, suthee.gene)
  					}
  					my.gen <- genome(granges.df)
  					my.chr <- as.character(unique(seqnames(granges.df)))
  					my.grtrack <- GeneRegionTrack(granges.df, genome = my.gen, chromosome = my.chr, name = isoform) #, transcriptAnnotation = "symbol")
  					my.itrack <- IdeogramTrack(genome = my.gen, chromosome = my.chr)
  					my.gtrack <- GenomeAxisTrack()
  					my.atrack <- AnnotationTrack(granges.df, name = "CpG")
  					grtrack.list[[m]]=my.grtrack
  					m=m+1
  				}

  				transcript.type="PacBio"

  			}
  			if(length(grtrack.list)>0){
  				pdf(file.path(dataDir, paste(suthee.gene, "_", transcript.type, ".isoforms.pdf", sep="")))
  				plotTracks(grtrack.list)
  				dev.off()
  			}
  		}
  	}
	}


	##-----------------------------------
	## main function module 2 : alignment
	##-----------------------------------

    ## compute neoantigens using the novel algorithm after the preparation steps
	  #' compute neoantigens using the novel algorithm after the preparation steps
	  #' @return no values will be returned
		neoantigenR.get.peptides<-function(){
	  alignment.indel.detailed.info.table	=	data.frame()
	  for(chosen.gene in target.gene.list){
  		input.file.name=file.path(dataDir, paste(chosen.gene, ".pacbio.dnaseq.protein.sequence.txt", sep=""))
  		if(!is.na(file.size(input.file.name)) & file.size(input.file.name)>0){
  			database.name=chosen.gene
  			matched.database.index=which(uniprot.protein.names.new==chosen.gene & uniprot.protein.swissprot.index==1) ## only swiss-prot ?
  			if(length(matched.database.index)>0){
  				mySequences.all <- readAAStringSet(input.file.name) #Biostrings::readAAStringSet
  				num.unique.isoforms=round(length(mySequences.all)/6)
  				forward.strand.isoform.idx=c(seq(from=1,to=num.unique.isoforms*6,by=6), seq(from=2,to=num.unique.isoforms*6,by=6),seq(from=3,to=num.unique.isoforms*6,by=6))	#1:(num.unique.isoforms*6)
  				for(n.isof in forward.strand.isoform.idx){
  					ref.sequence.name=names(attributes(mySequences.all[n.isof])$ranges)
  					if(length(grep("_R_", ref.sequence.name))>0 | length(grep("_F_", ref.sequence.name))>0){
  						ref.sequence.name=substr(ref.sequence.name, 1, nchar(ref.sequence.name)-4)
  					}
  					database.seq.isoforms=unlist(lapply(matched.database.index, function(idx){uniprot.human.database[[idx]][1]}))
  					for(ds.idx in 1:length(database.seq.isoforms)){
  						database.isoform.name=attributes(uniprot.human.database[[matched.database.index[ds.idx]]])$Annot
  						database.isoform.uniprotID=strsplit(database.isoform.name, "|",fixed=TRUE)[[1]][2]
  						reference.sequence=toString(mySequences.all[n.isof])
  						database.sequence=database.seq.isoforms[ds.idx]
  						mySequences.new=union(reference.sequence, database.sequence) # database.seq.isoforms)
  						myFirstAlignment=msa::msa(mySequences.new, order="input", type="protein")
  						reference.sequence=gsub("-","z", reference.sequence) # z for stop
  						pacbio.protein.sequences=reference.sequence
  						alignment.full.details.strings=msa::msaConvert(myFirstAlignment, type="bios2mds::align")
  						alignment.full.details.matrix=matrix("", length(alignment.full.details.strings), length(alignment.full.details.strings[[1]]))
  						for(i in 1:length(alignment.full.details.strings)){
  							alignment.full.details.matrix[i,]=alignment.full.details.strings[[i]]
  						}
  						alignment.matrix=alignment.full.details.matrix
  						identity.not.equal.vector=rep(0, ncol(alignment.matrix))
  						identity.equal.matrix=matrix(0, nrow(alignment.matrix)-1, ncol(alignment.matrix))
  						for(i in 1:ncol(alignment.matrix)){
  							if(length(unique(alignment.matrix[1:nrow(alignment.matrix),i]))==1){
  								identity.not.equal.vector[i]=0
  								identity.equal.matrix[1:(nrow(alignment.matrix)-1),i]=1
  							}else if(alignment.matrix[1,i]!="-" & length(unique(alignment.matrix[2:nrow(alignment.matrix),i]))==1 &
  								unique(alignment.matrix[2:nrow(alignment.matrix),i])[1]=="-"){
  								identity.not.equal.vector[i]=1
  							}else if(alignment.matrix[1,i]=="-" & length(unique(alignment.matrix[2:nrow(alignment.matrix),i]))==1 &
  								unique(alignment.matrix[2:nrow(alignment.matrix),i])[1]!="-"){
  								identity.not.equal.vector[i]=-1
  							}else{
  								identity.not.equal.vector[i]=2
  								identity.equal.matrix[which(alignment.matrix[2:nrow(alignment.matrix),i]==alignment.matrix[1,i]),i]=1
  							}
  						}
  						if(length(which(identity.not.equal.vector==0))>minimum.sequence.similarity*(min(c(nchar(reference.sequence), nchar(database.sequence))))){
  							k=1
  							mismatch.indel.regions=list()
  							prev.end=0;curr.start=0
  							cnt=0
  							new.common.seq=TRUE
  							for(j in 1:length(identity.not.equal.vector)){
  								if(identity.not.equal.vector[j]==0){
  									if(new.common.seq==TRUE){ # leading 0
  										cnt=1
  										curr.start=j-1  # a new start
  										new.common.seq=FALSE
  									}else{
  										cnt=cnt+1
  									}
  									if(prev.end>0 & cnt>=minimum.flanking.region.size &
  											(curr.start-prev.end)>=min.mismatch.size & (curr.start-prev.end)<=max.mismatch.size){  # delete the following condition 5/24/2017
  												## & (length(which(identity.not.equal.vector==0))/length(identity.not.equal.vector))>=minimum.sequence.similarity){ # minimum.alignment.match.ratio){
  										this.mismatch.indel.region=c(prev.end, curr.start)
  										left.flanking.seq=paste(alignment.matrix[1, (prev.end-minimum.flanking.region.size):(prev.end-1)], collapse="")
  										right.flanking.seq=paste(alignment.matrix[1, (curr.start+1):(curr.start+minimum.flanking.region.size)], collapse="")
  										mismatch.indel.region.seq=paste(alignment.matrix[1, prev.end:curr.start], collapse="")
  										mismatch.indel.region.seq=gsub("\\-", "", mismatch.indel.region.seq)  # could be all "---" if deletions
  										if(regexpr(mismatch.indel.region.seq, reference.sequence)[1]!=(-1) &
  											regexpr("-", left.flanking.seq)[1]<0 & regexpr("-", right.flanking.seq)[1]<0){
  											mismatch.indel.regions[[k]]=this.mismatch.indel.region;k=k+1
  										}
  										prev.end=0 # reset the prev.end, but keep the counter
  									}
  								}
  								if(identity.not.equal.vector[j]!=0 & cnt>0){
  									new.common.seq=TRUE
  									if(cnt>=minimum.flanking.region.size){  # a valid left flanking sequence, so a begining of mismatch/indel
  										prev.end=j
  										cnt=0     # reset counter
  									}
  								}
  							}
  							if(length(mismatch.indel.regions)>0){
  								for(pr in 1:length(mismatch.indel.regions)){
  									prev.end  =mismatch.indel.regions[[pr]][1]
  									curr.start=mismatch.indel.regions[[pr]][2]
  									indel.region.vals=identity.not.equal.vector[prev.end:curr.start]
  									type.indel.region=
  									if(length(unique(indel.region.vals))==1){
  										if(unique(indel.region.vals)==(-1)){
  											"deletion"
  										}else if(unique(indel.region.vals)==1){
  											"insertion"
  										}else{
  											"mismatch"
  										}
  									}else{
  										"mismatch"
  									}
  									indel.seq=
    									if(type.indel.region=="deletion"){
    									  paste(alignment.matrix[1, c((prev.end-minimum.flanking.region.size):(prev.end-1),(curr.start+1):(curr.start+minimum.flanking.region.size))], collapse="")
    									}else{
    										paste(alignment.matrix[1, prev.end:curr.start], collapse="")
    									}
  									indel.seq.larger=indel.seq  # extend both size for 10 AA for neoantigen prediction
  									if(type.indel.region!="deletion"){
  									  left.end=prev.end-10
  									  right.end=curr.start+10
  									  if(left.end<1){left.end=1}
  									  if(right.end>dim(alignment.matrix)[2]){right.end=dim(alignment.matrix)[2]}
  									  indel.seq.larger=paste(alignment.matrix[1, left.end:right.end], collapse="")
  									}
  									indel.region.seq.overlap.ratio=0
  									reference.indel.region.seq=""
  									if(type.indel.region=="mismatch"){
  										reference.indel.region.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")
  										reference.indel.region.seq=gsub("\\-", "", reference.indel.region.seq)
  										indel.region.seq.overlap.ratio=(1-round(length(which(indel.region.vals==2))/length(indel.region.vals),2))
  									}
  									if(type.indel.region=="deletion"){
  										reference.indel.region.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")
  									}
  									indel.seq.no.hypen=gsub("\\-", "", indel.seq)
  									if(nchar(indel.seq.no.hypen)/nchar(indel.seq)<0.200000){
  										indel.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")  # must be a deletion, any database isoform is good
  									}
  									indel.seq=gsub("\\-", "", indel.seq)
  									match.all.database.isoforms.index=
  									if(nchar(indel.seq)<1460){
  										attributes(regexpr(indel.seq, database.seq.isoforms))$match.length
  									}else{
  										1
  									}
  									if(validate.in.database.boolean==FALSE | (validate.in.database.boolean==TRUE & max(match.all.database.isoforms.index)<0)){
  									  chose.pacbio.transcript.id=as.vector(sample.pacbio.gencode.id.mapping[which(sample.pacbio.gencode.id.mapping[,2]==chosen.gene), 1])
  									  pacbio.transcript.index=ceiling(n.isof/6) # if(n.isof%%6==0){n.isof/6}else{round(n.isof/6)+1}
  									  chose.pacbio.transcript.id=strsplit(chose.pacbio.transcript.id[pacbio.transcript.index], ";")[[1]][1] # need global variable from function above Line-1246
  									  output.file.name=file.path(dataDir, paste(chosen.gene, ".", chose.pacbio.transcript.id, ".F", n.isof, ".Ref", ds.idx, ".pacbio.dnaseq.protein.msa.txt", sep=""))
  										sink(output.file.name)
  										print(myFirstAlignment, show="complete")
  										cat("\n")
  										sink()
  										cat(prev.end, " ", curr.start, " ", indel.seq, "\n")
  										if(type.indel.region=="deletion"){
  											curr.start=prev.end-1
  											prev.end=prev.end-minimum.flanking.region.size
  										}
  										if(type.indel.region=="mismatch"){
  											curr.start=prev.end+nchar(indel.seq)-1
  										}
  										indel.pacbio.info=as.vector(unlist(find.indels.location.by.location(chosen.gene, chose.pacbio.transcript.id, indel.seq, pacbio.protein.sequences,
  											gene.df.pacbio, reference.gff.info.chosen, pacbio.isof.ids, compute.query.exon.details, compute.database.exon.details)))
  										average.coverage=0

  										# prev.end is the location in the alignment matrix, not from its original sequence
  										indel.seq.pos=regexpr(indel.seq, reference.sequence)[1]
  										prev.end.new=indel.seq.pos
  										curr.start.new=prev.end.new + (curr.start-prev.end)

  										indel.seq.detailed.info=c(chosen.gene, ref.sequence.name, database.isoform.name, prev.end.new, curr.start.new, prev.end.new*3-2, curr.start.new*3, indel.seq, average.coverage, n.isof, reference.indel.region.seq, indel.region.seq.overlap.ratio, type.indel.region, database.isoform.uniprotID, indel.seq.larger, indel.pacbio.info)
  										if(length(indel.seq.detailed.info)<36){
  											indel.seq.detailed.info=c(indel.seq.detailed.info, rep("", 36-length(indel.seq.detailed.info)))
  										}

  										indel.seq.detailed.info=t(data.frame(indel.seq.detailed.info))
  										alignment.indel.detailed.info.table=rbind(alignment.indel.detailed.info.table, indel.seq.detailed.info) ## Error in match.names(clabs, names(xi)) :	names do not match previous names

  									}
  								}
  							}
  						}
  					}
  				}
  			}
  		}
	  }
	  if(nrow(alignment.indel.detailed.info.table)>0){
	    colnames(alignment.indel.detailed.info.table)=c("Gene","PacBioSequenceID","UniprotProteinID","ProteinStart","ProteinEnd","DNAStart","DNAEnd","IndelSequence","Coverage","Isof.Index", "reference.indel.region.seq", "indel.region.seq.overlap.ratio", "type", "uniprotID", "IndelSequenceAnchorSeq", "IndepPositioninGenome","Chrom","PacBio","region","PacBioStart","PacbioEnd","dot","strand","dot2","PacbioTranscriptID","ReferenceChrom","Reference_exon_start","Reference_exon_end", "dot3", "geneName", "strand2", "exontype","geneID","exonID","transcriptID","unkonwn")
	  }
	  #print(alignment.indel.detailed.info.table)
	  alignment.indel.detailed.info.table <<- alignment.indel.detailed.info.table
	  #print(alignment.indel.detailed.info.table)
	}



	##-----------------------------------
	## Postprocessing, outputing results
	##-----------------------------------

	## write the predicted neoantigens for final results
	#' write the predicted neoantigens for final results
	#' @return a predicted neoantigen file will be written
	neoantigenR.write<-function(){
  	alignment.indel.detailed.info.table=data.frame(alignment.indel.detailed.info.table)
  	alignment.indel.detailed.info.table=alignment.indel.detailed.info.table[!duplicated(alignment.indel.detailed.info.table),]
  	write.table(alignment.indel.detailed.info.table, file.path(dataDir, "putative_neoantigen_candidates.txt"), sep=";", quote=FALSE, row.names=TRUE, col.names=TRUE)

  	sink(file.path(dataDir, "putative_neoantigen_candidate_sequences.fasta"))
  	for(ik in 1:dim(alignment.indel.detailed.info.table)[1]){
  	  header=paste(">", as.vector(alignment.indel.detailed.info.table[ik, 1]), sep="")
  	  cat(header, "\n")
  	  indel.w.anchor.seq=as.vector(alignment.indel.detailed.info.table[ik, 'IndelSequenceAnchorSeq'])
  	  indel.w.anchor.seq=gsub("-", "", indel.w.anchor.seq)
  	  cat(indel.w.anchor.seq, "\n")
  	}
  	sink()
	}






