### pcawg.colour.palette.R #########################################################################
#
# Authors: Jennifer Aguiar & Constance Li (constance.li@oicr.on.ca)
#

### DESCRIPTION ####################################################################################
#
# Return standard PCAWG colour palettes. Case insensitive.
# 
# To see all available schemes, set:			scheme='all', return.scheme = FALSE
# To return all full schemes, set: 			scheme='all', return.scheme = TRUE
# To return specific full schemes, set: 	scheme=<wanted scheme>, return.scheme = TRUE
#	 Note: x will be ignored when scheme='all' OR return.scheme = TRUE 
# To return colours for specific values, 

### ARGUMENTS ######################################################################################
# x 				Chracter vector with terms to be mapped to colours. Ignored if scheme='all' or 
#						return.scheme=TRUE
# scheme 			String specifying desired colour scheme. To see all available schemes, use 
#						scheme='all', returns.scheme=FALSE
# fill.colour 		Unrecognized output will be filled with this colour. Default to 'slategrey'
# return.scheme 	TRUE/FALSE. Set to true to return full specified scheme. Set to false to map
#						x to colours 
# 

### MAIN ###########################################################################################

pcawg.colour.palette <- function(
	x = NULL,
	scheme = NULL,
	fill.colour = 'slategrey',
	return.scheme = FALSE
	) {

	# Define all colours 
	# Coding SNV mutation subtypes & consequences 
	nonsynonymous <- '#698B69'
	synonymous <- '#FFD700'
	stop.gain <- '#8B4789'
	stop.loss <- '#DA70D6'
	indel.frameshift <- '#FF8C00'
	indel.nonframeshift <- '#003366'
	splicing <- '#00CED1'
	# Non-coding SNV mutation subtypes, consequences & gene types
	non.coding <- '#A80015'
	promoter <- '#4C191E'
	enhancer <- '#7F000F'
	operator <- '#A84955'
	silencer <- '#E78A96'
	insulator <- '#FFC1C9'
	lncRNA <- '#331900'
	sncRNA <- '#594027'
	tRNA <- '#A87849'
	rRNA <- '#E7B98A'
	miRNA <- '#FFE0C1'
	utr5.utr3 <- '#1A1A1A'
	intronic <- '#4D4D4D'
	intergenic <- '#7F7F7F'
	telomeres <- '#B3B3B3'
	# Structral variant mutation subtypes
	cna.gain <- '#FF0000'
	cna.loss <- '#0000FF'
	inversion <- '#FFA500'
	transposition <- '#7300E7'
	translocation <- '#458B00'
	# Chromosomes
	chr1 <- '#DE47AB'
	chr2 <- '#72BE97'
	chr3 <- '#F7F797'
	chr4 <- '#7C749B'
	chr5 <- '#E85726'
	chr6 <- '#B395F8'
	chr7 <- '#DC8747'
	chr8 <- '#96D53D'
	chr9 <- '#DC85EE'
	chr10 <- '#7D32B3'
	chr11 <- '#88DB68'
	chr12 <- '#78AAF1'
	chr13 <- '#D9C6CA'
	chr14 <- '#336C80'
	chr15 <- '#F7CA44'
	chr16 <- '#32C7C7'
	chr17 <- '#D4C5F2'
	chr18 <- '#995493'
	chr19 <- '#F88B78'
	chr20 <- '#475ECC'
	chr21 <- '#E0BD8C'
	chr22 <- '#9E2800'
	chrX <- '#F2BBD2'
	chrY <- '#B6EBEA'
	# Sex 
	male <- '#B6EBEA'
	female <- '#F2BBD2'
	# Tumour stage 
	st.one <- '#FFFFFF'
	st.two <- '#FFFF00'
	st.three <- '#FFA500'
	st.four <- '#FF0000'
	st.one.two <- '#FFE4B5'
	st.one.three <- '#EEE8AA'
	st.two.one <- '#FFD700'
	st.two.three <- '#F4A460'
	# TNM Cat
	tnm.zero <- '#FFFFFF'
	tn.one <- '#FFD399'
	tn.two <- '#FFAE45'
	tn.three <- '#B87217'
	tn.four <- '#774607'
	m.one <- '#000000'
	tnm.x <- '#708090'
	# Grade 
	gr.one <- '#FFFFFF'
	gr.two <- '#9CF0FC'
	gr.three <- '#335FE5'
	gr.four <- '#003366'
	gr.well <- '#FFFFFF'
	gr.mod <- '#9CE750'
	gr.poor <- '#00CC00'
	gr.un <- '#005900'
	pr.three.three <- '#FFFFFF'
	pr.three.four <- '#FFFF00'
	pr.three.five <- '#CD2990'
	pr.four.three <- '#FFA500'
	pr.four.four <- '#FF0000'
	pr.four.five <- '#A52A2A'
	pr.five.three <- '#8B008B'
	pr.five.four <- '#0000CD'
	pr.five.five <- '#000000'
	# Primary or Met 
	primary <- '#FFFFFF'
	metastatic <- '#7217A5'
	# Generic 
	other <- '#E5E5E5'
	unknown <- '#708090'
	# Tumour Subtype 
	biliary.adenoca <- '#00CD66'
	bladder.tcc <- '#EEAD0E' 
	bone.benign <- '#F0EE60'
	bone.osteosarc <- '#FFD700'
	softtissue.leiomyo <- '#FFEC8B'
	softtissue.liposarc <- '#CDCB50'
	bone.epith <- '#ADAC44'
	breast.adenoca <- '#CD6090' 
	cervix.scc <- '#79CDCD' 
	cns.medullo <- '#D8BFD8'
	cns.piloastro <- '#B0B0B0' 
	cns.gbm <- '#3D3D3D' 
	cns.gbm.alt <- '#4A4A4A' 
	cns.oligo <- '#787878'
	colorect.adenoca <- '#191970' 
	eso.adenoca <- '#1E90FF'
	head.scc <- '#8B2323'
	kidney.rcc <- '#FF4500' 
	kidney.chrcc <- '#B32F0B' 
	liver.hcc <- '#006400' 
	lung.scc <- '#FDF5E6'
	lung.adenoca <- '#FFFFFF' 
	lymph.bnhl <- '#698B22'
	lymph.cll <- '#F4A35D'
	myeloid.mpn <- '#FFC100' 
	myeloid.aml <- '#CD6600' 
	ovary.adenoca <- '#008B8B' 
	panc.adenoca <- '#7A378B'
	panc.endocrine <- '#E066FF' 
	prost.adenoca <- '#87CEFA'
	skin.melanoma <- '#000000' 
	stomach.adenoca <- '#BFEFFF' 
	thy.adenoca <- '#9370DB'
	uterus.adenoca <- '#FF8C69'
	bone.cart <- '#DDCDCD' 
	breast.lobularca <- '#DDCDCD' 
	breast.lobularca.alt <- '#F095BD'
	breast.dcis <- '#DDCDCD'
	lymph.nos <- '#DDCDCD'
	lymph.nos.alt <- '#698B22'
	myeloid.mds <- '#DDCDCD'
	cervix.adenoca <- '#DDCDCD'

	#-----------------------------------------------------------------------------------------------
	# Some input checking & processing 
	if (class(x) == 'factor') {
		stop('x cannot be a factor: please coerce to character before passing')
	}
	# Some parameters override provided input x. 
	if ((scheme == 'all' || return.scheme) && length(x) != 0) {
		warning('Input x ignored when scheme = \'all\' OR return.scheme = TRUE. Returning all schemes')
	}
	scheme <- tolower(scheme)
	x.input <- tolower(x)
	x.input <- gsub('-', '.', x.input)
	if (return.scheme || scheme == 'all') {
		x.input <- NULL
	}

	colour.schemes <- list(
		coding.snv = list(
			levels = c(
				'nonsynonymous',
				'synonymous',
				'stopgain',
				'stoploss',
				'indel.frameshift',
				'indel.nonframeshift',
				'splicing'
				),
			colours = c(
				nonsynonymous,
				synonymous,
				stop.gain,
				stop.loss,
				indel.frameshift,
				indel.nonframeshift,
				splicing
				)
			),
		noncoding.snv = list(
			levels = c(
				'noncoding',
				'promoter',
				'enhancer',
				'operator',
				'silencer',
				'insulator',
				'lncrna',
				'sncrna',
				'trna',
				'rrna',
				'mirna',
				'utr5.utr3',
				'intronic',
				'intergenic',
				'telomeres'
				),
			colours = c(
				non.coding,
				promoter,
				enhancer,
				operator,
				silencer,
				insulator,
				lncRNA,
				sncRNA,
				tRNA,
				rRNA,
				miRNA,
				utr5.utr3,
				intronic,
				intergenic,
				telomeres
				)
			),
		# ADD TFBS
		structural.variants = list(
			levels = c(
				'cna.gain',
				'cna.loss',
				'inversion',
				'transposition',
				'translocation'
				),
			colours = c(
				cna.gain,
				cna.loss,
				inversion,
				transposition,
				translocation
				)
			),
		chromosomes = list(
			levels	= c(
				'1',
				'2',
				'3',
				'4',
				'5',
				'6',
				'7',
				'8',
				'9',
				'10',
				'11',
				'12',
				'13',
				'14',
				'15',
				'16',
				'17',
				'18',
				'19',
				'20',
				'21',
				'22',
				'x',
				'y'
				),
			colours = c(
				chr1,
				chr2,
				chr3,
				chr4,
				chr5,
				chr6,
				chr7,
				chr8,
				chr9,
				chr10,
				chr11,
				chr12,
				chr13,
				chr14,
				chr15,
				chr16,
				chr17,
				chr18,
				chr19,
				chr20,
				chr21,
				chr22,
				chrX,
				chrY
				)
			),
		sex = list(
			levels = c(
				'male',
				'female'
				),
			colours = c(
				male,
				female
				)
			),
		stage.arabic = list(
			levels = c(
				'1',
				'2',
				'3',
				'4'
				),
			colours = c(
				st.one,
				st.two,
				st.three,
				st.four
				)
			),
		stage.roman = list(
			levels = c(
				'i',
				'i.ii',
				'i.iii',
				'ii',
				'ii.i',
				'ii.iii',
				'iii',
				'iv'
				),
			colours = c(
				st.one,
				st.one.two,
				st.one.three,
				st.two,
				st.two.one,
				st.two.three,
				st.three,
				st.four
				)
			),
		t.category = list(
			levels	= c(
				'0',
				'1',
				'2',
				'3',
				'4',
				'x'
				),
			colours = c(
				tnm.zero,
				tn.one,
				tn.two,
				tn.three,
				tn.four,
				tnm.x
				)
			),
		n.category = list(
			levels	= c(
				'0',
				'1',
				'2',
				'3',
				'4',
				'x'
				),
			colours = c(
				tnm.zero,
				tn.one,
				tn.two,
				tn.three,
				tn.four,
				tnm.x
				)
			),
		m.category = list(
			levels = c(
				'0',
				'1',
				'x'
				),
			colours = c(
				tnm.zero,
				m.one,
				tnm.x
				)
			),
		grade = list(
			levels = c(
				'G1',
				'G2',
				'G3',
				'G4'
				),
			colours = c(
				gr.one,
				gr.two,
				gr.three,
				gr.four
				)
			),
		grade.word = list(
			levels = c(
				'well.differentiated',
				'moderately.differentiated',
				'poorly.differentiated',
				'undifferentiated'
				),
			colours = c(
				gr.well,
				gr.mod,
				gr.poor,
				gr.un
				)
			),
		prostate.grade = list(
			levels	= c(
				'3+3',
				'3+4',
				'3+5',
				'4+3',
				'4+4',
				'4+5',
				'5+3',
				'5+4',
				'5+5'
				),
			colours = c(
				pr.three.three,
				pr.three.four,
				pr.three.five,
				pr.four.three,
				pr.four.four,
				pr.four.five,
				pr.five.three,
				pr.five.four,
				pr.five.five
				)
			),
		primary.met = list(
			levels = c(
				'primary',
				'metastatic'
				),
			colours = c(
				primary,
				metastatic
				)
			),
		tumour.subtype = list(
			levels	= c(
				'biliary.adenoca',
				'bladder.tcc',
				'bone.benign',
				'bone.osteosarc',
				'softtissue.leiomyo', 
				'softtissue.liposarc',
				'bone.epith',
				'breast.adenoca',
				'cervix.scc',
				'cns.medullo',
				'cns.piloastro',
				'cns.gbm',
				'cns.gbm.alt',
				'cns.oligo',
				'colorect.adenoca',
				'eso.adenoca',
				'head.scc',
				'kidney.rcc',
				'kidney.chrcc',
				'liver.hcc',
				'lung.scc',
				'lung.adenoca',
				'lymph.bnhl',
				'lymph.cll',
				'myeloid.mpn',
				'myeloid.aml',
				'ovary.adenoca',
				'panc.adenoca',
				'panc.endocrine',
				'prost.adenoca',
				'skin.melanoma',
				'stomach.adenoca',
				'thy.adenoca',
				'uterus.adenoca',
				'bone.cart',
				'breast.lobularca',
				'breast.lobularca.alt',
				'breast.dcis',
				'lymph.nos',
				'lymph.nos.alt',
				'myeloid.mds',
				'cervix.adenoca'
				),
			colours = c(
				biliary.adenoca,
				bladder.tcc,
				bone.benign,
				bone.osteosarc,
				softtissue.leiomyo,
				softtissue.liposarc,
				bone.epith,
				breast.adenoca,
				cervix.scc,
				cns.medullo,
				cns.piloastro,
				cns.gbm,
				cns.gbm.alt,
				cns.oligo,
				colorect.adenoca,
				eso.adenoca,
				head.scc,
				kidney.rcc,
				kidney.chrcc,
				liver.hcc,
				lung.scc,
				lung.adenoca,
				lymph.bnhl,
				lymph.cll,
				myeloid.mpn,
				myeloid.aml,
				ovary.adenoca,
				panc.adenoca,
				panc.endocrine,
				prost.adenoca,
				skin.melanoma,
				stomach.adenoca,
				thy.adenoca,
				uterus.adenoca,
				bone.cart,
				breast.lobularca,
				breast.lobularca.alt,
				breast.dcis,
				lymph.nos,
				lymph.nos.alt,
				myeloid.mds,
				cervix.adenoca
				)
			),
		organ.system = list(
			levels = c(
				'biliary',
				'bladder',
				'bone.softtissue',
				'breast',
				'cervix',
				'cns',
				'colon.rectum',
				'esophagus',
				'head.neck',
				'kidney',
				'liver',
				'lung',
				'lymphoid',
				'myeloid',
				'ovary',
				'pancreas',
				'prostate',
				'skin',
				'stomach',
				'thyroid',
				'uterus'
				),
			colours = c(
				biliary.adenoca,
				bladder.tcc,
				softtissue.leiomyo,
				breast.adenoca,
				cervix.scc,
				cns.oligo,
				colorect.adenoca,
				eso.adenoca,
				head.scc,
				kidney.rcc,
				liver.hcc,
				lung.scc,
				lymph.bnhl,
				myeloid.aml,
				ovary.adenoca,
				panc.adenoca,
				prost.adenoca,
				skin.melanoma,
				stomach.adenoca,
				thy.adenoca,
				uterus.adenoca
				)
			)
		)
	
	# Error if wanted scheme doesn't match existing schemes
	if (is.null(colour.schemes[[scheme]]) && scheme != 'all'){
		stop('Scheme not found!')
	}
	# Return full specified schemes if return.scheme is TRUE
	if (return.scheme & 'all' == scheme) {
		return(colour.schemes)
		} else if (return.scheme & 'all' != scheme) {
			return(colour.schemes[[scheme]])
		} else if (!return.scheme & 'all' == scheme) {
			return(names(colour.schemes))
		}
	
	# Form output colours 
	matched <- match(x.input, colour.schemes[[scheme]]$levels);
	x.colours <- colour.schemes[[scheme]]$colours[matched]
	names(x.colours) <- colour.schemes[[scheme]]$levels[matched]

	# Deal with unrecognized input by setting to fill colour (slategray by default)
	if (any(is.na(x.colours))) {
		warning('Unrecognized input value for x. Default to fill.colour.')
	}
	x.colours[which(is.na(x.colours))] <- fill.colour;

	return(x.colours)
}

# Copyright (c) 2016 Ontario Institute for Cancer Research 

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
