##running R script
library(crossword)
output_folder="crossword_output_01"
iterations=2
input_folder=system.file("extdata/peanut_input",package="crossword")
im_type="svg"
by_chromosomes=TRUE
heterozygous=FALSE
gff="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/peanut.gff"
chr_stat="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/chr_siz.txt"
chr_length="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/chr_len.txt"
window_size=200000
input="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/parental_genotypes.hapmap"
outcross=0.1
gen2phy="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/parental_genotypes.loc"
homo=TRUE
phenotyping_method="QTN_random"
qtn=20
h2=0.6
heritability_mode="average"
effect_distribution="equal"
min_qtn_freq=0
highest_P=NA
lowest_P=NA
high_to_low_percentage=0
input_effects="/home/korani/R/x86_64-pc-linux-gnu-library/3.4/crossword/extdata/peanut_input/qtn_effects.efc"
biased_selection=TRUE
dominant=FALSE
I="start"
I="start"
dir.create(output_folder)
parental_genotypes = get_parental_genotypes(input,gen2phy,homo)
gen2phy = physical2genomic(gff,chr_stat,,window_size)
###############################################
###iteration# 1################################
l=1
dir.create(paste0(output_folder,"/iteration_",l))
qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution)
vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode="average")
loci = parental_genotypes[[2]][as.character(qtn_effect$QTN),]
write.table(cbind(qtn_effect,loci),paste0(output_folder,"/iteration_",l,"/selected_QTNs"),col.names=TRUE,row.names=FALSE,quote =FALSE)
pop="pop1"
pop1=crossword::cross(P1="grg",P2=c("tr","tg"),N=10,id=pop,chr_length=chr_length,parental_genotypes=parental_genotypes,isdh=FALSE)
pop="pop2"
pop2=advance(pop=pop1,F=5,level="individual",clevel="individual",id=pop,chr_length=chr_length)
pop="pop3"
pop3=create_families(pop=pop2,S=5,id=pop,chr_length=chr_length)
pheno = phenotype(pop=pop3,qtn_effect=qtn_effect,tgv_only=FALSE,vr=vr,parental_genotypes=parental_genotypes)
percentage=FALSE
pop6=select(pop=pop3,N=3,level="individual",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_1/pop3")
pop9=combine_populations(c(pop1,pop2))
haplotypes_out(pop=pop2,output="crossword_output_01/iteration_1/haplo1_out",parental_genotypes=parental_genotypes)
pop5=haplotypes_in(input="crossword_output_01/iteration_1/haplo1_out")
pop10=dh(pop=pop1)
pop="pop2"
ss="crossword_output_01/iteration_1/pop2"
draw_haplotypes(haplotypes=pop2,output_folder=ss,parental_genotypes=parental_genotypes,heterozygous=heterozygous,by_chromosomes=by_chromosomes,im_type=im_type)
pop="pop3"
ss="crossword_output_01/iteration_1/pop3_plots"
draw_population(pop=pop3,file=ss,parental_genotypes=parental_genotypes,im_type=im_type)
pheno = phenotype(pop=pop2,qtn_effect=qtn_effect,tgv_only=FALSE,vr=vr,parental_genotypes=parental_genotypes)
pheno2 = pheno
pheno3A = pheno
geno2=get_genotypes(parental_genotypes=parental_genotypes,pop=pop2)
ms2=create_marker_set(genotypes=geno2,N=200,MAF=0.1)
pop="pop4"
pop4=create_families(pop=pop2,S=20,id=pop,chr_length=chr_length)
geno3=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4)
ms3=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4,pre_selected_markers=ms2)
ms4=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4,pre_selected_markers=ms2)
genotypes_out(genotypes=ms3,output="crossword_output_01/iteration_1/pop4_marker_set",level=NA,pop=NA,A_parent=NA,B_parent=NA,parental_genotypes=parental_genotypes)
pheno3=genomic_prediction(train_geno=ms2,train_pheno=pheno2,predict_geno_pop=ms3,method="rrBLUP",level="cross",train_pop=pop2)
pheno=pheno3
percentage=FALSE
gselection1=select(pop=pop4,N=3,level="individual",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_1/pop4")
pheno=pheno3
percentage=FALSE
gselection2=select(pop=pop4,N=1,level="family",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_1/pop4")
phenotype_out(pheno=pheno3A,file="crossword_output_01/iteration_1/pheno3_predicted",im_type=im_type)
pop3_b=mas(pop=pop3,marker=list(c("Aradu.A04_106615237","C",0.5),c("Aradu.A03_7819678","T",0.5)),level="individual",parental_genotypes=parental_genotypes)
pop="pop10"
pop="pop10"
pop10=create_population(c("tg","tr","grg"),chr_length=chr_length,parental_genotypes=parental_genotypes)
pop="pop11"
pop10=create_population(c("tg","tr","grg"),chr_length=chr_length,parental_genotypes=parental_genotypes)
###############################################
###iteration# 2################################
l=2
dir.create(paste0(output_folder,"/iteration_",l))
qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution)
vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode="average")
loci = parental_genotypes[[2]][as.character(qtn_effect$QTN),]
write.table(cbind(qtn_effect,loci),paste0(output_folder,"/iteration_",l,"/selected_QTNs"),col.names=TRUE,row.names=FALSE,quote =FALSE)
pop="pop1"
pop1=crossword::cross(P1="grg",P2=c("tr","tg"),N=10,id=pop,chr_length=chr_length,parental_genotypes=parental_genotypes,isdh=FALSE)
pop="pop2"
pop2=advance(pop=pop1,F=5,level="individual",clevel="individual",id=pop,chr_length=chr_length)
pop="pop3"
pop3=create_families(pop=pop2,S=5,id=pop,chr_length=chr_length)
pheno = phenotype(pop=pop3,qtn_effect=qtn_effect,tgv_only=FALSE,vr=vr,parental_genotypes=parental_genotypes)
percentage=FALSE
pop6=select(pop=pop3,N=3,level="individual",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_2/pop3")
pop9=combine_populations(c(pop1,pop2))
haplotypes_out(pop=pop2,output="crossword_output_01/iteration_2/haplo1_out",parental_genotypes=parental_genotypes)
pop5=haplotypes_in(input="crossword_output_01/iteration_2/haplo1_out")
pop10=dh(pop=pop1)
pop="pop2"
ss="crossword_output_01/iteration_2/pop2"
draw_haplotypes(haplotypes=pop2,output_folder=ss,parental_genotypes=parental_genotypes,heterozygous=heterozygous,by_chromosomes=by_chromosomes,im_type=im_type)
pop="pop3"
ss="crossword_output_01/iteration_2/pop3_plots"
draw_population(pop=pop3,file=ss,parental_genotypes=parental_genotypes,im_type=im_type)
pheno = phenotype(pop=pop2,qtn_effect=qtn_effect,tgv_only=FALSE,vr=vr,parental_genotypes=parental_genotypes)
pheno2 = pheno
pheno3A = pheno
geno2=get_genotypes(parental_genotypes=parental_genotypes,pop=pop2)
ms2=create_marker_set(genotypes=geno2,N=200,MAF=0.1)
pop="pop4"
pop4=create_families(pop=pop2,S=20,id=pop,chr_length=chr_length)
geno3=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4)
ms3=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4,pre_selected_markers=ms2)
ms4=get_genotypes(parental_genotypes=parental_genotypes,pop=pop4,pre_selected_markers=ms2)
genotypes_out(genotypes=ms3,output="crossword_output_01/iteration_2/pop4_marker_set",level=NA,pop=NA,A_parent=NA,B_parent=NA,parental_genotypes=parental_genotypes)
pheno3=genomic_prediction(train_geno=ms2,train_pheno=pheno2,predict_geno_pop=ms3,method="rrBLUP",level="cross",train_pop=pop2)
pheno=pheno3
percentage=FALSE
gselection1=select(pop=pop4,N=3,level="individual",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_2/pop4")
pheno=pheno3
percentage=FALSE
gselection2=select(pop=pop4,N=1,level="family",method="top",pheno=pheno,percentage=percentage,im_type=im_type,file="crossword_output_01/iteration_2/pop4")
phenotype_out(pheno=pheno3A,file="crossword_output_01/iteration_2/pheno3_predicted",im_type=im_type)
pop3_b=mas(pop=pop3,marker=list(c("Aradu.A04_106615237","C",0.5),c("Aradu.A03_7819678","T",0.5)),level="individual",parental_genotypes=parental_genotypes)
pop="pop10"
pop="pop10"
pop10=create_population(c("tg","tr","grg"),chr_length=chr_length,parental_genotypes=parental_genotypes)
pop="pop11"
pop10=create_population(c("tg","tr","grg"),chr_length=chr_length,parental_genotypes=parental_genotypes)
