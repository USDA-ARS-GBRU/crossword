get_genotypes <- function(parental_genotypes,pop,pre_selected_markers)
{
	require(Rcpp)
	if(missing(pop))
	{
	    stop("ERROR: population should be provided")
	}else
	{
	    pop = pop
	}
	if(missing(pre_selected_markers))
	{
	    pre_selected_markers = NA
	}else
	{
	    pre_selected_markers = pre_selected_markers
	}
	haplotypes=pop
	haplo=pop[[2]]
	genotypes = parental_genotypes$genotypes
	genotypes2 = apply(genotypes,1,paste0,collapse="")
	geno = rownames(parental_genotypes$genotypes)
	geno1 = gsub('_.*$','',geno)
	geno2 =  gsub('^.*_','',geno)
	chr = names(haplotypes[[2]][[1]])
	A1 = matrix(nrow = length(rownames(genotypes)), ncol = 2*length(haplo))
	header1 = ""
    cppFunction('CharacterMatrix fun2(List parental_genotypes, List haplotypes, CharacterVector geno1, CharacterVector geno2,CharacterVector chr, CharacterVector genotypes2){
		//get individual IDs from scheme (ind)
        DataFrame scheme;
        scheme = haplotypes[0];
        CharacterVector ind = scheme[0];
        //get haplotypes
        List haplo = haplotypes[1];     
        //get genotypes (genotypes) , snp IDs (geno) and loci (loci)
        CharacterMatrix genotypes = parental_genotypes[0];
        CharacterVector geno = rownames(parental_genotypes[0]);
        DataFrame LOCI;
        LOCI = parental_genotypes[1];
        NumericVector loci = LOCI[0];   
        List mat;
        List pat;
        int m;
        int p;    
        NumericVector matA;
        NumericVector matL;
        NumericVector patA;
        NumericVector patL;
        CharacterVector SNPs;
        std::string SNP1;
        std::string SNP2;
        std::string SNP1_2;
        int Y = ind.length();
        int X = geno.length();
        List INDs;
        List INDS_chr;
        CharacterVector chr_x;
        CharacterMatrix output(X,Y);
        for (int y=0;y<Y;y++)
        {
            INDs = haplo[y];
            for (int x = 0 ;x<X;x++)
            {
               chr_x = geno1[x];
               INDS_chr=INDs[chr_x];
			   INDS_chr = INDS_chr[0];
               mat = INDS_chr["mat"];
			   pat = INDS_chr["pat"];
			   matL = Rcpp::wrap(mat["locations"]);
			   matA = Rcpp::wrap(mat["alleles"]);
			   patL = Rcpp::wrap(pat["locations"]);
			   patA = Rcpp::wrap(pat["alleles"]);        
               m = matL.length();
               p = patL.length();
			   SNPs = genotypes2[x];
               for (int i =0;i<m;i++)
               {
                    if(loci[x] <= matL[i])
                    {      
						SNP1 = (Rcpp::as<std::string>(SNPs)).substr(matA[i]-1,1);
                        break;        
                    }   
               }
               for (int i =0;i<p;i++)
               {
                    if(loci[x] <= patL[i])
                    {
                        SNP2 = (Rcpp::as<std::string>(SNPs)).substr(patA[i]-1,1);
                        break;
                    }
               }
               SNP1_2 = SNP1 + SNP2;
			   output(x,y) = SNP1_2;//Rcpp::as<std::string>(SNP1_2);
            }
        }
        return(output);}')
	A1 = fun2(parental_genotypes,haplotypes,geno1,geno2,chr,genotypes2)
	rownames(A1) = geno
	colnames(A1) = names(haplo)
	scheme = haplotypes[[1]]
    ids = as.character(scheme[scheme[,2]== scheme[nrow(scheme),2],]$id)
    A1 = as.data.frame(A1[,ids])
    colnames(A1) = ids   
    if(!is.na(pre_selected_markers)[[1]])
    {
        if(is.character(pre_selected_markers))
        {
            A1 = get_marker_intersection(genotypes=list(genotypes=A1,gen2phy=parental_genotypes$gen2phy),input_geno_file=pre_selected_markers)[[1]]
        }else
        {
            A1 = get_marker_intersection(genotypes=list(genotypes=A1,gen2phy=parental_genotypes$gen2phy),geno=pre_selected_markers)[[1]]
        }
    } 
    B = parental_genotypes[[2]]
    B2 = as.data.frame(B[rownames(A1),])
    rownames(B2) = rownames(A1)
    colnames(B2) = "gen_loc"    
	return(list(genotypes=A1,gen2phy=B2))
}

