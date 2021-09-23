#!/usr/local/bin/python

#usage:
#python denovo.py [conventional|all] [child|all] input.vcf.gz child_name output.denovo.vcf.gz

#ex:
#python denovo.py conventional child variant_calls/trio1.vcf.gz HG00404 denovo_calls/trio1.denovo.vcf.gz

#child/all : output denovo vcf file contains only the child / child and the parents

#conventional/all : only the loci where parents carry reference allele and the child carries at least one alternative 
# allele is considered denovo / all loci where the child carries an allele not present in parents are considered denovo

import sys
import gzip

novelty=sys.argv[1]
outtype=sys.argv[2]
infile=sys.argv[3]
child=sys.argv[4]
outfile=sys.argv[5]


t= gzip.open(infile,"r")
denovo= gzip.open(outfile, "w")


for line in t:
    if line.startswith(b'#CHROM'):
        denovo.write(b'# De novo variants are subsetted using denovo.py')
        denovo.write(str.encode('\n#'+str(sys.argv)+"\n"))
        header=line.decode().strip().split()
        child_idx=header.index(child)
        parent_idx=[9,10,11]
        parent_idx.remove(child_idx)

        if  outtype=="all":
            denovo.write(line)

        elif outtype=="child":
            line_list=line.decode().strip().split()
            out_idx=[0,1,2,3,4,5,6,7,8]+[child_idx]
            outline_list=[line_list[i] for i in out_idx]
            outline="\t".join(outline_list)
            denovo.write(str.encode(outline+"\n"))
        else:
            raise ValueError('The output file can contain either \"child\" or \"all\" individuals. Check the second argument')

    elif line.startswith(b'#') :
        denovo.write(line)
    
    else:
        line_list=line.decode().strip().split()
        if line_list[6] != 'PASS' and line_list[6] !='.':
            continue

        else:
            alleles_c=[line_list[child_idx][i] for i in [0,2]]
            alleles_p1=[line_list[parent_idx[0]][i] for i in [0,2]]
            alleles_p2=[line_list[parent_idx[1]][i] for i in [0,2]]

            if novelty== "conventional" and outtype=="all":
                if alleles_p1 == ["0","0"] and alleles_p2 == ["0","0"] and alleles_c != ["0","0"]:
                    outline="\t".join(line_list)
                    denovo.write(str.encode(outline+"\n"))
            
            elif novelty== "conventional" and outtype=="child":
                if alleles_p1 == ["0","0"] and alleles_p2 == ["0","0"] and alleles_c != ["0","0"]:
                    outline_list=[line_list[i] for i in out_idx]
                    outline="\t".join(outline_list)
                    denovo.write(str.encode(outline+"\n"))

            elif novelty== "all" and outtype=="all":
                p=alleles_p1+alleles_p2
                if (alleles_c[0] not in p) or (alleles_c[1] not in p) :
                    outline="\t".join(line_list)
                    denovo.write(str.encode(outline+"\n"))
        
            elif novelty== "all" and outtype=="child":
                p=alleles_p1+alleles_p2
                if (alleles_c[0] not in p) or (alleles_c[1] not in p) :
                    outline_list=[line_list[i] for i in out_idx]
                    outline="\t".join(outline_list)
                    denovo.write(str.encode(outline+"\n"))

            else:
                raise ValueError('De novo vairant type can be either \"conventional\" or \"all\". Check the first argument')

t.close()
denovo.close()
