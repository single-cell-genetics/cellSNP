# This module is to parse the gtf file, whoes format can be found at:
# www.ensembl.org/info/website/upload/gff.html

# Note: here we need the order of feature lines like this: gene --> transcript
# --> exon 1 --> exon 2 ..., which is not required for the gtf format. But this
# assumption is usually right for the gtf file downloaded from Ensembl database.

import sys
import subprocess
import numpy as np

class Transcript:
    def __init__(self, chrom, strand, start, stop, tran_id, tran_name="*", 
        biotype="*"):
        """a general purpose transcript object with the basic information.
        """
        self.chrom  = chrom
        self.strand = strand
        self.start  = int(start)
        self.stop   = int(stop)
        self.tranID = tran_id
        self.exons  = np.zeros((0,2), "int")
        self.seglen = None
        self.tranL  = 0
        self.exonNum = 0
        self.biotype = biotype
        self.tranName = tran_name
        

    def add_exon(self, chrom, strand, start, stop):
        if strand != self.strand or chrom != self.chrom:
            print("The exon has different chrom or strand to the transcript.")
            return
        _exon = np.array([start, stop], "int").reshape(1,2)
        self.exons = np.append(self.exons, _exon, axis=0)
        self.exons = np.sort(self.exons, axis=0)
        self.tranL += abs(int(stop) - int(start) + 1)
        self.exonNum += 1


        self.seglen = np.zeros(self.exons.shape[0] * 2 - 1, "int")
        self.seglen[0] = self.exons[0,1]-self.exons[0,0] + 1
        for i in range(1, self.exons.shape[0]):
            self.seglen[i*2-1] = self.exons[i,0]-self.exons[i-1,1] - 1
            self.seglen[i*2] = self.exons[i,1]-self.exons[i,0] + 1

        if ["-","-1","0",0,-1].count(self.strand) > 0:
            self.seglen = self.seglen[::-1]

class Gene:
    def __init__(self, chrom, strand, start, stop, gene_id, gene_name="*",
        biotype="*"):
        """
        """
        self.chrom  = chrom
        self.strand = strand
        self.start  = int(start)
        self.stop   = int(stop)
        self.geneID = gene_id
        self.trans  = []
        self.tranNum = 0
        self.biotype = biotype
        self.geneName = gene_name
        
    def add_transcipt(self, transcript):
        self.trans.append(transcript)
        self.tranNum += 1

    def get_gene_info(self):
        RV = [self.geneID, self.geneName, self.chrom, self.strand, self.start,
              self.stop, self.biotype]
        _trans = []
        for t in self.trans:
            _trans.append(t.tranID)
        RV.append(",".join(_trans))
        return RV

    def add_premRNA(self):
        _tran = Transcript(self.chrom, self.strand, self.start, self.stop, 
                           self.geneID+".p", self.geneName, self.biotype)
        _tran.add_exon(self.chrom, self.strand, self.start, self.stop)
        self.trans.append(_tran)
        self.tranNum += 1
        
    def get_exon_max_num(self):
        exonMax = 0
        for _tran in self.trans:
            exonMax = max(exonMax, _tran.exonNum)
        return exonMax

    def gene_ends_update(self):
        for t in self.trans:
            self.start = min(self.start, np.min(t.exons))
            self.stop  = max(self.stop,  np.max(t.exons))


def parse_attribute(attStr, default="*", 
    ID_tags="ID,gene_id,transcript_id,mRNA_id",
    Name_tags="Name,gene_name,transcript_name,mRNA_name",
    Type_tags="Type,gene_type,gene_biotype,biotype",
    Parent_tags="Parent"):
    """
    Parse attributes in GTF or GFF3

    Parameters
    ----------
    attStr: string
        String containing attributes either in GTF or GFF3 format.
    default: string
        default value for ID, Name, Type and Parent.
    ID_tags: string
        Multiple tags for ID. Use comma for delimit. 
        If multiple tags found, use the last one.
    Name_tags: string
        Multiple tags for Name. Use comma for delimit. 
        If multiple tags found, use the last one.
    Type_tags: string
        Multiple tags for Type. Use comma for delimit. 
        If multiple tags found, use the last one.
    Parent_tags: string
        Multiple tags for Parent. Use comma for delimit. 
        If multiple tags found, use the last one.

    Returns
    -------
    RV: library of string
        Library of all tags, always including ID, Name, Type, Parenet.
    """
    RV = {}
    RV["ID"] = default
    RV["Name"] = default
    RV["Type"] = default
    RV["Parent"] = default
    ID_tags = ID_tags.split(",")
    Name_tags = Name_tags.split(",")
    Type_tags = Type_tags.split(",")
    Parent_tags = Parent_tags.split(",")

    attList = attStr.rstrip().split(";")
    for att in attList:
        while len(att) > 0 and att[0] == " ": 
            att = att[1:]
        if len(att) == 0: 
            continue
        if att.find("=") > -1:
            _att = att.split("=") #GFF3
        else:
            _att = att.split(" ") #GTF

        if len(_att) < 2:
            print("Can't pase this attribute: %s" %att)
            continue

        if _att[1][0] == '"':
            _att[1] = _att[1].split('"')[1]

        if ID_tags.count(_att[0]) == 1:
            RV["ID"] = _att[1]
        elif Name_tags.count(_att[0]) == 1:
            RV["Name"] = _att[1]
        elif Type_tags.count(_att[0]) == 1:
            RV["Type"] = _att[1]
        elif Parent_tags.count(_att[0]) == 1:
            RV["Parent"] = _att[1]
        else: RV[_att[0]] = _att[1]

    return RV

def load_genes(anno_file, comments="#,>", geneTag="gene", 
        tranTag="transcript,mRNA", exonTag="exon"):
    """
    Load genes from gtf or gff3 file.

    Parameters
    ----------
    anno_file: str
        path for the annotation file in GTF or GFF3 format.
    comments: string
        Multiple comments. Use comma for delimit. 
    geneTag: string
        Multiple tags for gene. Use comma for delimit. 
    tranTag: string
        Multiple tags for transcript. Use comma for delimit. 
    exonTag: string
        Multiple tags for exon. Use comma for delimit. 

    Return
    ------
    genes: list of ``pyseqlib.Gene``
        a list of loaded genes
    """

    if anno_file.endswith(".gz") or anno_file.endswith(".gzip"):
        import gzip
        anno_in = gzip.open(anno_file, "rb")
    else:
        anno_in = open(anno_file, "r")

    geneTag = geneTag.split(",")
    tranTag = tranTag.split(",")
    exonTag = exonTag.split(",")
    comments = comments.split(",")
    
    genes = []
    _gene = None
    for _line in anno_in:
        try:
            _line = _line.decode('utf-8')
        except AttributeError:
            pass
        if comments.count(_line[0]):
            continue
            
        aLine = _line.split("\t")
        if len(aLine) < 8:
            continue
        elif geneTag.count(aLine[2]) == 1:
            if _gene is not None: 
                genes.append(_gene)

            RVatt = parse_attribute(aLine[8], ID_tags="ID,gene_id",
                Name_tags="Name,gene_name")
            _gene = Gene(aLine[0], aLine[6], aLine[3], aLine[4],
                RVatt["ID"], RVatt["Name"], RVatt["Type"])

        elif tranTag.count(aLine[2]) == 1:
            RVatt = parse_attribute(aLine[8],ID_tags="ID,transcript_id,mRNA_id",
                Name_tags="Name,transcript_name,mRNA_name")
            _tran  = Transcript(aLine[0], aLine[6], aLine[3], aLine[4],
                RVatt["ID"], RVatt["Name"], RVatt["Type"])

            if _gene is not None:
                _gene.add_transcipt(_tran)
            else:
                print("Gene is not ready before transcript.")

        elif exonTag.count(aLine[2]) == 1:
            if aLine[0] != _gene.trans[-1].chrom:
                print("Exon from a different chrom of transcript.")
                continue
            if aLine[6] != _gene.trans[-1].strand:
                print("Exon from a different strand of transcript.")
                continue
            if _gene is not None and len(_gene.trans) > 0:
                _gene.trans[-1].add_exon(aLine[0], aLine[6], aLine[3], aLine[4])
                # _gene.gene_ends_update()
            else:
                print("Gene or transcript is not ready before exon.")
    
    anno_in.close()
    if _gene is not None: 
        genes.append(_gene)

    return genes


def save_genes(out_file, genes, atype="GFF3", tags="gene,mRNA,exon"):
    """Save genes into file in GFF3 or GTF format.

    Parameters
    ----------
    out_file: string
        full path for output file
    genes: list of ``pyseqlib.Gene``
        a list of loaded genes
    atype: string
        type of output file: GTF or GFF3
    tags: string
        tags for gene, mRNA, exon, respectively. Use comma for delimit
    """
    if out_file.endswith(".gz"):
        out_file = out_file[:-3]
    elif out_file.endswith(".gzip"):
        out_file = out_file[:-5]
    fid = open(out_file, "w")
    fid.writelines("#%s file produced by savegene.\n" %atype)

    tags = tags.split(",")
    aLine = [".", ".", ",", "0", "0", ".", "+", ".", ""]
    for g in genes:
        aLine[0] = g.chrom
        aLine[2] = tags[0]
        aLine[3] = str(g.start)
        aLine[4] = str(g.stop)
        aLine[6] = g.strand
        if atype.upper() == "GFF3":
            aLine[8] = "ID=%s;gene_id=%s" %(g.geneID, g.geneID)
            if g.geneName != "*" and  g.geneName != "#":
                aLine[8] += ";gene_name=%s" %g.geneName
            if g.biotype != "*" and  g.biotype != "#":
                aLine[8] += ";gene_type=%s" %g.biotype
        else:
            aLine[8] = "gene_id \"%s\"" %g.geneID
            if g.geneName != "*" and  g.geneName != "#":
                aLine[8] += "; gene_name \"%s\"" %g.geneName
            if g.biotype != "*" and  g.biotype != "#":
                aLine[8] += "; gene_type \"%s\"" %g.biotype
        fid.writelines("\t".join(aLine) + "\n")

        for t in g.trans:
            aLine[0] = t.chrom
            aLine[2] = tags[1]
            aLine[3] = str(t.start)
            aLine[4] = str(t.stop)
            aLine[6] = t.strand
            if atype.upper() == "GFF3":
                aLine[8] = "ID=%s;Parent=%s" %(t.tranID, g.geneID)
            else:
                aLine[8] = "gene_id \"%s\"; transcript_id \"%s\"" %(g.geneID, 
                    t.tranID)
            fid.writelines("\t".join(aLine) + "\n")

            exons = t.exons
            for i in range(exons.shape[0]):
                aLine[0] = t.chrom
                aLine[2] = tags[2]
                aLine[3] = str(exons[i,0])
                aLine[4] = str(exons[i,1])
                aLine[6] = t.strand
                if atype.upper() == "GFF3":
                    aLine[8] = "ID=%s.%d;Parent=%s" %(t.tranID, i+1, t.tranID)
                else:
                    aLine[8] = "gene_id \"%s\"; transcript_id \"%s\"" %(g.geneID, 
                        t.tranID)
                fid.writelines("\t".join(aLine) + "\n")
    fid.close()
    
    bashCommand = "gzip -f %s" %(out_file) 
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = pro.communicate()[0]
        

def savegene(**kwargs):
    print("This function is deprecated, please use save_genes")
    return save_genes(**kwargs)

def loadgene(**kwargs):
    print("This function is deprecated, please use load_genes")
    return save_genes(**kwargs)