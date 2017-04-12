#! /usr/bin/python
import sys
import os
import argparse
import subprocess
import datetime
from Bio import SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def get_ave_read_length(filename):
    """
    Parameters
    ----------
    filename : str
        Name of file BAM/SAM file.

    Returns
    -------
    int
        Average read length from stats.

    """
    p1 = subprocess.Popen(["samtools", "stats", filename], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["cut", "-f", "2-"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(["grep", "average length"], stdin=p2.stdout, stdout=subprocess.PIPE)
    p4 = subprocess.Popen(["cut", "-f", "2-"], stdin=p3.stdout, stdout=subprocess.PIPE)
    return int(p4.stdout.read().strip())

def is_paired_end(filename):
    """
    Parameters
    ----------
    filename : str
        Name of BAM/SAM file.

    Returns
    -------
    boolean
        True if reads are paired-end, False otherwise. 

    """
    return int(subprocess.check_output(["samtools", "view", "-c", "-f", "1", filename], universal_newlines=True).strip())

def run_TCR_pipe(genome, output, bam, unmapped, bases, strand, num_iterations, threshold_score, 
    min_overlap, rsem, bowtie2, single_cell, path, sumF, lowQ, samtools, trim):
    """
    Parameters
    ----------

    Returns
    -------

    """
    print("~" * 100)
    print(args.genome)

    check_parameters(genome, strand, single_cell, path, sumF)
    if single_cell == True:
        # TODO: Fix this, won't work for SE
        sys.exit(0)
    if path == "./":
        path = os.getcwd()
    if not path.endswith("/"):
        path = path + "/"
    final_stat_dict = {}
    tcrFout = open(sumF + ".TCRs.txt", "w")
    tcrFout.write("cell\tChain\tStatus\tRank of TCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in TCR\t#reads in CDR3\t#reads in V\t#reads in J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n")

    for cell_folder in os.listdir(path):
        full_path = path + cell_folder + "/"
        if os.path.exists(full_path) and os.path.isdir(full_path):
            print(str(datetime.datetime.now()) + " Working on: " + cell_folder)
            found, nbam, nunmapped, noutput = format_files(full_path, bam, unmapped, output)
            if not found:
                print(str(datetime.datetime.now()) +
                    " There is not a bam or unmapped file in\n"
                    "this folder, moving to the next folder.")
            else:
                curr_folder = os.path.abspath(os.path.dirname(sys.argv[0])) + '/'
                reconstruction = curr_folder + '/vdj.alignment'
                if genome == "hg38":
                    extra = "id.name"
                else:
                    extra = "gene.id"
                fasta = curr_folder + "Data/{0}/{0}.TCR.fa".format(genome)
                bed = curr_folder + "Data/{0}/{0}.TCR.bed".format(genome)
                mapping = curr_folder + "Data/{0}/{0}.{1}.mapping.TCR.txt".format(genome, extra)
                aaF = curr_folder + "Data/{0}/{0}.TCR.conserved.AA.txt".format(genome)
                ref_ind = curr_folder + "Data/{0}/index/{0}".format(genome)
                run_single_cell(fasta, bed, noutput, nbam, nunmapped, mapping, bases, strand, 
                              reconstruction, aaF , num_iterations, threshold_score,
                              min_overlap, rsem, bowtie2, lowQ, samtools, ref_ind, trim)

                add_cell_to_TCR_sum(cell_folder, noutput, tcrFout)
                final_stat_dict = add_to_stat_dict(noutput, cell_folder, final_stat_dict)
    sumFout = open(sumF + ".summary.txt", "w")
    sumFout.write("sample\talpha\tbeta\n")
    for cell in sorted(final_stat_dict):
        fout = cell + "\t" + final_stat_dict[cell]["alpha"] + "\t" + final_stat_dict[cell]["beta"] + "\n"
        sumFout.write(fout)
    sumFout.close()


def add_cell_to_TCR_sum(cell_folder, noutput, tcrFout):
    """
    Adds the entries in a cell folder's .summary.txt file 
    to the main .TCRs.txt file, which will contain information
    from all files.

    Parameters
    ----------
    cell_folder : str
    noutput : str
    tcrFout : file

    Returns
    -------

    """
    print("*" * 100)
    print(type(cell_folder), type(noutput), type(tcrFout))
    print(cell_folder, noutput, tcrFout)
    if os.path.isfile(noutput + ".summary.txt"):
        curr_out = open(noutput + ".summary.txt", "r")
        # move past header
        curr_out.readline()

        # first data line
        l = curr_out.readline()
        while l != "":
            newL = cell_folder + "\t" + l
            tcrFout.write(newL)
            l = curr_out.readline()
        curr_out.close()


def get_chain_message(stat, msg):
    """
    Parameters
    ----------
    stat : str
    msg : str

    Returns
    -------
    str  

    """
    if stat == "Productive":
        return stat
    # Gunjan - Why does unproductive take priority? Sometimes, 
    # unproductive for when only J found. 
    if stat.startswith("Unproductive") and msg != "Productive":
        return "Unproductive"
    if stat == "Failed reconstruction - reached maximum number of iterations" and msg == "None":
        return "Failed - reconstruction didn\'t converge"
    if stat == "Failed reconstruction - V and J segment do not overlap" and msg == "None":
        return "Failed - V and J reconstruction don\'t overlap"


def update_chain_message(junctions):
    """
    Parameters
    ----------
    junctions : str

    Returns
    -------
    str

    """
    if os.path.isfile(junctions) == True and os.stat(junctions).st_size != 0:
        return "Failed - found V and J segments but wasn\'t able to extend them"
    else:
        return "Failed - didn\'t find any V and J segments in original mapping"


def add_to_stat_dict(noutput, cell_folder, final_stat_dict):
    """
    Parameters
    ----------
    noutput : 
    cell_folder : 
    final_stat_dict :

    Returns
    -------
    final_stat_dict : 

    """
    if cell_folder in final_stat_dict:
        print("Error! {} appear more than once in final stat dictionary".format(cell_folder))
    failed_message = "Failed - found V and J segments but wasn\'t able to extend them"
    final_stat_dict[cell_folder] = {"alpha": failed_message, "beta": failed_message}
    alphaJunc = noutput + '.alpha.junctions.txt'
    betaJunc = noutput + '.beta.junctions.txt'
    if os.path.isfile(noutput + ".summary.txt"):
        currOut = open(noutput + ".summary.txt", "r")
        msgA = "None"
        msgB = "None"
        currOut.readline()
        l = currOut.readline()
        while l != '':
            lArr = l.strip('\n').split('\t')
            chain = lArr[0]
            stat = lArr[1]
            
            # only update the message if get_chain_message() doesn't return None
            if chain == "alpha":
                msgA = get_chain_message(stat, msgA) or msgA
            else:
                msgB = get_chain_message(stat, msgB) or msgB
            l = currOut.readline()


        currOut.close()
        if msgA == 'None':
            msgA = update_chain_message(alphaJunc)
        if msgB == 'None':
            msgB = update_chain_message(betaJunc)
    else:
        # Can else case be added here, to generalize with update_chain_message()?
        if os.path.isfile(betaJunc) == True:
            # What about else cases here???
            if os.stat(betaJunc).st_size == 0:
                msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(alphaJunc) == True):
            if os.stat(alphaJunc).st_size == 0:
                msgA = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgA = 'Failed - didn\'t find any V and J segments in original mapping'

    final_stat_dict[cell_folder]["alpha"] = msgA
    final_stat_dict[cell_folder]["beta"] = msgB
    return final_stat_dict


def format_path(full_path, file_name):
    """
    Parameters
    ----------
    full_path : 
    file_name : 

    Returns
    -------
    str 

    """
    if file_name.startswith("/"):
        return full_path + file_name[1:]
    if file_name.startswith("./"):
        return full_path + file_name[2:]
    return full_path + file_name


def format_files(full_path, bam, unmapped, output):
    """
    Parameters
    ----------
    full_path :
    bam : 
    unmapped : 
    output : 

    Returns
    -------
    found : 
    nbam : 
    nunmapped : 
    noutput :

    """
    found = True
    nbam = format_path(full_path, bam)
    nunmapped = format_path(full_path, unmapped)
    if os.path.isfile(nunmapped) and os.path.isfile(nbam):
        noutput = make_output_dir(output, full_path)
    else:
        noutput = output
        found = False
    return found, nbam, nunmapped, noutput


def make_output_dir(output, full_path):
    """
    Parameters
    ----------
    output : 
    full_path : 

    Returns
    -------
    noutput : 

    """
    noutput = format_path("", output)
    if output.endswith('/'):
        noutput = output[:-1]
    if output.find('/') != -1:
        outArr = noutput.split('/')
        currPath = full_path
        for i in range(len(outArr)-1):
            currPath = currPath + outArr[i] + '/'
            if not os.path.exists(currPath):
                os.makedirs(currPath)
    noutput = full_path + noutput
    return noutput


def run_single_cell(fasta, bed, output, bam, unmapped, mapping, bases, strand, 
                  reconstruction, aaF, num_iterations, threshold_score, 
                  min_overlap, rsem, bowtie2, lowQ, samtools, ref_ind, trim):
    """
    Parameters
    ----------

    Returns
    -------
    None

    """
    idNameDict = make_id_name_dict(mapping)
    fastaDict = make_fasta_dict(fasta)
    vdjDict = make_VDJ_bed_dict(bed, idNameDict)

    ave_length = get_ave_read_length(unmapped)
    paired_end = is_paired_end(bam)
    print("average length " + str(ave_length))

    if paired_end:
        print(str(datetime.datetime.now()) + " Pre-processing alpha chain")
        unDictAlpha = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, 
                                   bases, 'A', strand, lowQ)
        print(str(datetime.datetime.now()) + " Pre-processing beta chain")
        unDictBeta = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, 
                                  bases, 'B', strand, lowQ)
    else:
        print(str(datetime.datetime.now()) + " Pre-processing alpha chain")
        print(str(datetime.datetime.now()) + " Pre-processing beta chain")
        analyze_chain_single_end(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 
                              strand, lowQ, bowtie2, ref_ind, trim)    
        unDictAlpha = write_unmapped_reads_to_dict_SE(unmapped)
        unDictBeta = unDictAlpha
    
    print(str(datetime.datetime.now()) + " Reconstructing alpha chains")
    subprocess.call([reconstruction, output + '.alpha.mapped.and.unmapped.fa', 
                     output + '.alpha.junctions.txt', output + '.reconstructed.junctions.alpha.fa', 
                     str(num_iterations), str(threshold_score), str(min_overlap)])
    print(str(datetime.datetime.now()) + " Reconstructing beta chains")
    subprocess.call([reconstruction, output + '.beta.mapped.and.unmapped.fa', 
                     output + '.beta.junctions.txt', output + '.reconstructed.junctions.beta.fa', 
                     str(num_iterations), str(threshold_score), str(min_overlap)])
    
    print(str(datetime.datetime.now()) + " Creating full TCR sequencing")
    fullTcrFileAlpha = output + '.alpha.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.alpha.fa'
    create_TCR_full_output(fastaDict, tcrF, fullTcrFileAlpha, bases, idNameDict)
    fullTcrFileBeta = output + '.beta.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.beta.fa'
    create_TCR_full_output(fastaDict, tcrF, fullTcrFileBeta , bases, idNameDict)
    
    print(str(datetime.datetime.now()) + " Running RSEM to quantify expression of all possible isoforms")
    outDirInd = output.rfind('/')
    if outDirInd != -1:
        outDir = output[:outDirInd+1]
    else:
        outDir = os.getcwd()
    run_rsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools)
    pick_final_isoforms(fullTcrFileAlpha, fullTcrFileBeta, output)
    bestAlpha = output + '.alpha.full.TCRs.bestIso.fa'
    bestBeta = output + '.beta.full.TCRs.bestIso.fa'
    
    print(str(datetime.datetime.now()) + " Finding productive CDR3")
    aaDict = make_AA_Dict(aaF)
    fDictAlpha, fDictBeta = {}, {}
    if os.path.isfile(bestAlpha):
        fDictAlpha = find_CDR3(bestAlpha, aaDict, fastaDict)
    if os.path.isfile(bestBeta):
        fDictBeta = find_CDR3(bestBeta, aaDict, fastaDict)

    betaRsemOut = output + '.beta.rsem.out.genes.results'
    alphaRsemOut = output + '.alpha.rsem.out.genes.results'
    alphaBam = output + '.alpha.rsem.out.transcript.sorted.bam'
    betaBam = output + '.beta.rsem.out.transcript.sorted.bam'
    print(str(datetime.datetime.now()) + " Writing results to summary file")

    # no modifications to unDictAlpha/unDictBeta after this point,
    # just checking its contents
    make_single_cell_output_file(fDictAlpha, fDictBeta, output, betaRsemOut, alphaRsemOut, 
                             alphaBam, betaBam, fastaDict, unDictAlpha, unDictBeta, idNameDict)


def analyze_chain_single_end(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, strand, 
                          lowQ, bowtie2, ref_ind, trim):
    """
    Parameters
    ----------

    Returns
    -------
    None

    """
    mappedReadsDictAlpha = dict()
    mappedReadsDictBeta = dict()
 
    fastq = unmapped + ".fq"
    fastq2 = bam + ".fq"

    sam = fastq + '.trimmed.sam'
    temp1 = sam + '.temp1'
    temp2 = sam + '.temp2'

    if bowtie2 != '':
        if bowtie2.endswith('/'):
            bowtieCall = bowtie2 + 'bowtie2'
        else:
            bowtieCall = bowtie2 + '/bowtie2'
    else:
        bowtieCall = 'bowtie2'

    if not os.path.isfile(fastq):
        fastq_file = open(fastq, "w+")
        subprocess.call(["samtools", "bam2fq", unmapped], stdout=fastq_file)
        fastq_file.close()
    if not os.path.isfile(fastq2):
        fastq2_file = open(fastq2, "w+")
        subprocess.call(["samtools", "fastq", bam], stdout=fastq2_file)
        fastq2_file.close()

    # need fastq format of unmapped.bam here
    # change this to use trim length as a parameter too
    subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x', ref_ind, '-U', fastq, '-S', temp1, "--trim3", str(trim)])
    subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x', ref_ind, '-U', fastq, '-S', temp2, "--trim5", str(trim)])
    subprocess.call(["samtools", "merge", "-f", sam, temp1, temp2])

    sam_filtered = output + ".trimmed.filtered.sam"
    samF = open(sam_filtered, "w")
    subprocess.call(["samtools", "view", "-b", "-F", "4", sam], stdout=samF)
    samF.close()
    subprocess.call(["rm", temp1, temp2])

    print("DONE WITH BOWTIE!!!!")

    if os.path.isfile(sam_filtered):
        mappedReadsDictAlpha = find_reads_and_segments(sam_filtered, mappedReadsDictAlpha, idNameDict, 'A')
        mappedReadsDictBeta = find_reads_and_segments(sam_filtered, mappedReadsDictBeta, idNameDict, 'B')

    # adding originally mapped reads to dictionaries
    sam_mapped = output + ".sorted.sam"
    subprocess.call(["samtools", "view", "-h", "-o", sam_mapped, bam])
    if os.path.isfile(sam_mapped):
        mappedReadsDictAlpha = find_reads_and_segments(sam_mapped, mappedReadsDictAlpha, idNameDict, 'A')
        mappedReadsDictBeta = find_reads_and_segments(sam_mapped, mappedReadsDictBeta, idNameDict, 'B')

    alphaOut = output + '.alpha.junctions.txt'
    alphaOutReads = output + '.alpha.mapped.and.unmapped.fa'
    betaOutReads = output + '.beta.mapped.and.unmapped.fa'
    betaOut = output + '.beta.junctions.txt'
    write_junction_file_SE(mappedReadsDictAlpha, idNameDict, alphaOut, fastaDict, bases, 'alpha')
    write_junction_file_SE(mappedReadsDictBeta, idNameDict, betaOut, fastaDict, bases, 'beta')

    write_reads_file_SE(mappedReadsDictAlpha, alphaOutReads, fastq, fastq2)
    write_reads_file_SE(mappedReadsDictBeta, betaOutReads, fastq, fastq2)


def write_reads_file_SE(mappedReadsDict, outReads, fastq, fastq2):
    """
    Parameters
    ----------
    mappedReadsDict : 
    outReads : 
    fastq : 
    fastq2 : 

    Returns
    -------
    None

    """
    seen = []
    if fastq.endswith('.gz'):
        subprocess.call(['gunzip', fastq])
        newFq = fastq.replace('.gz','')
    else:
        newFq = fastq
    out = open(outReads, 'w')
    fqF = open(newFq, 'rU')
    for record in SeqIO.parse(fqF, 'fastq'):
        if record.id in mappedReadsDict and record.seq not in seen:
            newRec = SeqRecord(record.seq, id = record.id, description = '')
            SeqIO.write(newRec,out,'fasta')
            # seen.append(record.seq)
    fqF.close()
    fqF2 = open(fastq2, 'rU')
    for record in SeqIO.parse(fqF2, 'fastq'):
        if record.id in mappedReadsDict and record.seq not in seen:
            newRec = SeqRecord(record.seq, id = record.id, description = '')
            SeqIO.write(newRec,out,'fasta')
            # seen.append(record.seq)
    fqF2.close()
    out.close()
    if fastq.endswith('.gz'):
        subprocess.call(['gzip',newFq])



def write_junction_file_SE(mappedReadsDict,idNameDict, output, fastaDict, bases, chain):
    """
    Parameters
    ----------
    mappedReadsDict : 
    idNameDict : 
    output : 
    fastaDict : 


    Returns
    -------

    """
    out = open(output, 'w')
    vSegs = []
    jSegs = []
    cSegs = []
    for read in mappedReadsDict:
        for seg in mappedReadsDict[read]:
            if idNameDict[seg].find('V') != -1:
                if seg not in vSegs:
                    vSegs.append(seg)
            elif idNameDict[seg].find('J') != -1:
                if seg not in jSegs:
                    jSegs.append(seg)
            elif idNameDict[seg].find('C') != -1:
                if seg not in cSegs:
                    cSegs.append(seg)
            else:
                print "Error! not V/J/C in fasta dict"
    if len(vSegs) == 0:
        print "Did not find any V segments for " + chain + " chain"
    else:
        if len(cSegs) == 0:
            print "Did not find any C segments for " + chain + " chain"
            cSegs = ['NA']
        if len(jSegs) == 0:
            print "Did not find any J segments for " + chain + " chain"
            jSegs = ['NA']
        for vSeg in vSegs:
            for jSeg in jSegs:
                for cSeg in cSegs:
                    add_segment_to_junction_file_SE(vSeg,jSeg,cSeg,out,fastaDict, bases, idNameDict)
    out.close()



def add_segment_to_junction_file_SE(vSeg,jSeg,cSeg,out,fastaDict, bases, idNameDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    vSeq = fastaDict[vSeg]
    if jSeg != 'NA':
        jName = idNameDict[jSeg]
        jSeq = fastaDict[jSeg]
    else:
        jSeq = ''
        jName = 'NoJ'
    if cSeg != 'NA':
        cName = idNameDict[cSeg]
        cSeq = fastaDict[cSeg]
    else:
        cName = 'NoC'
        cSeq = ''
    jcSeq = jSeq + cSeq
    lenSeg = min(len(vSeq),len(jcSeq))
    print("-------------", len(vSeq),len(jcSeq), "-----------")
    if bases != -10:
        if lenSeg < bases:
            sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, taking the length' \
                    'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
            sys.stdout.flush()
        else:
            lenSeg = bases
    jTrim = jcSeq[:lenSeg]
    vTrim = vSeq[-1*lenSeg:]
    junc = vTrim + jTrim
    recordName = vSeg + '.' + jSeg + '.' + cSeg + '(' + idNameDict[vSeg] + '-' + jName + '-' + cName + ')'
    record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
    SeqIO.write(record,out,'fasta')


def find_reads_and_segments(samF, mappedReadsDict, idNameDict, chain):
    """
    Parameters
    ----------
    samF : str
        Name of samfile outputed by Bowtie2.
    mappedReadsDict : dict
        Initially empty.
    idNameDict : dict
        Contains (sequence ID, name) pairs i.e. ('ENST00000390536': 'TRAJ1').
    chain : str
        Either "A" or "B".

    Returns
    -------
    mappedReadsDict : dict
        Contains (readName, [segment]) pairs i.e. 
        ('NS500531:34:H2TMNBGXX:1:11108:15464:19131': ['ENST00000390428']).
    """
    samFile = pysam.AlignmentFile(samF,'r')
    readsIter = samFile.fetch(until_eof = True)
    for read in readsIter:
        if read.is_unmapped == False:
            seg = samFile.getrname(read.reference_id)
            if seg in idNameDict:
                if idNameDict[seg].find(chain) != -1:
                    readName = read.query_name
                    if readName not in mappedReadsDict:
                        mappedReadsDict[readName] = []
                    if seg not in mappedReadsDict[readName]:
                        if chain in idNameDict[seg]:
                            mappedReadsDict[readName].append(seg)
    samFile.close()
    return mappedReadsDict



def make_single_cell_output_file(alphaDict, betaDict, output, betaRsem, alphaRsem, alphaBam, betaBam, fastaDict,
                             unDictAlpha, unDictBeta, idNameDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    outF = open(output + '.summary.txt', 'w')
    outF.write('Chain\tStatus\tRank of TCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in TCR\t#reads in CDR3\t#reads in V\t#reads in J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n')
    if (len(alphaDict) > 0):
        write_chain(outF, 'alpha',alphaDict,alphaRsem, alphaBam, fastaDict,unDictAlpha, output, idNameDict)
    if (len(betaDict) > 0):
        write_chain(outF,'beta',betaDict, betaRsem, betaBam, fastaDict, unDictBeta, output, idNameDict)
    outF.close()


def write_chain(outF, chain,cdrDict,rsemF, bamF, fastaDict, unDict, output, idNameDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    writtenArr = []
    if os.path.exists(rsemF):
        noRsem = False
        (rsemDict,unRsemDict) = make_rsem_dict(rsemF, cdrDict)
    else:
        noRsem = True
    for tcr in cdrDict:
        jStart = -1
        cdrInd = -1
        cInd = -1
        if cdrDict[tcr]['stat'] == 'Productive':
            isProd = True
        else:
            isProd = False
        fLine = chain + '\t' + cdrDict[tcr]['stat'] + '\t'
        if noRsem:
            rank = 'NA'
        else:
            rank = get_rank(tcr, rsemDict, unRsemDict, isProd, noRsem)
        fLine += str(rank) + '\t'
        nameArr = tcr.split('.')
        fLine += nameArr[0] + '\t' + nameArr[1] + '\t' + nameArr[2] + '\t'
        fLine += cdrDict[tcr]['CDR3 NT'] + '\t' + cdrDict[tcr]['CDR3 AA'] + '\t'
        fullSeq = cdrDict[tcr]['Full Seq'].upper()
        if not noRsem:
            totalCount = find_counts_in_region(bamF, 0, len(fullSeq), tcr)
            fLine += str(totalCount) + '\t'
            cSeq = fastaDict[nameArr[5]].upper()
            cInd = fullSeq.find(cSeq)
            if cInd == -1:
                sys.stderr.write(str(datetime.datetime.now()) + 'Error! could not find C segment sequence in the full sequence\n')
                sys.stderr.flush()
            cCounts = find_counts_in_region(bamF, cInd, len(fullSeq), tcr)
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                cdrInd = fullSeq.find(cdrDict[tcr]['CDR3 NT'].upper())
            else:
                cdrInd = -1
            if ((cdrInd == -1) & (cdrDict[tcr]['CDR3 NT'] != 'NA')):
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! Cound not find CDR3 NT sequence in the full sequence\n')
                sys.stderr.flush()
            if cdrInd != -1:
                cdrCounts = find_counts_in_region(bamF, cdrInd, cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr)
                jStart = cdrInd + len(cdrDict[tcr]['CDR3 NT'])

                jCounts = find_counts_in_region(bamF, jStart, cInd, tcr)
                vCounts = find_counts_in_region(bamF, 0, cdrInd, tcr)
                fLine += str(cdrCounts) + '\t' + str(vCounts) + '\t' + str(jCounts) + '\t' + str(cCounts) + '\t'
            else:
                fLine += 'NA\tNA\tNA\t' + str(cCounts) + '\t'
            vId = nameArr[3]
            jId = nameArr[4]
            cId = nameArr[5]
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                (unDictRatioCDR, unCDRcount) = get_un_dict_ratio(bamF, cdrInd , cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr, unDict)
                (unDictRatioALL, unAllcount) = get_un_dict_ratio(bamF, 0 , len(fullSeq), tcr, unDict)
                fLine += str(unDictRatioALL) + '\t' + str(unAllcount) + '\t' + str(unDictRatioCDR) + '\t' + str(unCDRcount) + '\t'
            else:
                fLine += 'NA\tNA\tNA\tNA\t'
            writtenArr.append(vId)
            writtenArr.append(jId)
            fLine += vId + '\t' + jId + '\t' + cId + '\n'
        else:
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' + nameArr[3] + '\t' + nameArr[4] + '\t' + nameArr[5] + '\n'
            #print fLine
        outF.write(str(fLine))
    write_failed_reconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict )

def write_failed_reconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    recF = output + '.reconstructed.junctions.' + chain + '.fa'
    if os.path.isfile(recF):
        f = open(recF, 'rU')
        segDict = dict()
        for tcrRecord in SeqIO.parse(f, 'fasta'):
            tcrSeq = str(tcrRecord.seq)
            if tcrSeq.find('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN') != -1:
                status = 'Failed reconstruction - reached maximum number of iterations'
                segDict = add_segments_to_dict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
            elif tcrSeq.find('NNNN') != -1:
                status = 'Failed reconstruction - V and J segment do not overlap'
                segDict = add_segments_to_dict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
        f.close()
        if len(segDict) > 0:
            write_seg_dict(segDict, outF, chain)


def write_seg_dict(segDict, outF, chain):
    """
    Parameters
    ----------

    Returns
    -------

    """
    for seg in segDict:
        currDict = segDict[seg]
        pairs = ''
        for pair in currDict['pairs']:
            pairs += pair + '.'
        pairs = pairs[:-1]
        if currDict['len'] > 0:
            fLine = chain + '\t' + currDict['status'] + '\t'
            rank = find_curr_rank(segDict, seg, currDict['len'])
            fLine += str(rank) + '\t'
            if currDict['type'] == 'V':
                fLine += currDict['name'] + '\t' + 'paired with: ' + pairs + '\t'
            else:
                fLine +=  'paired with: ' + pairs + '\t' + currDict['name'] + '\t'
            fLine += 'NA\t' + currDict['seq'] + '\tNA\tNA\tNA\t'
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\t'
            if currDict['type'] == 'V':
                fLine += seg + '\tNA\tNA\n'
            else:
                fLine += 'NA\t' + seg + '\tNA\n'
        #print fLine
            outF.write(str(fLine))


def find_curr_rank(segDict, seg, currLen):
    """
    Parameters
    ----------

    Returns
    -------

    """
    rank = 1
    for s in segDict:
        if s != seg:
            if segDict[s]['len'] > currLen:
                rank += 1
    return rank

def add_segments_to_dict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    head = tcrRecord.id
    headArr = head.split('.')
    vId = headArr[0]
    jId = headArr[1].split('(')[0]
    currSeqArr = tcrRecord.seq.split('N')
    vSeq = currSeqArr[0]
    jSeq = currSeqArr[-1]
    minLen = min(len(fastaDict[vId]),len(fastaDict[jId]))
    tupArr = [(vId, vSeq),(jId, jSeq)]
    for i in range(0,len(tupArr)):
        (id,seq) = tupArr[i]
        if id not in writtenArr:
            if id in segDict:
                if i == 0:
                    if idNameDict[jId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[jId])
                    if str(segDict[id]['seq'][-20:]) != str(seq[-20:]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same V-segment %s\n' % id)
                        sys.stderr.flush()
                else:
                    if idNameDict[vId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[vId])
                    if str(segDict[id]['seq'][:20]) != str(seq[:20]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same J-segment %s\n' % id)
                        sys.stderr.flush()
            else:
                segDict[id] = dict()
                segDict[id]['status'] = status
                segDict[id]['seq'] = seq
                segDict[id]['len'] = len(seq) - minLen
                segDict[id]['pairs'] = []

                if i == 0:
                    segDict[id]['type'] = 'V'
                    segDict[id]['pairs'].append(idNameDict[jId])
                else:
                    segDict[id]['type'] = 'J'
                    segDict[id]['pairs'].append(idNameDict[vId])
                segDict[id]['name'] = idNameDict[id]

    return segDict


def get_rank(tcr, rsemDict, unRsemDict, isProd, noRsem):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if isProd:
        currDict = rsemDict
    else:
        currDict = unRsemDict
    if not noRsem:
        currCount = currDict[tcr]
        rank = 1
        for rec in currDict:
            if rec != tcr:
                if unRsemDict[rec] > currCount:
                    rank += 1
        return rank
    else:
        return 'NA'


def get_un_dict_ratio(bamF, start, end, tcr, unDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    unMappedCount = 0
    usedArr = []
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in usedArr:
            usedArr.append(newName)
            if newName in unDict:
                unMappedCount += 1
    mappedFile.close()
    return (float(float(unMappedCount)/len(unDict)), unMappedCount)


def find_counts_in_region(bamF, start, end, tcr):
    """
    Parameters
    ----------

    Returns
    -------

    """
    readsArr = []
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in readsArr:
            readsArr.append(newName)
    mappedFile.close()
    counts = len(readsArr)
    return counts

def make_rsem_dict(rsemF, cdrDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    fDict = dict()
    unDict = dict()
    f = open(rsemF,'r')
    f.readline()
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        name = lArr[1]
        if name in cdrDict:
            if cdrDict[name]['stat'] == 'Productive':
                fDict[name] = float(lArr[4])
            unDict[name] = float(lArr[4])
        l = f.readline()
    f.close()
    return (fDict,unDict)

def find_CDR3(fasta, aaDict, vdjFaDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    f = open(fasta, 'rU')
    fDict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        if record.id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! sane name for two fasta entries %s\n' % record.id)
            sys.stderr.flush()
        else:
            idArr = record.id.split('.')
            vSeg = idArr[0]
            jSeg = idArr[1]
            if ((vSeg in aaDict) & (jSeg in aaDict)):
                currDict = find_V_and_J_aa_map(aaDict[vSeg],aaDict[jSeg],record.seq)
            else:
                if vSeg in aaDict:
                    newVseg = aaDict[vSeg]
                else:
                    vId = idArr[3]
                    currSeq = vdjFaDict[vId]
                    newVseg = get_best_V_aa(Seq(currSeq))
                if jSeg in aaDict:
                    newJseg = aaDict[jSeg]
                else:
                    jId = idArr[4]
                    currSeq= vdjFaDict[jId]
                    newJseg = get_best_J_aa(Seq(currSeq))
                currDict = find_V_and_J_aa_map(newVseg,newJseg,record.seq)
            fDict[record.id] = currDict
    f.close()
    return fDict


# either FXXG OR GXG pattern that J ends with
def get_best_J_aa(currSeq):
    """
    Parameters
    ----------

    Returns
    -------

    """
    firstSeq = get_NT_seq(currSeq)
    secondSeq = get_NT_seq(currSeq[1:])
    thirdSeq = get_NT_seq(currSeq[2:])
    pos = 10
    seq = ''
    found = False
    for s in [firstSeq, secondSeq, thirdSeq]:
        tempSeq = s[:8]
        indF = tempSeq.find('F')
        if indF != -1:
            if s[indF+3] == 'G':
                found = True
                if indF < pos:
                    pos = indF
                    seq = s[indF:]
        indG = tempSeq.find('G')
        if (indG != -1):
            if found == False:
                if s[indG+2] == 'G':
                    if indG < pos:
                        found = True
                        seq = s[indG:]
                        pos = indG
    if ((found == False) & (indF != -1)):
        seq = s[indF:]
    if seq != '':
        return seq
    else:
        return firstSeq

# if anything missing in data file, 
# go backwards till last C found
def get_best_V_aa(currSeq):
    """
    Parameters
    ----------

    Returns
    -------

    """
    firstSeq = get_NT_seq(currSeq)
    secondSeq = get_NT_seq(currSeq[1:])
    thirdSeq = get_NT_seq(currSeq[2:])
    pos = 10
    seq = ''
    for s in [firstSeq, secondSeq, thirdSeq]:
        #print "S: " + s
        tempSeq = s[-8:]
        #print "tempSeq: " + tempSeq
        ind = tempSeq.find('C')
        stopInd = tempSeq.find('*')
        #print "Ind: " + str(ind)
        if ((ind != -1) & (stopInd == -1)):
            #print "inside the ind"
            if ind < pos:
                goodRF = is_good_RF(s)
                if goodRF:
                    pos = ind
                    seq = s[:-8+ind + 1]
    if seq != '':
        return seq
    else:
        return firstSeq

# looking for something that starts with M, no * inside, ends with C
def is_good_RF(s):
    """
    Parameters
    ----------

    Returns
    -------

    """
    mInd = s.find('M')
    if mInd == -1:
        return False
    stopInd = s.find('*')
    if stopInd == -1:
        return True
    stopIndNext = s[stopInd+1:].find('*')
    while stopIndNext != -1:
        stopInd = stopIndNext + stopInd + 1
        stopIndNext = s[stopInd+1:].find('*')
        mInd = s[stopInd+1:].find('M')
        mInd = mInd + stopInd + 1
    if mInd != -1:
        return True
    else:
        return False



def find_V_and_J_aa_map(vSeg,jSeg,fullSeq):
    """
    Parameters
    ----------

    Returns
    -------

    """
    fDict = dict()
    # 3 possible frames
    firstSeq = get_NT_seq(fullSeq)
    secondSeq = get_NT_seq(fullSeq[1:])
    thirdSeq = get_NT_seq(fullSeq[2:])
    ntArr = [fullSeq, fullSeq[1:],fullSeq[2:]]
    aaSeqsArr = [firstSeq, secondSeq, thirdSeq]
    cdrArr = []
    posArr = []
    # look at each possible frame to see 
    # which has V and J in same frame
    for aaSeq in aaSeqsArr:
        (cdr, pos) = get_CDR3(aaSeq, vSeg,jSeg)
        cdrArr.append(cdr)
        posArr.append(pos)
    maxLen = 0
    bestCDR = ''
    bestSeq = ''
    hasStop = False
    bestPos = -1
    bestCDRnt = ''
    foundGood = False
    vPos = -1
    jPos = -1
    for i in range(0,3):
        if posArr[i] != -1:
            if ((cdrArr[i] != 'Only J') & (cdrArr[i] != 'Only V')):
                if len(cdrArr[i]) > maxLen:
                    if cdrArr[i].find('*') == -1:
                        foundGood = True
                        bestCDR = cdrArr[i]
                        bestPos = posArr[i]
                        maxLen = len(cdrArr[i])
                        bestSeq = ntArr[i]
                    else:
                        if maxLen == 0:
                            foundGood = True
                            bestPos = posArr[i]
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            hasStop = True
                else:
                    if hasStop == True:
                        # asterisk present means there's a stop codon
                        if cdrArr[i].find('*') == -1:
                            foundGood = True
                            bestPos = posArr[i]
                            hasStop = False
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
            else:
                if not foundGood:
                    if (cdrArr[i] == 'Only J'):
                        jPos = posArr[i]-i
                    elif (cdrArr[i] == 'Only V'):
                        vPos = posArr[i]-i
    if ((vPos != -1) & (jPos != -1) & (not foundGood)):
        bestCDRnt = fullSeq[3*vPos:3*jPos]
        bestCDR = 'NA'
    elif bestPos != -1:
        bestCDRnt = bestSeq[3*bestPos : 3*bestPos+3*len(bestCDR)]
    if bestCDR.find('*') != -1:
        stat = 'Unproductive - stop codon'
    else:
        stat = 'Productive'
    if maxLen == 0:
        if (('Only J' in cdrArr) & ('Only V' in cdrArr)):
            stat = 'Unproductive - Frame shift'
        else:
            if (('Only J' not in cdrArr) & ('Only V' in cdrArr)):
                stat = 'Unproductive - found only V segment'
            elif (('Only J' in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - found only J segment'
            elif (('Only J' not in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - didn\'t find V and J segment'
            else:
                stat = 'Unproductive'
            bestCDR = 'NA'
            bestCDRnt = 'NA'
    fDict['stat'] = stat
    fDict['CDR3 AA'] = bestCDR
    fDict['CDR3 NT'] = bestCDRnt
    fDict['Full Seq'] = fullSeq
    return fDict


# Is get_CDR3() essentially calling .find()?

# annotates the actual junction
# V and J have a conserved amino acid
# stuff in between those conserved amino acid 
# is CDR3
def get_CDR3(aaSeq, vSeq, jSeq):
    """
    Parameters
    ----------

    Returns
    -------

    """
    minDist = 14
    pos = -1
    for i in range(0,len(aaSeq) - len(vSeq) + 1):
        subAA = aaSeq[i:i+len(vSeq)]
        if len(subAA) != len(vSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong sub length\n')
            sys.stderr.flush()
        dist = 0
        for k in range(0,len(vSeq)):
            if vSeq[k] != subAA[k]:
                dist += 1
        if ((dist < minDist) & (subAA.endswith('C'))):
            minDist = dist
            pos = i + len(vSeq)
    jPos = -1
    minDistJ = 2
    for j in range(pos+1, len(aaSeq) - len(jSeq) + 1):
        subAA = aaSeq[j: j + len(jSeq)]
        if len(subAA) != len(jSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0,len(jSeq)):
            if jSeq[m] != subAA[m]:
                dist += 1
        if (dist <= minDistJ):
            if is_legal(subAA):
                jPos = j
                minDistJ = dist
    if pos == -1:
        if jPos != -1:
            return('Only J', jPos)
        else:
            return('No V/J found', -1)
    else:
        if jPos == -1:
            return('Only V',pos)
    return(aaSeq[pos:jPos], pos)


def is_legal(subAA):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if ((subAA[0] == 'F') & (subAA[3] == 'G')):
        return True
    if ((subAA[1] == 'G') & (subAA[3] == 'G')):
        return True
    return False


def get_NT_seq(fullSeq):
    """
    Parameters
    ----------

    Returns
    -------

    """
    mod = len(fullSeq) % 3
    if mod != 0:
        fSeq = fullSeq[:-mod].translate()
    else:
        fSeq = fullSeq.translate()
    return fSeq


def make_AA_Dict(aaF):
    """
    Parameters
    ----------

    Returns
    -------

    """
    fDict = dict()
    f = open(aaF,'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        if lArr[0] in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! %s appear twice in AA file\n' % lArr[0])
            sys.stderr.flush()
        fDict[lArr[0]] = lArr[1]
        l = f.readline()
    f.close()
    return fDict


def pick_final_isoforms(fullTcrFileAlpha, fullTcrFileBeta, output):
    """
    Parameters
    ----------

    Returns
    -------

    """
    pick_final_isoform_chain(fullTcrFileAlpha, output + '.alpha.full.TCRs.bestIso.fa', output + '.alpha.rsem.out.genes.results')
    pick_final_isoform_chain(fullTcrFileBeta, output + '.beta.full.TCRs.bestIso.fa', output + '.beta.rsem.out.genes.results')


def pick_final_isoform_chain(fullTCRfa, newFasta, rsemF):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if os.path.isfile(fullTCRfa):
        f = open(fullTCRfa, 'rU')
        outFa = open(newFasta, 'w')
        fastaDict = dict()
        byVJDict = dict()
        for record in SeqIO.parse(f,'fasta'):
            if record.id in fastaDict:
            #record.id = record.id + '_2'
                sys.stderr.write(str(datetime.datetime.now()) + 'error! same name for two fasta entries %s\n' % record.id)
                sys.stderr.flush()
            fastaDict[record.id] = record.seq
            onlyVJrec = str(record.id)
            idArr = onlyVJrec.strip('\n').split('.')
            vjStr = idArr[0] + '.' + idArr[1]
            if vjStr not in byVJDict:
                byVJDict[vjStr] = []
            byVJDict[vjStr].append(record.id)
        for vjStr in byVJDict:

            if len(byVJDict[vjStr]) == 1:
                cId = byVJDict[vjStr][0]
                cSeq = fastaDict[cId]
                newRec = SeqRecord(cSeq, id = cId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
            else:
            #print vjStr
            #print byVJDict[vjStr]
                bestId = find_best_C(byVJDict[vjStr], rsemF)
            #print "best: " + bestId
                bSeq = fastaDict[bestId]
                newRec = SeqRecord(bSeq, id = bestId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
        outFa.close()
        f.close()


def find_best_C(vjArr, rsemF):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if (os.path.exists(rsemF)):
        f = open(rsemF, 'r')
        f.readline()
        l = f.readline()
        bestSeq = 'name'
        maxCount = 0.0
        while l != '':
            lArr = l.strip('\n').split('\t')
            if lArr[0] in vjArr:
                currCount = float(lArr[4])
                if currCount > maxCount:
                    bestSeq = lArr[0]
                    maxCount = currCount
            l = f.readline()
        f.close()
        if bestSeq == 'name':
            return vjArr[0]
        return bestSeq
    else:
        return vjArr[0]


def run_rsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if samtools != '':
        if samtools[-1] != '/':
            rsem += '/'
    rsemIndDir = outDir + 'rsem_ind'
    if os.path.exists(rsemIndDir) == False:
        os.makedirs(rsemIndDir)
    if rsem != '':
        if rsem[-1] != '/':
            rsem += '/'
    if bowtie2 != '':
        if bowtie2[-1] != '/':
            bowtie2 += '/'
    if os.path.exists(fullTcrFileAlpha):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2,
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-path',bowtie2, '--bowtie2-mismatch-rate', '0.0' , output + '.alpha.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-mismatch-rate', '0.0', output + '.alpha.mapped.and.unmapped.fa', rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        unsortedBam = output + '.alpha.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for alpha chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.alpha.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])

    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any alpha chains, not running RSEM on alpha\n")
        sys.stdout.flush()
    if os.path.exists(fullTcrFileBeta):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2', '--bowtie2-path',
                             bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa', output + '.beta.R2.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2',
                              '--bowtie2-mismatch-rate', '0.0', output + '.beta.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        unsortedBam = output + '.beta.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for beta chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.beta.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])
    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any beta chains, not running RSEM on beta\n")
        sys.stdout.flush()


def create_TCR_full_output(fastaDict, tcr, outName, bases, mapDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    tcrF = open(tcr, 'rU')
    found = False
    ffound = False
    for tcrRecord in SeqIO.parse(tcrF, 'fasta'):
        tcrSeq = str(tcrRecord.seq)
        if tcrSeq.find('NNNNN') == -1:
            if ffound == False:
                ffound = True
                outF = open(outName, 'a+')
            idArr = tcrRecord.id.split('.')

            vEns = idArr[0]
            jEns = idArr[1].split('(')[0]
            vSeq = fastaDict[vEns]
            jSeq = fastaDict[jEns]
            vSeqTrim = ''
            jSeqTrim = ''
            print("**************", len(vSeq), len(jSeq), "**************")
            if bases == -10:
                bases = min(len(vSeq), len(jSeq))
            found = False
            for i in reversed(range(20,bases)):
                juncStart = tcrSeq[:i]
                vInd = vSeq.find(juncStart)
                if (vInd != -1):
                    found = True
                    vSeqTrim = vSeq[:vInd]
                    break
            if found == False:
                vSeqTrim = vSeq[:-bases]
            found = False
            for j in reversed(range(20,bases)):
                juncEnd = tcrSeq[-j:]
                jInd = jSeq.find(juncEnd)
                if (jInd != -1):
                    found = True
                    jSeqTrim = jSeq[jInd + j:]
                    break
            if found == False:
                jSeqTrim = jSeq[bases:]
            # Add TRBC or TRAC
            cArr = []
            if (str(tcrRecord.id).find('TRB')!= -1):
                for ens in mapDict:
                    if mapDict[ens].find('TRBC') != -1:
                        cArr.append(ens)
            elif (str(tcrRecord.id).find('TRA')!= -1):
                for ens in mapDict:
                    if mapDict[ens].find('TRAC') != -1:
                        cArr.append(ens)
            else:
                sys.stderr.write(str(datetime.datetime.now()) + " Error! no TRBC or TRAC\n")
                sys.stderr.flush()
            for ens in cArr:
                cSeq = fastaDict[ens]
                newSeq = vSeqTrim + tcrSeq + jSeqTrim + cSeq
                newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + mapDict[ens] + '.' + vEns + '.' + jEns + '.' + ens
                record = SeqRecord(Seq(newSeq,IUPAC.ambiguous_dna), id = newId, description = '')
                SeqIO.write(record,outF,'fasta')
    tcrF.close()
    if found == True:
        outF.close()


def add_mapped_pairs_to_seq_dict(seqDict, bam, out, lowQ, alignedDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    firstDict = dict()
    secondDict = dict()
    for name in seqDict:
        if seqDict[name][0] == '0':
            firstDict[name] = '1'
        if seqDict[name][1] == '1':
            secondDict[name] = '1'
        if ((seqDict[name][1] == '1') & (seqDict[name][0] == '0')):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! empty record insdie seqDict\n')
            sys.stderr.flush()
    f = pysam.AlignmentFile(bam,"rb")
    readsIter = f.fetch()
    for read in readsIter:
        name = read.query_name
        pos = -1
        if ((read.query_name in firstDict) & (read.is_read1)):
            pos = 0
            check = '0'
        elif ((read.query_name in secondDict) & (read.is_read2)):
            pos = 1
            check = '1'
        if pos != -1:
            qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
            if read.is_read1:
                currName = name + '\\1'
            else:
                currName = name + '\\2'
            if lowQ:
                if read.mate_is_reverse:
                    if read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                else:
                    if not read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                if currName not in alignedDict:
                    alignedDict[currName] = str(rSeq)
                    record = SeqRecord(rSeq, id = currName, description = '')
                    SeqIO.write(record,out,'fasta')
                else:
                    if str(alignedDict[currName]) != str(rSeq):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different sequences in alignedDict\n' % currName)
                        sys.stderr.flush()
            if read.is_reverse:
                qSeq = qSeq.reverse_complement()
            if seqDict[name][pos] != check:
                if str(qSeq) != str(seqDict[name][pos]):
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different mapped sequences not in the V/J region\n' % name)
                    sys.stderr.flush()
            seqDict[name][pos] = qSeq
    f.close()
    return seqDict

def write_unmapped_reads_to_dict_SE(unmapped):
    """
    Parameters
    ----------

    Returns
    -------

    """
    un_dict = {}
    f = pysam.AlignmentFile(unmapped,"rb")
    readsIter = f.fetch(until_eof = True)
    for read in readsIter:
        name = read.query_name
        un_dict[name] = '1'
    return un_dict


# Create a dict {'Alpha':{'C':[bed],'V':[bed],'J':[bed]}, 'Beta':{'C':[],'V':[],'J':[]}}
def make_VDJ_bed_dict(bed,idNameDict):
    """
    Parameters
    ----------

    Returns
    -------

    """
    fDict = {'Alpha':{'C':[],'V':[],'J':[]}, 'Beta':{'C':[],'V':[],'J':[]}}
    f = open(bed, 'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        gID = lArr[3]
        gName = idNameDict[gID]
        chain = ''
        if gName.startswith('TRA'):
            chain = 'Alpha'
        elif gName.startswith('TRB'):
            chain = 'Beta'
        else:
            print(str(datetime.datetime.now()) + ' Error! {} name is not alpha or beta chain, ignoring it'.format(gName))
        if gName.find('C') != -1:
            fDict[chain]['C'].append(l)
        elif gName.find('V') != -1:
            fDict[chain]['V'].append(l)
        elif gName.find('J') != -1:
            fDict[chain]['J'].append(l)
        l = f.readline()
    f.close()
    return fDict


# Creates a dictionary of ENSEMBL ID -> fasta sequence
def make_fasta_dict(fasta):
    """
    Parameters
    ----------

    Returns
    -------

    """
    inF = open(fasta,'rU')
    fastaDict = dict()
    for record in SeqIO.parse(inF, 'fasta'):
        fastaDict[record.id] = str(record.seq)
    inF.close()
    return fastaDict


# Creates a dictionary of ENSEMBL ID -> Gene name
def make_id_name_dict(mapping):
    """
    Parameters
    ----------

    Returns
    -------

    """
    f = open(mapping, 'r')
    fDict = dict()
    linesArr = f.read().split('\n')
    f.close()
    for line in linesArr:
        lineArr = line.split('\t')
        id = lineArr[0]
        name = lineArr[1]
        ind = name.find('Gene:')
        if ind != -1:
            name = name[ind+len('Gene:'):]
        if id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s appear twice in mapping file\n' % id)
            sys.stderr.flush()
        fDict[id] = name
    return fDict


def check_parameters(genome, strand, single_cell, path, sumF):
    """
    Parameters
    ----------

    Returns
    -------

    """
    if genome != 'hg38' and genome != 'mm10' and genome != 'hg19' and genome != 'mm10_ncbi':
        sys.exit("-genome only accept one of the following: mm10, mm10_ncbi, hg38, hg19")
    if strand.lower() not in ['none','minus','plus']:
        sys.exit("-strand should be one of: none, minus, plus")
    if not single_cell:
        if path == '':
            sys.exit("when running on multiple cells you must include the -path parameter")
        if sumF == '':
            sys.exit("when running on multiple cells you must include the -sumF parameter")
        if not os.path.isdir(path):
            sys.exit("%s path does not exists. Please check your -path parameter and run again" % path)

def main(args):
    run_TCR_pipe(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome','-g','-G', help='Alignment genome. Currently supported: mm10 and hg38', required=True)
    parser.add_argument('-single_cell', help='add if you are only running on a single cell. If so,'
                                                        'it will ignore -path and -subpath arguments', action='store_true')
    parser.add_argument('-lowQ', help='add if you want to add \"low quality\" reads as input to the reconstruction '
                                                        'algorithm', action='store_true')
    parser.add_argument('-path','-p','-P', help='The path for the data directory. Assumes that every subdirectory'
                                                    'is a single cell', default='')
    parser.add_argument('-sumF', help='prefix for summary outputs', default='')
    parser.add_argument('-bowtie2','-bw','-BW', help='Path to bowtie2. If not used assumes that bowtie2 is in the'
                                                 'default path', default = '')
    parser.add_argument('-rsem','-RSEM', help='Path to rsem. If not used assumes that rsem is in the'
                                                'default path', default = '')
    parser.add_argument('-strand', help='Strand of the right most read in genomic coordinates. Options are: [minus, plus, '
                                        'none]. Defualt is minus', default = 'minus')
    parser.add_argument('-output','-out','-o','-O', help='output prefix, relative to /path/singleCellFolder', required=True)
    parser.add_argument('-bam', help='Input bam alignment file, relative to /path/singleCellFolder/ if working on multiple files', default = './tophat_output/picard_output/sorted.bam')
    parser.add_argument('-unmapped','-u','-U', help='bam file of the unmapped reads, relative to /path/singleCellFolder/', default = './tophat_output/unmapped.bam')
    parser.add_argument('-bases','-b','-B', help='Number of bases to take from each V and J segments, default is min(len(V), len(J) ', type=int, default=-10)
    parser.add_argument('-iterations','-iter','-i','-I', help='Number of iterations for the reconstruction'
                                                              'algorithm, default is 20', type=int, default=20)
    parser.add_argument('-samtools', help='Path to samtools. If not used assumes that samtools is in the default path', default = '')
    parser.add_argument('-score','-sc','-SC', help='Alignment score threshold. Default is 15', type=int, default=15)
    parser.add_argument('-overlap','-ol','-OL', help='Number of minimum bases that overlaps V and J ends,'
                                                              'default is 10', type=int, default=10)
    parser.add_argument('-trim', help='The number of bases to trim from each end of a sigle-end read.', type=int, default=25)
    args = parser.parse_args()
    main(args)
