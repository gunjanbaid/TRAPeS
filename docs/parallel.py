from joblib import Parallel, delayed
import multiprocessing

def runTCRpipe(genome, output, bam, unmapped, bases, strand, numIterations,thresholdScore, minOverlap, rsem, bowtie2, singleCell, path, sumF, lowQ, samtools, trim):
    checkParameters(genome, strand, singleCell, path, sumF)
    if singleCell == True:
        # TODO: Fix this, won't work for SE
        #runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore, minOverlap,
        #          rsem, bowtie2, lowQ, singleEnd, fastq, trimmomatic, transInd)
        sys.exit(0)
    if path == './':
        path = os.getcwd()
    if not path.endswith('/'):
        path = path + '/'
    finalStatDict = dict()
    tcrFout = open(sumF + '.TCRs.txt','w')
    opened = False
 
    def processCellWrapper(cellFolder):
        return processCell(genome, output, bam, unmapped, bases, strand, numIterations,
        thresholdScore, minOverlap, rsem, bowtie2, singleCell, path, sumF, lowQ, samtools, 
        top, by_exp, read_overlap, one_side, trim, cellFolder)
 
    # parallel processing of each cell 
    pool = ThreadPool(processes=20)
    pool.map(processCellWrapper, (cellFolder for cellFolder in os.listdir(path)))
    pool.close()

    for cellFolder in os.listdir(path):
    	opened = addCellToTCRsum(cellFolder, noutput, opened, tcrFout)
    	finalStatDict = addToStatDict(noutput, cellFolder, finalStatDict)
    sumFout = open(sumF + '.summary.txt','w')
    sumFout.write('sample\talpha\tbeta\n')
    for cell in sorted(finalStatDict):
        fout = cell + '\t' + finalStatDict[cell]['alpha'] + '\t' + finalStatDict[cell]['beta'] + '\n'
        sumFout.write(fout)
    sumFout.close()


def process_cell(genome, output, bam, unmapped, bases, strand, num_iterations, threshold_score, 
    min_overlap, rsem, bowtie2, single_cell, path, sum_f, low_q, samtools, top, by_exp, 
    read_overlap, one_side, trim, cellFolder):
    for cell_folder in os.listdir(path):
        full_path = path + cell_folder + '/'
        if ((os.path.exists(full_path)) & (os.path.isdir(full_path))):
            sys.stdout.write(str(datetime.datetime.now()) + " Working on: " + cell_folder + '\n')
            sys.stdout.flush()
            (found, nbam, nunmapped, noutput) = utils.format_files(full_path, bam, unmapped, output)
            if not found:
                sys.stderr.write(str(datetime.datetime.now()) + 
                " There is not a bam or unmapped file in this folder, moving to the next folder\n")
                sys.stderr.flush()
            else:
                curr_folder = os.path.abspath(os.path.dirname(sys.argv[0])) + '/'
                reconstruction = curr_folder + '/vdj.alignment'
                if genome == 'hg38':
                    fasta = curr_folder + 'data/hg38/hg38.TCR.fa'
                    bed = curr_folder + 'data/hg38/hg38.TCR.bed'
                    mapping = curr_folder + 'data/hg38/hg38.id.name.mapping.TCR.txt'
                    aa_f = curr_folder + 'data/hg38/hg38.TCR.conserved.AA.txt'
                    ref_ind = curr_folder + 'data/hg38/index/hg38'
                if genome == 'mm10':
                    fasta = curr_folder + 'data/mm10/mm10.TCR.fa'
                    bed = curr_folder + 'data/mm10/mm10.TCR.bed'
                    mapping = curr_folder + 'data/mm10/mm10.gene.id.mapping.TCR.txt'
                    aa_f = curr_folder + 'data/mm10/mm10.conserved.AA.txt'
                    ref_ind = curr_folder + 'data/mm10/index/mm10'
                if genome == 'mm10_ncbi':
                    fasta = curr_folder + 'data/mm10_ncbi/mm10.TCR.fa'
                    bed = curr_folder + 'data/mm10_ncbi/mm10.TCR.bed'
                    mapping = curr_folder + 'data/mm10_ncbi/mm10.gene.id.mapping.TCR.txt'
                    aa_f = curr_folder + 'data/mm10_ncbi/mm10.conserved.AA.txt'
                    ref_ind = curr_folder + 'data/mm10_ncbi/index/mm10'
                if genome == 'hg19':
                    fasta = curr_folder + 'data/hg19/hg19.TCR.fa'
                    bed = curr_folder + 'data/hg19/hg19.TCR.bed'
                    mapping = curr_folder + 'data/hg19/hg19.gene.id.mapping.TCR.txt'
                    aa_f = curr_folder + 'data/hg19/hg19.conserved.AA.txt'
                    ref_ind = curr_folder + 'data/hg19/index/hg19'

                process_single_cell.run_single_cell(fasta, bed, noutput, nbam, nunmapped, mapping,
                    bases, strand, reconstruction, aa_f, num_iterations, threshold_score,
                    min_overlap, rsem, bowtie2, low_q, samtools, top, by_exp,
                    read_overlap, one_side, ref_ind, trim)
                opened = write_output_files.add_cell_to_tcr_sum(cell_folder, noutput, opened, tcr_fout)
                final_stat_dict = write_output_files.add_to_stat_dict(noutput, cell_folder, final_stat_dict)