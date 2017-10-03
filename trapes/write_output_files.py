"""
This module contains functions used to write output files containing summary statistics for individual cells and
summary statistics for all cells.

"""

import os
from Bio import SeqIO


def make_single_cell_output_file(alpha_dict, beta_dict, output, beta_rsem, alpha_rsem, alpha_bam, beta_bam, fasta_dict,
                                 un_dict_alpha, un_dict_beta, id_name_dict):
    out_f = open(output + '.summary.txt', 'w')
    out_f.write(
        'Chain\tStatus\tRank of TCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in TCR\t#reads in CDR3\t#reads in V\t#reads in '
        'J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the '
        'reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n')
    if (len(alpha_dict) > 0):
        write_chain(out_f, 'alpha', alpha_dict, alpha_rsem, alpha_bam, fasta_dict, un_dict_alpha, output, id_name_dict)
    if (len(beta_dict) > 0):
        write_chain(out_f, 'beta', beta_dict, beta_rsem, beta_bam, fasta_dict, un_dict_beta, output, id_name_dict)
    out_f.close()


def write_chain(out_f, chain, cdr_dict, rsem_f, bam_f, fasta_dict, un_dict, output, id_name_dict):
    written_arr = []
    if os.path.exists(rsem_f):
        no_rsem = False
        (rsem_dict, un_rsem_dict) = make_rsem_dict(rsem_f, cdr_dict)
    else:
        no_rsem = True
    for tcr in cdr_dict:
        j_start = -1
        cdr_ind = -1
        c_ind = -1
        if cdr_dict[tcr]['stat'] == 'Productive':
            is_prod = True
        else:
            is_prod = False
        f_line = chain + '\t' + cdr_dict[tcr]['stat'] + '\t'
        if no_rsem:
            rank = 'NA'
        else:
            rank = get_rank(tcr, rsem_dict, un_rsem_dict, is_prod, no_rsem)
        f_line += str(rank) + '\t'
        name_arr = tcr.split('.')
        f_line += name_arr[0] + '\t' + name_arr[1] + '\t' + name_arr[2] + '\t'
        f_line += cdr_dict[tcr]['CDR3 NT'] + '\t' + cdr_dict[tcr]['CDR3 AA'] + '\t'
        full_seq = cdr_dict[tcr]['Full Seq'].upper()
        if not no_rsem:
            total_count = find_counts_in_region(bam_f, 0, len(full_seq), tcr)
            f_line += str(total_count) + '\t'
            c_name = name_arr[5]
            while c_name.endswith('_2'):
                c_name = c_name[:-2]
            c_seq = fasta_dict[c_name].upper()
            c_ind = full_seq.find(c_seq)
            if c_ind == -1:
                sys.stderr.write(
                    str(datetime.datetime.now()) + 'Error! could not find C segment sequence in the full sequence\n')
                sys.stderr.flush()
                c_counts = 'NA'
            else:
                c_counts = find_counts_in_region(bam_f, c_ind, len(full_seq), tcr)
            if cdr_dict[tcr]['CDR3 NT'] != 'NA':
                cdr_ind = full_seq.find(cdr_dict[tcr]['CDR3 NT'].upper())
            else:
                cdr_ind = -1
            if ((cdr_ind == -1) & (cdr_dict[tcr]['CDR3 NT'] != 'NA')):
                sys.stderr.write(
                    str(datetime.datetime.now()) + ' Error! Cound not find CDR3 NT sequence in the full sequence\n')
                sys.stderr.flush()
            if cdr_ind != -1:
                cdr_counts = find_counts_in_region(bam_f, cdr_ind, cdr_ind + len(cdr_dict[tcr]['CDR3 NT']), tcr)
                j_start = cdr_ind + len(cdr_dict[tcr]['CDR3 NT'])
                if c_ind != -1:
                    j_counts = find_counts_in_region(bam_f, j_start, c_ind, tcr)
                else:
                    j_counts = find_counts_in_region(bam_f, j_start, j_start + 50, tcr)
                v_counts = find_counts_in_region(bam_f, 0, cdr_ind, tcr)
                f_line += str(cdr_counts) + '\t' + str(v_counts) + '\t' + str(j_counts) + '\t' + str(c_counts) + '\t'
            else:
                f_line += 'NA\tNA\tNA\t' + str(c_counts) + '\t'
            v_id = name_arr[3]
            j_id = name_arr[4]
            c_id = name_arr[5]
            if cdr_dict[tcr]['CDR3 NT'] != 'NA':
                (un_dict_ratio_cdr, un_cdr_count) = get_un_dict_ratio(bam_f, cdr_ind,
                                                                      cdr_ind + len(cdr_dict[tcr]['CDR3 NT']), tcr,
                                                                      un_dict)
                (un_dict_ratio_all, un_all_count) = get_un_dict_ratio(bam_f, 0, len(full_seq), tcr, un_dict)
                f_line += str(un_dict_ratio_all) + '\t' + str(un_all_count) + '\t' + str(
                    un_dict_ratio_cdr) + '\t' + str(
                    un_cdr_count) + '\t'
            else:
                f_line += 'NA\tNA\tNA\tNA\t'
            written_arr.append(v_id)
            written_arr.append(j_id)
            f_line += v_id + '\t' + j_id + '\t' + c_id + '\n'
        else:
            f_line += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' + name_arr[3] + '\t' + name_arr[4] + '\t' + \
                      name_arr[
                          5] + '\n'

            # print fLine

        out_f.write(str(f_line))
    write_failed_reconstructions(out_f, chain, written_arr, output, id_name_dict, fasta_dict)


def make_rsem_dict(rsem_f, cdr_dict):
    f_dict = dict()
    un_dict = dict()
    f = open(rsem_f, 'r')
    f.readline()
    l = f.readline()
    while l != '':
        l_arr = l.strip('\n').split('\t')
        name = l_arr[1]
        if name in cdr_dict:
            if cdr_dict[name]['stat'] == 'Productive':
                f_dict[name] = float(l_arr[4])
            un_dict[name] = float(l_arr[4])
        l = f.readline()
    f.close()
    return (f_dict, un_dict)


def get_rank(tcr, rsem_dict, un_rsem_dict, is_prod, no_rsem):
    if is_prod:
        curr_dict = rsem_dict
    else:
        curr_dict = un_rsem_dict
    if not no_rsem:
        curr_count = curr_dict[tcr]
        rank = 1
        for rec in curr_dict:
            if rec != tcr:
                if un_rsem_dict[rec] > curr_count:
                    rank += 1
        return rank
    else:
        return 'NA'


def find_counts_in_region(bam_f, start, end, tcr):
    reads_arr = []
    mapped_file = pysam.AlignmentFile(bam_f, "rb")
    reads_iter = mapped_file.fetch(tcr, start, end)
    for read in reads_iter:
        if read.is_read1:
            new_name = read.query_name + '_1'
        else:
            new_name = read.query_name + '_2'
        if new_name not in reads_arr:
            reads_arr.append(new_name)
    mapped_file.close()
    counts = len(reads_arr)
    return counts


def get_un_dict_ratio(bam_f, start, end, tcr, un_dict):
    un_mapped_count = 0
    used_arr = []
    mapped_file = pysam.AlignmentFile(bam_f, "rb")
    reads_iter = mapped_file.fetch(tcr, start, end)
    for read in reads_iter:
        if read.is_read1:
            new_name = read.query_name + '_1'
        else:
            new_name = read.query_name + '_2'
        if new_name not in used_arr:
            used_arr.append(new_name)
            if new_name in un_dict:
                un_mapped_count += 1
    mapped_file.close()
    return (float(float(un_mapped_count) / len(un_dict)), un_mapped_count)


def write_failed_reconstructions(out_f, chain, written_arr, output, id_name_dict, fasta_dict):
    rec_f = output + '.reconstructed.junctions.' + chain + '.fa'
    if os.path.isfile(rec_f):
        f = open(rec_f, 'rU')
        seg_dict = dict()
        for tcr_record in SeqIO.parse(f, 'fasta'):
            tcr_seq = str(tcr_record.seq)
            if tcr_seq.find('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN') != -1:
                status = 'Failed reconstruction - reached maximum number of iterations'
                seg_dict = add_segments_to_dict(seg_dict, status, written_arr, tcr_record, id_name_dict, fasta_dict)
            elif tcr_seq.find('NNNN') != -1:
                status = 'Failed reconstruction - V and J segment do not overlap'
                seg_dict = add_segments_to_dict(seg_dict, status, written_arr, tcr_record, id_name_dict, fasta_dict)
        f.close()
        if len(seg_dict) > 0:
            write_seg_dict(seg_dict, out_f, chain)


def add_segments_to_dict(seg_dict, status, written_arr, tcr_record, id_name_dict, fasta_dict):
    head = tcr_record.id
    head_arr = head.split('.')
    v_id = head_arr[0]
    j_id = head_arr[1].split('(')[0]
    curr_seq_arr = tcr_record.seq.split('N')
    v_seq = curr_seq_arr[0]
    j_seq = curr_seq_arr[-1]
    min_len = min(len(fasta_dict[v_id]), len(fasta_dict[j_id]))
    tup_arr = [(v_id, v_seq), (j_id, j_seq)]
    for i in range(0, len(tup_arr)):
        (id, seq) = tup_arr[i]
        if id not in written_arr:
            if id in seg_dict:
                if i == 0:
                    if id_name_dict[j_id] not in seg_dict[id]['pairs']:
                        seg_dict[id]['pairs'].append(id_name_dict[j_id])
                    if str(seg_dict[id]['seq'][-20:]) != str(seq[-20:]):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same '
                                                       'V-segment %s\n' % id)
                        sys.stderr.flush()
                else:
                    if id_name_dict[v_id] not in seg_dict[id]['pairs']:
                        seg_dict[id]['pairs'].append(id_name_dict[v_id])
                    if str(seg_dict[id]['seq'][:20]) != str(seq[:20]):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same '
                                                       'J-segment %s\n' % id)
                        sys.stderr.flush()
            else:
                seg_dict[id] = dict()
                seg_dict[id]['status'] = status
                seg_dict[id]['seq'] = seq
                seg_dict[id]['len'] = len(seq) - min_len
                seg_dict[id]['pairs'] = []

                if i == 0:
                    seg_dict[id]['type'] = 'V'
                    seg_dict[id]['pairs'].append(id_name_dict[j_id])
                else:
                    seg_dict[id]['type'] = 'J'
                    seg_dict[id]['pairs'].append(id_name_dict[v_id])
                seg_dict[id]['name'] = id_name_dict[id]

    return seg_dict


def write_seg_dict(seg_dict, out_f, chain):
    for seg in seg_dict:
        curr_dict = seg_dict[seg]
        pairs = ''
        for pair in curr_dict['pairs']:
            pairs += pair + '.'
        pairs = pairs[:-1]
        if curr_dict['len'] > 0:
            f_line = chain + '\t' + curr_dict['status'] + '\t'
            rank = find_curr_rank(seg_dict, seg, curr_dict['len'])
            f_line += str(rank) + '\t'
            if curr_dict['type'] == 'V':
                f_line += curr_dict['name'] + '\t' + 'paired with: ' + pairs + '\t'
            else:
                f_line += 'paired with: ' + pairs + '\t' + curr_dict['name'] + '\t'
            f_line += 'NA\t' + curr_dict['seq'] + '\tNA\tNA\tNA\t'
            f_line += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\t'
            if curr_dict['type'] == 'V':
                f_line += seg + '\tNA\tNA\n'
            else:
                f_line += 'NA\t' + seg + '\tNA\n'
                # print fLine
            out_f.write(str(f_line))


def find_curr_rank(seg_dict, seg, curr_len):
    rank = 1
    for s in seg_dict:
        if s != seg:
            if seg_dict[s]['len'] > curr_len:
                rank += 1
    return rank


def add_cell_to_tcr_sum(cell_folder, noutput, opened, tcr_f_out):
    if os.path.isfile(noutput + '.summary.txt'):
        curr_out = open(noutput + '.summary.txt', 'r')
        if not opened:
            opened = True
            head = curr_out.readline()
            head = 'cell\t' + head
            tcr_f_out.write(head)
        else:
            curr_out.readline()
        l = curr_out.readline()
        while l != '':
            new_l = cell_folder + '\t' + l
            tcr_f_out.write(new_l)
            l = curr_out.readline()
        curr_out.close()
    return opened


def add_to_stat_dict(noutput, cell_folder, final_stat_dict):
    if cell_folder in final_stat_dict:
        print "Error! %s appear more than once in final stat dictionary" % cell_folder
    final_stat_dict[cell_folder] = {'alpha': 'Failed - found V and J segments but wasn\'t able to extend them',
                                    'beta': 'Failed - found V and J segments but wasn\'t able to extend them'}
    if os.path.isfile(noutput + '.summary.txt'):
        curr_out = open(noutput + '.summary.txt', 'r')
        msg_a = 'None'
        msg_b = 'None'
        curr_out.readline()
        l = curr_out.readline()
        while l != '':
            l_arr = l.strip('\n').split('\t')
            chain = l_arr[0]
            stat = l_arr[1]
            if stat == 'Productive':
                if chain == 'alpha':
                    msg_a = 'Productive'
                else:
                    msg_b = 'Productive'
            elif stat == 'Productive (no 118 PHE found)':
                if chain == 'alpha':
                    msg_a = 'Productive (no 118 PHE found)'
                else:
                    msg_b = 'Productive (no 118 PHE found)'
            elif stat.startswith('Unproductive'):
                if chain == 'alpha':
                    if msg_a != 'Productive':
                        msg_a = 'Unproductive'
                else:
                    if msg_b != 'Productive':
                        msg_b = 'Unproductive'
            elif stat.startswith('Failed reconstruction'):
                if stat == 'Failed reconstruction - reached maximum number of iterations':
                    if chain == 'alpha':
                        if msg_a == 'None':
                            msg_a = 'Failed - reconstruction didn\'t converge'
                    else:
                        if msg_b == 'None':
                            msg_b = 'Failed - reconstruction didn\'t converge'
                elif stat == 'Failed reconstruction - V and J segment do not overlap':
                    if chain == 'alpha':
                        if msg_a == 'None':
                            msg_a = 'Failed - V and J reconstruction don\'t overlap'
                    else:
                        if msg_b == 'None':
                            msg_b = 'Failed - V and J reconstruction don\'t overlap'
            l = curr_out.readline()
        curr_out.close()
        if msg_a == 'None':
            alpha_junc = noutput + '.alpha.junctions.txt'
            if (os.path.isfile(alpha_junc) == True):
                if os.stat(alpha_junc).st_size == 0:
                    msg_a = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msg_a = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msg_a = 'Failed - didn\'t find any V and J segments in original mapping'
        if msg_b == 'None':
            beta_junc = noutput + '.beta.junctions.txt'
            if (os.path.isfile(beta_junc) == True):
                if os.stat(beta_junc).st_size == 0:
                    msg_b = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msg_b = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msg_b = 'Failed - didn\'t find any V and J segments in original mapping'

    else:
        beta_junc = noutput + '.beta.junctions.txt'
        alpha_junc = noutput + '.alpha.junctions.txt'
        if os.path.isfile(beta_junc) == True:
            if os.stat(beta_junc).st_size == 0:
                msg_b = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msg_b = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(alpha_junc) == True):
            if os.stat(alpha_junc).st_size == 0:
                msg_a = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msg_a = 'Failed - didn\'t find any V and J segments in original mapping'
    final_stat_dict[cell_folder]['alpha'] = msg_a
    final_stat_dict[cell_folder]['beta'] = msg_b
    return final_stat_dict
