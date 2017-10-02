"""
This module contains functions used to process each individual cell.

"""

import sys
import os
import subprocess
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import paired_end
import write_output_files
import utils


def run_single_cell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, 
    aa_f, num_iterations, threshold_score, min_overlap, rsem, bowtie2, low_q, samtools, top, 
    by_exp, read_overlap, one_side, ref_ind, trim):
    id_name_dict = make_id_name_dict(mapping)
    fasta_dict = make_fasta_dict(fasta)
    vdj_dict = make_vdj_bed_dict(bed, id_name_dict)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing alpha chain\n")
    sys.stdout.flush()
    if utils.is_paired_end(bam):
        un_dict_alpha = paired_end.analyze_chain(fasta_dict, vdj_dict, output, bam, unmapped, id_name_dict, bases, 'A',strand, low_q, top, by_exp, read_overlap)
        sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing beta chain\n")
        sys.stdout.flush()
        un_dict_beta = paired_end.analyze_chain(fasta_dict, vdj_dict, output, bam, unmapped, id_name_dict, bases, 'B',strand, low_q, top, by_exp, read_overlap)
    else:
        print(str(datetime.datetime.now()) + " Pre-processing alpha chain")
        print(str(datetime.datetime.now()) + " Pre-processing beta chain")
        single_end.analyze_chain_single_end(fasta_dict, vdj_dict, output, bam, unmapped, 
            id_name_dict, bases, strand, lowQ, bowtie2, ref_ind, trim)    
        un_dict_alpha = single_end.write_unmapped_reads_to_dict_se(unmapped)
        un_dict_beta = dict(un_dict_beta)

    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing beta chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.beta.mapped.and.unmapped.fa', output + '.beta.junctions.txt', output + '.reconstructed.junctions.beta.fa', str(num_iterations), str(threshold_score), str(min_overlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing alpha chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.alpha.mapped.and.unmapped.fa', output + '.alpha.junctions.txt', output + '.reconstructed.junctions.alpha.fa', str(num_iterations), str(threshold_score), str(min_overlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Creating full TCR sequencing\n")
    sys.stdout.flush()
    full_tcr_file_alpha = output + '.alpha.full.TCRs.fa'
    tcr_f = output + '.reconstructed.junctions.alpha.fa'
    (c_seq, c_name, c_id) = utils.get_c_info(vdj_dict['Alpha']['C'][0], id_name_dict, fasta_dict)
    create_tcr_full_output(fasta_dict, tcr_f, full_tcr_file_alpha, bases, id_name_dict, c_seq, c_name, c_id, one_side)
    full_tcr_file_beta = output + '.beta.full.TCRs.fa'
    tcr_f = output + '.reconstructed.junctions.beta.fa'
    (c_seq, c_name, c_id) = utils.get_c_info(vdj_dict['Beta']['C'][0], id_name_dict, fasta_dict)
    create_tcr_full_output(fasta_dict, tcr_f, full_tcr_file_beta, bases, id_name_dict, c_seq, c_name, c_id, one_side)
    sys.stdout.write(str(datetime.datetime.now()) + " Running RSEM to quantify expression of all possible isoforms\n")
    sys.stdout.flush()
    out_dir_ind = output.rfind('/')
    if out_dir_ind != -1:
        out_dir = output[:out_dir_ind + 1]
    else:
        out_dir = os.getcwd()
    if paired_end:
        paired_end.run_rsem(out_dir, rsem, bowtie2, full_tcr_file_alpha, full_tcr_file_beta, output, samtools)
    else:
        single_end.runRsemSE(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools)
    pick_final_isoforms(full_tcr_file_alpha, full_tcr_file_beta, output)
    best_alpha = output + '.alpha.full.TCRs.bestIso.fa'
    best_beta = output + '.beta.full.TCRs.bestIso.fa'
    sys.stdout.write(str(datetime.datetime.now()) + " Finding productive CDR3\n")
    sys.stdout.flush()
    aa_dict = make_aa_dict(aa_f)
    if os.path.isfile(best_alpha):
        f_dict_alpha = find_cdr3(best_alpha, aa_dict, fasta_dict)
    else:
        f_dict_alpha = dict()
    if os.path.isfile(best_beta):
        f_dict_beta = find_cdr3(best_beta, aa_dict, fasta_dict)
    else:
        f_dict_beta = dict()
    beta_rsem_out = output + '.beta.rsem.out.genes.results'
    alpha_rsem_out = output + '.alpha.rsem.out.genes.results'
    alpha_bam = output + '.alpha.rsem.out.transcript.sorted.bam'
    beta_bam = output + '.beta.rsem.out.transcript.sorted.bam'
    sys.stdout.write(str(datetime.datetime.now()) + " Writing results to summary file\n")
    sys.stdout.flush()
    write_output_files.make_single_cell_output_file(f_dict_alpha, f_dict_beta, output, 
        beta_rsem_out, alpha_rsem_out, alpha_bam, beta_bam, fasta_dict,
        un_dict_alpha, un_dict_beta, id_name_dict)


# Creates a dictionary of ENSEMBL ID -> Gene name
def make_id_name_dict(mapping):
    f = open(mapping, 'r')
    f_dict = dict()
    lines_arr = f.read().split('\n')
    f.close()
    for line in lines_arr:
        line_arr = line.split('\t')
        id = line_arr[0]
        name = line_arr[1]
        ind = name.find('Gene:')
        if ind != -1:
            name = name[ind + len('Gene:'):]
        if id in f_dict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s appear twice in mapping file\n' % id)
            sys.stderr.flush()
        f_dict[id] = name
    return f_dict


# Creates a dictionary of ENSEMBL ID -> fasta sequence
def make_fasta_dict(fasta):
    inF = open(fasta, 'rU')
    fasta_dict = dict()
    for record in SeqIO.parse(inF, 'fasta'):
        fasta_dict[record.id] = str(record.seq)
    inF.close()
    return fasta_dict


# Create a dict {'Alpha':{'C':[bed],'V':[bed],'J':[bed]}, 'Beta':{'C':[],'V':[],'J':[]}}
def make_vdj_bed_dict(bed, id_name_dict):
    f_dict = {'Alpha': {'C': [], 'V': [], 'J': []}, 'Beta': {'C': [], 'V': [], 'J': []}}
    f = open(bed, 'r')
    l = f.readline()
    while l != '':
        l_arr = l.strip('\n').split('\t')
        g_id = l_arr[3]
        g_name = id_name_dict[g_id]
        chain = ''
        if (g_name.startswith('TRA')):
            chain = 'Alpha'
        elif (g_name.startswith('TRB')):
            chain = 'Beta'
        else:
            sys.stderr.write(
                str(datetime.datetime.now()) + ' Error! %s name is not alpha or beta chain, ignoring it\n' % g_name)
            sys.stderr.flush()
        if g_name.find('C') != -1:
            f_dict[chain]['C'].append(l)
        elif g_name.find('V') != -1:
            f_dict[chain]['V'].append(l)
        elif g_name.find('J') != -1:
            f_dict[chain]['J'].append(l)
        l = f.readline()
    f.close()
    return f_dict


def make_aa_dict(aa_f):
    f_dict = dict()
    f = open(aa_f, 'r')
    l = f.readline()
    while l != '':
        l_arr = l.strip('\n').split('\t')
        if l_arr[0] in f_dict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! %s appear twice in AA file\n' % l_arr[0])
            sys.stderr.flush()
        f_dict[l_arr[0]] = l_arr[1]
        l = f.readline()
    f.close()
    return f_dict


def create_tcr_full_output(fasta_dict, tcr, out_name, bases, map_dict, c_seq, c_name, c_id, one_side):
    tcr_f = open(tcr, 'rU')
    found = False
    f_found = False
    rec_name_arr = []
    for tcr_record in SeqIO.parse(tcr_f, 'fasta'):
        added_c = False
        tcr_seq = str(tcr_record.seq)
        if tcr_seq.find('NNNNN') == -1:
            if f_found == False:
                f_found = True
                out_f = open(out_name, 'w')
            id_arr = tcr_record.id.split('.')
            v_ens = id_arr[0]
            j_ens = id_arr[1].split('(')[0]
            v_seq = fasta_dict[v_ens]
            j_seq = fasta_dict[j_ens]
            rec_name_arr = write_record(tcr_record, tcr_seq, added_c, v_ens, j_ens, v_seq, j_seq, map_dict, bases,
                                        c_seq, c_id,
                                        c_name, out_f, fasta_dict, rec_name_arr)
        elif one_side:
            cur_seq = tcr_seq.split('NNNN')[0]
            j_seg = find_best_j_for_seq(cur_seq, fasta_dict, map_dict)
            if j_seg != 'NA':
                if f_found == False:
                    f_found = True
                    out_f = open(out_name, 'w')
                id_arr = tcr_record.id.split('.')
                v_ens = id_arr[0]
                v_seq = fasta_dict[v_ens]
                for j_ens in j_seg:
                    j_seq = fasta_dict[j_ens]
                    rec_name_arr = write_record(tcr_record, cur_seq, added_c, v_ens, j_ens, v_seq, j_seq, map_dict,
                                                bases, c_seq,
                                                c_id, c_name, out_f, fasta_dict, rec_name_arr)
    tcr_f.close()
    if found == True:
        out_f.close()


def find_best_j_for_seq(cur_seq, fasta_dict, id_name_dict):
    j_arr_old = find_js_per_len(cur_seq, fasta_dict, id_name_dict, 20)
    if len(j_arr_old) == 0:
        return 'NA'
    for x in range(21, len(cur_seq)):
        new_arr = find_js_per_len(cur_seq, fasta_dict, id_name_dict, x)
        if len(new_arr) == 0:
            return j_arr_old
        else:
            j_arr_old = new_arr
    print 'Found a full J segment as the V/J junction, ignoring this reconstruction'
    return 'NA'


def find_js_per_len(cur_seq, fasta_dict, id_name_dict, trim):
    f_arr = []
    for seq in fasta_dict:
        if id_name_dict[seq].find('J') != -1:
            j_seq = fasta_dict[seq]
            len_j = len(j_seq)
            for i in range(0, len_j):
                if ((i + trim) <= len_j):
                    trim_j = j_seq[i:i + trim]
                    if cur_seq.find(trim_j) != -1:
                        if seq not in f_arr:
                            f_arr.append(seq)
                            break
    return f_arr


def write_record(tcr_record, tcr_seq, added_c, v_ens, j_ens, v_seq, j_seq, map_dict, bases, c_seq, c_id, c_name, out_f,
                 fasta_dict,
                 rec_name_arr):
    v_seq_trim = ''
    j_seq_trim = ''
    if bases == -10:
        bases = min(len(v_seq), len(j_seq))
    elif bases > len(j_seq):
        j_seq = j_seq + c_seq
        added_c = True
    found = False
    for i in reversed(range(20, bases)):
        junc_start = tcr_seq[:i]
        v_ind = v_seq.find(junc_start)
        if (v_ind != -1):
            found = True
            v_seq_trim = v_seq[:v_ind]
            break
    if found == False:
        v_seq_trim = v_seq[:-bases]
    found = False
    for j in reversed(range(20, bases)):
        junc_end = tcr_seq[-j:]
        j_ind = j_seq.find(junc_end)
        if (j_ind != -1):
            found = True
            j_seq_trim = j_seq[j_ind + j:]
            break
    if found == False:
        j_seq_trim = j_seq[bases:]
    # Add TRBC or TRAC
    c_arr = []
    if (str(tcr_record.id).find('TRB') != -1):
        for ens in map_dict:
            if map_dict[ens].find('TRBC') != -1:
                c_arr.append(ens)
    elif (str(tcr_record.id).find('TRA') != -1):
        for ens in map_dict:
            if map_dict[ens].find('TRAC') != -1:
                c_arr.append(ens)
    else:
        sys.stderr.write(str(datetime.datetime.now()) + " Error! no TRBC or TRAC\n")
        sys.stderr.flush()
    if not added_c:
        for ens in c_arr:
            c_seq = fasta_dict[ens]
            new_seq = v_seq_trim + tcr_seq + j_seq_trim + c_seq
            new_id = map_dict[v_ens] + '.' + map_dict[j_ens] + '.' + map_dict[
                ens] + '.' + v_ens + '.' + j_ens + '.' + ens
            while new_id in rec_name_arr:
                new_id += '_2'
            rec_name_arr.append(new_id)
            record = SeqRecord(Seq(new_seq, IUPAC.ambiguous_dna), id=new_id, description='')
            SeqIO.write(record, out_f, 'fasta')
    else:
        new_seq = v_seq_trim + tcr_seq + j_seq_trim
        new_id = map_dict[v_ens] + '.' + map_dict[j_ens] + '.' + c_name + '.' + v_ens + '.' + j_ens + '.' + c_id
        while new_id in rec_name_arr:
            new_id += '_2'
        rec_name_arr.append(new_id)
        record = SeqRecord(Seq(new_seq, IUPAC.ambiguous_dna), id=new_id, description='')
        SeqIO.write(record, out_f, 'fasta')
    return rec_name_arr


def pick_final_isoforms(full_tcr_file_alpha, full_tcr_file_beta, output):
    pick_final_isoform_chain(full_tcr_file_alpha, output + '.alpha.full.TCRs.bestIso.fa',
                             output + '.alpha.rsem.out.genes.results')
    pick_final_isoform_chain(full_tcr_file_beta, output + '.beta.full.TCRs.bestIso.fa',
                             output + '.beta.rsem.out.genes.results')


def pick_final_isoform_chain(fullTCRfa, newFasta, rsemF):
    if os.path.isfile(fullTCRfa):
        f = open(fullTCRfa, 'rU')
        out_fa = open(newFasta, 'w')
        fasta_dict = dict()
        by_vj_dict = dict()
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in fasta_dict:
                # record.id = record.id + '_2'
                sys.stderr.write(
                    str(datetime.datetime.now()) + 'error! same name for two fasta entries %s\n' % record.id)
                sys.stderr.flush()
            fasta_dict[record.id] = record.seq
            only_vj_rec = str(record.id)
            id_arr = only_vj_rec.strip('\n').split('.')
            vj_str = id_arr[0] + '.' + id_arr[1]
            if vj_str not in by_vj_dict:
                by_vj_dict[vj_str] = []
            by_vj_dict[vj_str].append(record.id)
        for vj_str in by_vj_dict:

            if len(by_vj_dict[vj_str]) == 1:
                c_id = by_vj_dict[vj_str][0]
                c_seq = fasta_dict[c_id]
                new_rec = SeqRecord(c_seq, id=c_id, description='')
                SeqIO.write(new_rec, out_fa, 'fasta')
            else:
                # print vjStr
                # print byVJDict[vjStr]
                best_id = find_best_c(by_vj_dict[vj_str], rsemF)
                # print "best: " + bestId
                b_seq = fasta_dict[best_id]
                new_rec = SeqRecord(b_seq, id=best_id, description='')
                SeqIO.write(new_rec, out_fa, 'fasta')
        out_fa.close()
        f.close()


def find_best_c(vj_arr, rsem_f):
    if (os.path.exists(rsem_f)):
        f = open(rsem_f, 'r')
        f.readline()
        l = f.readline()
        best_seq = 'name'
        max_count = 0.0
        while l != '':
            l_arr = l.strip('\n').split('\t')
            if l_arr[0] in vj_arr:
                curr_count = float(l_arr[4])
                if curr_count > max_count:
                    best_seq = l_arr[0]
                    max_count = curr_count
            l = f.readline()
        f.close()
        if best_seq == 'name':
            return vj_arr[0]
        return best_seq
    else:
        return vj_arr[0]


def find_cdr3(fasta, aa_dict, vdj_fa_dict):
    f = open(fasta, 'rU')
    f_dict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        if record.id in f_dict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! same name for two fasta entries %s\n' % record.id)
            sys.stderr.flush()
        else:
            id_arr = record.id.split('.')
            v_seg = id_arr[0]
            j_seg = id_arr[1]
            if ((v_seg in aa_dict) & (j_seg in aa_dict)):
                curr_dict = find_v_and_j_aa_map(aa_dict[v_seg], aa_dict[j_seg], record.seq)
            else:
                if v_seg in aa_dict:
                    new_vseg = aa_dict[v_seg]
                else:
                    v_id = id_arr[3]
                    curr_seq = vdj_fa_dict[v_id]
                    new_vseg = get_best_v_aa(Seq(curr_seq))
                if j_seg in aa_dict:
                    new_jseg = aa_dict[j_seg]
                else:
                    j_id = id_arr[4]
                    curr_seq = vdj_fa_dict[j_id]
                    new_jseg = get_best_j_aa(Seq(curr_seq))
                curr_dict = find_v_and_j_aa_map(new_vseg, new_jseg, record.seq)
            f_dict[record.id] = curr_dict
    f.close()
    return f_dict


def get_nt_seq(full_seq):
    mod = len(full_seq) % 3
    if mod != 0:
        f_seq = full_seq[:-mod].translate()
    else:
        f_seq = full_seq.translate()
    return f_seq


def find_v_and_j_aa_map(v_seg, j_seg, full_seq):
    f_dict = dict()
    first_seq = get_nt_seq(full_seq)
    second_seq = get_nt_seq(full_seq[1:])
    third_seq = get_nt_seq(full_seq[2:])
    nt_arr = [full_seq, full_seq[1:], full_seq[2:]]
    aa_seqs_arr = [first_seq, second_seq, third_seq]
    cdr_arr = []
    pos_arr = []
    f_ph_arr = []
    for aa_seq in aa_seqs_arr:
        (cdr, pos, cur_ph) = get_cdr3(aa_seq, v_seg, j_seg)
        cdr_arr.append(cdr)
        pos_arr.append(pos)
        f_ph_arr.append(cur_ph)
    max_len = 0
    best_cdr = ''
    best_seq = ''
    has_stop = False
    best_pos = -1
    best_cdr_nt = ''
    found_good = False
    v_pos = -1
    j_pos = -1
    f_ph = False
    for i in range(0, 3):
        if pos_arr[i] != -1:
            if ((cdr_arr[i] != 'Only J') & (cdr_arr[i] != 'Only V')):
                if len(cdr_arr[i]) > max_len:
                    if cdr_arr[i].find('*') == -1:
                        found_good = True
                        best_cdr = cdr_arr[i]
                        best_pos = pos_arr[i]
                        max_len = len(cdr_arr[i])
                        best_seq = nt_arr[i]
                        f_ph = f_ph_arr[i]
                    else:
                        if max_len == 0:
                            found_good = True
                            best_pos = pos_arr[i]
                            best_cdr = cdr_arr[i]
                            max_len = len(cdr_arr[i])
                            best_seq = nt_arr[i]
                            has_stop = True
                            f_ph = f_ph_arr[i]
                else:
                    if has_stop == True:
                        if cdr_arr[i].find('*') == -1:
                            found_good = True
                            best_pos = pos_arr[i]
                            has_stop = False
                            best_cdr = cdr_arr[i]
                            max_len = len(cdr_arr[i])
                            best_seq = nt_arr[i]
                            f_ph = f_ph_arr[i]
            else:
                if not found_good:
                    f_ph = f_ph_arr[i]
                    if (cdr_arr[i] == 'Only J'):
                        j_pos = pos_arr[i] - i
                    elif (cdr_arr[i] == 'Only V'):
                        v_pos = pos_arr[i] - i
    if ((v_pos != -1) & (j_pos != -1) & (not found_good)):
        best_cdr_nt = full_seq[3 * v_pos:3 * j_pos]
        best_cdr = 'NA'
    elif best_pos != -1:
        best_cdr_nt = best_seq[3 * best_pos: 3 * best_pos + 3 * len(best_cdr)]
    if best_cdr.find('*') != -1:
        stat = 'Unproductive - stop codon'
    elif f_ph:
        stat = 'Productive (no 118 PHE found)'
    else:
        stat = 'Productive'
    if max_len == 0:
        if (('Only J' in cdr_arr) & ('Only V' in cdr_arr)):
            stat = 'Unproductive - Frame shift'
        else:
            if (('Only J' not in cdr_arr) & ('Only V' in cdr_arr)):
                stat = 'Unproductive - found only V segment'
            elif (('Only J' in cdr_arr) & ('Only V' not in cdr_arr)):
                stat = 'Unproductive - found only J segment'
            elif (('Only J' not in cdr_arr) & ('Only V' not in cdr_arr)):
                stat = 'Unproductive - didn\'t find V and J segment'
            else:
                stat = 'Unproductive'
            best_cdr = 'NA'
            best_cdr_nt = 'NA'
    f_dict['stat'] = stat
    f_dict['CDR3 AA'] = best_cdr
    f_dict['CDR3 NT'] = best_cdr_nt
    f_dict['Full Seq'] = full_seq
    return f_dict


def get_cdr3(aa_seq, v_seq, j_seq):
    min_dist = 14
    pos = -1
    for i in range(0, len(aa_seq) - len(v_seq) + 1):
        sub_aa = aa_seq[i:i + len(v_seq)]
        if len(sub_aa) != len(v_seq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong sub length\n')
            sys.stderr.flush()
        dist = 0
        for k in range(0, len(v_seq)):
            if v_seq[k] != sub_aa[k]:
                dist += 1
        if ((dist < min_dist) & (sub_aa.endswith('C'))):
            min_dist = dist
            pos = i + len(v_seq)
    j_pos = -1
    min_dist_j = 4
    cur_ph = False
    for j in range(pos + 1, len(aa_seq) - len(j_seq) + 1):
        sub_aa = aa_seq[j: j + len(j_seq)]
        if len(sub_aa) != len(j_seq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0, len(j_seq)):
            if j_seq[m] != sub_aa[m]:
                dist += 1
        if (dist <= min_dist_j):
            if is_legal(sub_aa):
                j_pos = j
                min_dist_j = dist
                cur_ph = False
            else:
                if dist < min_dist_j:
                    cur_ph = True
                    j_pos = j
                    min_dist_j = dist
    if pos == -1:
        if j_pos != -1:
            return ('Only J', j_pos, cur_ph)
        else:
            return ('No V/J found', -1, cur_ph)
    else:
        if j_pos == -1:
            return ('Only V', pos, cur_ph)
    return (aa_seq[pos:j_pos], pos, cur_ph)


# Checks that the conserved amino acids remain
def is_legal(sub_aa):
    if (len(sub_aa) < 4):
        return False
    if (sub_aa[0] == 'F'):
        if ((sub_aa[1] == 'G') | (sub_aa[3] == 'G')):
            return True
    if ((sub_aa[1] == 'G') & (sub_aa[3] == 'G')):
        return True
    if ((sub_aa[0] == 'G') & (sub_aa[2] == 'G')):
        return True
    if sub_aa.startswith('GR'):
        return True
    if sub_aa.startswith('SR'):
        return True
    return False


def get_best_j_aa(curr_seq):
    first_seq = get_nt_seq(curr_seq)
    second_seq = get_nt_seq(curr_seq[1:])
    third_seq = get_nt_seq(curr_seq[2:])
    pos = 10
    seq = ''
    found = False
    for s in [first_seq, second_seq, third_seq]:
        temp_seq = s[:8]
        ind_f = temp_seq.find('F')
        if ind_f != -1:
            if ((ind_f + 3) <= len(s)):
                if s[ind_f + 3] == 'G':
                    found = True
                    if ind_f < pos:
                        pos = ind_f
                        seq = s[ind_f:]
        ind_g = temp_seq.find('G')
        if (ind_g != -1):
            if found == False:
                if ((ind_g + 2) <= len(s)):
                    if s[ind_g + 2] == 'G':
                        if ind_g < pos:
                            found = True
                            seq = s[ind_g:]
                            pos = ind_g
    if ((found == False) & (ind_f != -1)):
        seq = s[ind_f:]
    if seq != '':
        return seq
    else:
        return first_seq


def get_best_v_aa(curr_seq):
    first_seq = get_nt_seq(curr_seq)
    second_seq = get_nt_seq(curr_seq[1:])
    third_seq = get_nt_seq(curr_seq[2:])
    pos = 10
    seq = ''
    for s in [first_seq, second_seq, third_seq]:
        # print "S: " + s
        temp_seq = s[-8:]
        # print "tempSeq: " + tempSeq
        ind = temp_seq.find('C')
        stop_ind = temp_seq.find('*')
        # print "Ind: " + str(ind)
        if ((ind != -1) & (stop_ind == -1)):
            # print "inside the ind"
            if ind < pos:
                good_rf = is_good_rf(s)
                if good_rf:
                    pos = ind
                    seq = s[:-8 + ind + 1]
    if seq != '':
        return seq
    else:
        return first_seq


def is_good_rf(s):
    m_ind = s.find('M')
    if m_ind == -1:
        return False
    stop_ind = s.find('*')
    if stop_ind == -1:
        return True
    stop_ind_next = s[stop_ind + 1:].find('*')
    while stop_ind_next != -1:
        stop_ind = stop_ind_next + stop_ind + 1
        stop_ind_next = s[stop_ind + 1:].find('*')
        m_ind = s[stop_ind + 1:].find('M')
        m_ind = m_ind + stop_ind + 1
    if m_ind != -1:
        return True
    else:
        return False
