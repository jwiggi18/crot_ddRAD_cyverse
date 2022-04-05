"""
rsw.py

from Gamble et al. 2015 Molecular Biology and Evolution

This script executes a particular workflow on a set of files.  There are
approximately four steps.

Step 1:



default behavior:
    input files will be found in the current directory
    the input files will be named
        male_reads_1.fastq
        male_reads_2.fastq
        female_reads_1.fastq
        female_reads_2.fastq

  print:
    # of sex-specific tags identified by RADtools
    # of confirmed sex-specific tags

"""
import argparse
import os.path


def rsw(working_dir_path,
        allele_counts_file_name,
        male_reads_1_file_name,
        male_reads_2_file_name,
        female_reads_1_file_name,
        female_reads_2_file_name,
        rad_tools_output_file_name,
        max_row_cluster_tag_count,
        n_minus,
        seq_pattern_sex):

    print('working directory: {}'.format(working_dir_path))

    #rad_tools_output_file_path = os.path.join(working_dir_path, rad_tools_output_file_name)
    rad_tools_output_file_path = os.path.join(working_dir_path, rad_tools_output_file_name)
    allele_counts_file_path = os.path.join(working_dir_path, allele_counts_file_name)

    #  Step 0
    #  Accumulate allele counts per ClusterID from the RADtools output file.
    #
    write_allele_counts_file(rad_tools_output_file_path, allele_counts_file_path)

    #  Step 1
    #  The first step works with a RAD-Seq ?? output file that looks something like this:
    #    ClusterID  ClusterTags  SegPattern  Tag    tg1436_male  tg1437_female  tg1541_female  tg1578_male
    #    1          1            1--1        TGCAG  10           NA             NA              7
    #    2          1            1--1        TGCAG   9           NA             NA              9
    #    3          1            1111        TGCAG  10            6             12             18
    #    4          4            1-11        TGCAG  12           NA             10              5
    #    4          4            111-        TGCAG  20            5              6             NA
    #    4          4            1111        TGCAG  14            6              7             14
    #    4          4            1111        TGCAG  14            7              6             13
    #
    #  The goal of Step 1 is to extract rows with ClusterTags = less that some integer (set by the adjusting 
    #  the -m command) and single-sex SeqPattern (that can be adjusted using the -n command).  The argument
    #  seq_pattern_sex might look like 'mffm'.  This indicates the first and last individuals are male and
    #  the second and third individuals are female.  Rows with SeqPattern 1---, ---1, or 1--1 have tags found
    #  only in the male individuals, while rows with SeqPattern -1--, --1-, or -11- have tags found only in
    #  the female individuals.  There is a minimum number of individuals needed to consider a row as having
    #  a single-sex tag.  By default this minimum number is one fewer than the total number of individuals
    #  of a single sex.  Confusing? <perhaps, I will work on this>  For example, if we have a total of three 
    #  males and four females then by default a tag is considered a putative male-specific tag if it is 
    #  found in at least 3-1=2 male individuals. Similarly, a tag is considered a putative female-specific 
    #  tag if it is found in at least 4-1=3 individuals.

    print('maximum cluster tags per row: {}'.format(max_row_cluster_tag_count))

    # individuals_count_table is a dictionary with two items like this:
    #   {
    #     'f': 4,
    #     'm': 3
    #   }
    # in the example above there are 4 female individuals and 3 male individuals
    # this dictionary records the total number of male and female individuals
    # the values in this dictionary do not change
    individuals_count_table = {'f': 0, 'm': 0}
    for sex in seq_pattern_sex:
        individuals_count_table[sex] += 1

    # print the number of individuals of each sex
    print('thresholds for sex-specific tags:')
    for (sex, count) in individuals_count_table.items():
        sex_specific_threshold = count - n_minus
        if sex_specific_threshold < 1:
            raise Exception('sex-specific threshold for {} individuals is less than one'.format(sex))
        print('  {} of {} {} individual(s)'.format(count - n_minus, count, sex))

    row_data_by_sex = accumulate_row_data_by_sex(
        rad_tools_output_file_path,
        seq_pattern_sex,
        max_row_cluster_tag_count,
        n_minus,
        individuals_count_table
    )
    write_tags_file('female_specific_tags.txt', row_data_by_sex['f'])
    write_tags_file('male_specific_tags.txt', row_data_by_sex['m'])
    print('found {} male-specific tags'.format(len(row_data_by_sex['m'])))
    print('found {} female-specific tags'.format(len(row_data_by_sex['f'])))

    #
    #  Step 2
    #

    # need to read FASTQ files
    # for each male-specific Tag check the female_reads_1.fastq file for a match
    # keep unmatched tags in confirmed_male_specific_tags
    female_reads_1_file_path = os.path.join(working_dir_path, female_reads_1_file_name)

    male_tags_in_female_reads_counts = count_matching_tags(row_data_by_sex['m'], female_reads_1_file_path)
    write_tag_counts_file('confirmed_male_specific_tags.txt', male_tags_in_female_reads_counts)
    male_reads_1_file_path = os.path.join(working_dir_path, male_reads_1_file_name)

    female_tags_in_male_reads_counts = count_matching_tags(row_data_by_sex['f'], male_reads_1_file_path)
    write_tag_counts_file('confirmed_female_specific_tags.txt', female_tags_in_male_reads_counts)

    #
    #  Step 3
    #  For each sex-specific tag find matching reads in the corresponding forward
    #  reads file.
    #

    male_specific_forward_read_table = write_sex_specific_forward_read_files(
        male_tags_in_female_reads_counts,
        male_reads_1_file_path,
        'confirmed_male_specific_reads_{}_1.fastq'
    )
    female_specific_forward_read_table = write_sex_specific_forward_read_files(
        female_tags_in_male_reads_counts,
        female_reads_1_file_path,
        'confirmed_female_specific_reads_{}_1.fastq'
    )

    #
    # step 4
    #
    male_reads_2_file_path = os.path.join(working_dir_path, male_reads_2_file_name)
    write_sex_specific_paired_read_files(
        male_specific_forward_read_table,
        male_reads_2_file_path,
        'confirmed_male_specific_reads_{}_2.fastq'
    )

    female_reads_2_file_path = os.path.join(working_dir_path, female_reads_2_file_name)
    write_sex_specific_paired_read_files(
        female_specific_forward_read_table,
        female_reads_2_file_path,
        'confirmed_female_specific_reads_{}_2.fastq'
    )


#
# This function handles step 0
# Accumulate allele counts for each ClusterID and write them
# to a new file.
#
def write_allele_counts_file(rad_tools_output_file_path, allele_counts_file_path):
    cluster_id_totals = {}
    with open(rad_tools_output_file_path) as rad_tools_output:
        # read the header
        header_row = rad_tools_output.readline()
        header_row_values = header_row.strip().split()
        tag_index = header_row_values.index('Tag')
        allele_count_start_index = tag_index + 1
        for line in rad_tools_output.readlines():
            row_values = line.strip().split()
            row_cluster_id = int(row_values[0])
            allele_counts = []
            for count_string in row_values[allele_count_start_index:]:
                if count_string == 'NA':
                    allele_counts.append(0)
                else:
                    allele_counts.append(int(count_string))

            if row_cluster_id not in cluster_id_totals:
                # concatenate the first four columns of the current row with
                # the allele_counts list
                # insert the first of the row values for this cluster id
                cluster_id_totals[row_cluster_id] = row_values[:allele_count_start_index] + allele_counts
            else:
                # add allele_counts to the existing counts
                cluster_id_row_values = cluster_id_totals[row_cluster_id]
                for i, allele_count in enumerate(allele_counts):
                    cluster_id_row_values[allele_count_start_index+i] += allele_count

    with open(allele_counts_file_path, 'w') as allele_counts_file:
        allele_counts_file.write(header_row)
        for cluster_id in sorted(cluster_id_totals.keys()):
            allele_counts_file.write(' '.join([str(x) for x in cluster_id_totals[cluster_id]]))
            allele_counts_file.write('\n')


#
# This function handles step 1
#
def accumulate_row_data_by_sex(
        rad_tools_output_file_path,
        seq_pattern_sex,
        max_row_cluster_tag_count,
        n_minus,
        individuals_count_table):
    row_data_by_sex = {'m': {}, 'f': {}}
    # open the RAD tools output file
    with open(rad_tools_output_file_path) as rad_tools_output:
        # read the header
        # ClusterID ClusterTags SeqPattern Tag <name 1> ... <name N>
        header_row = rad_tools_output.readline()
        header_row_values = header_row.strip().split()
        # read all remaining lines
        # keep lines with ClusterTags = x (set using the -m command) that are also single-sex
        for line in rad_tools_output.readlines():
            # line is a string that looks like this
            #   '1 1 --11 ATGC NA NA 6 9\n'
            # row_values is a list of strings such as
            #   [ '1', '1', '--11', 'ATGC', 'NA', 'NA', '6', '9' ]
            row_values = line.strip().split()
            row_cluster_tag_count = int(row_values[1])
            row_individuals_with_tag = row_values[2]
            # zip('--11', 'mffm') returns a combination of the two input lists
            # that looks like this:
            #   [('-', 'm'),('-', 'f'),('1', 'f'), ('1', 'm')]
            # construct a dictionary row_sex_table that looks like this:
            #   {'f': 1, 'm': 1}
            # only enter a key in this dictionary if there is at least one individual of
            # the corresponding sex with the given tag
            row_sex_table = {}
            for (present, sex) in zip(row_individuals_with_tag, seq_pattern_sex):
                if present == '1':
                    if sex not in row_sex_table:
                        row_sex_table[sex] = 0
                    row_sex_table[sex] += 1
            # the current row is a single-sex row if row_sex_table has only one element
            #print('row cluster tag count: {}'.format(row_cluster_tag_count))
            #print('row sex table: {}'.format(row_sex_table))
            if (max_row_cluster_tag_count >= row_cluster_tag_count) and len(row_sex_table) == 1:
                #print('row accepted')
                # there is only one element in row_sex_table and we can get it with popitem()
                # sex is either 'm' or 'f'
                (sex, individuals_with_tag_count) = row_sex_table.popitem()
                individuals_count = individuals_count_table[sex]
                if (individuals_count - n_minus) < 1:
                    # this is a degenerate situation where we have something like
                    #   individuals_count = 3
                    #   n_minus = 3
                    # in cases like this tags will never be sex-specific
                    pass
                elif individuals_with_tag_count < (individuals_count - n_minus):
                    # in this case we have something like
                    #   individuals_with_tag_count = 2
                    #   individuals_count = 4
                    #   n_minus = 1
                    # so we will not count this tag as being sex-specific since
                    # this tag does not meet the threshold of individuals for a
                    # putative single-sex tag
                    pass
                else:
                    # in this case the individuals_with_tag_count meets the
                    # threshold for a putative sex-specific tag so we will
                    # store the row data in row_data_by_sex
                    # dict(zip(header_row_values, row_values)) looks like this:
                    #   {
                    #      'ClusterID'     : '1',
                    #      'ClusterTags'   : '1',
                    #      'SeqPattern'    : '--11',
                    #      'Tag'           : 'ATGC',
                    #      'tg1436_male'   : 'NA',
                    #      'tg1437_female' : 'NA',
                    #      'tg1541_female' : '6',
                    #      'tg1578_male'   : '9'
                    #  }
                    #print('header row values: {}'.format(header_row_values))
                    #print('row values: {}'.format(row_values))
                    row_data = dict(zip(header_row_values, row_values))
                    row_data_by_sex[sex][row_data['Tag']] = row_data
                    # row_data_by_sex looks like
                    #   {
                    #      'm': {
                    #              'ACTG': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'ACTT': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'ACTA': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'ACTC': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #           }
                    #      'f': {
                    #              'GCTG': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'GCTT': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'GCTA': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #              'GCTC': {'ClusterID': ..., 'ClusterTags': ..., ...}
                    #           }
                    #   }
            else:
                pass
    return row_data_by_sex


def write_tags_file(file_path, row_data):
    with open(file_path, 'w') as tags_file:
        for (tag, data) in row_data.items():
            tags_file.write(tag)
            tags_file.write('\n')


def write_tag_counts_file(file_path, tag_counts):
    with open(file_path, 'w') as tag_counts_file:
        for (tag, count) in tag_counts.items():
            tag_counts_file.write('{:<5} {}\n'.format(count, tag))


# the FASTQ files can be very large
# we need to deal with them one record at a time
# so we'll use a fancy generator 'read_fastq'
def count_matching_tags(sex_specific_tags, opposite_sex_fastq_file_path):
    tag_counter = {tag: 0 for tag in sex_specific_tags.keys()}
    for (fqid, seq, plus, quality) in read_fastq(opposite_sex_fastq_file_path):
        for tag in tag_counter.keys():
            if tag in seq:
                #print('found tag {} in sequence {} from file {}'.format(tag, seq, opposite_sex_fastq_file_path))
                tag_counter[tag] += 1
            else:
                pass
    #print('{}'.format(tag_counter))
    return tag_counter


# this function handles step 3
# also return a dictionary useful for step 4 that looks like this:
#   {
#     <forward reads output file path>: (n, {dictionary of forward reads keyed by FASTQ id}),
#     <forward reads output file path>: (n, {dictionary of forward reads keyed by FASTQ id}),
#     <forward reads output file path>: (n, {dictionary of forward reads keyed by FASTQ id})
#   }
#
#  forward read:
#  @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633 1:N:0:ATTCCT
#  AGAGTTGCAGGCCACACACACAAATTCAAAGAGTCCACAACCGCCAGGGCTTGCAGTTGAGTGCTGAGCGCCGAGTGTAGCCCGGAGCGTATCTTACCAG
#  +
#  @@?DDBDDFHHHFIADHH1FEG;EHHIGGGI@F?FHHCHB?DFGIEE;@166=EHHBCCBD7;C@C@;;>B:;;@8<>4:@:>@<9<>B58<A@CDCCCC
#
#  paired read:
#  @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633 2:N:0:ATTCCT
#  GACTGCCAGGGCTTGCAGCTGAGCGCCAAGAGCCTGGTAAGATACGCTCCGGGCTACACTCGGCGCTCAGCACTCAACTGCAAGCCCTGGCGGTTGTGGA
#  +
#  B?@=DDDDFD::CFGGBGGBAHH9?DHG>=GCG)9DE*BF@DHH)-<1@F;BHB>B;@A3>=?B@?@BB?CCBCDACC>34@CDA<898AB95)5<9?##
def write_sex_specific_forward_read_files(tag_counts, forward_reads_file_path, output_file_name_template):
    # first make a dictionary with sex-specific tags as keys and empty lists as values
    sex_specific_forward_reads = {
        tag: []
        for (tag, count)
        in tag_counts.items()
        if count == 0
    }
    # keep a dictionary of tags found in the forward read file like this
    #   {
    #      tag : [list of 4-tuples],
    #      tag : [list of 4-tuples], ...
    #   }
    for (fqid, seq, plus, quality) in read_fastq(forward_reads_file_path):
        for tag in sex_specific_forward_reads.keys():
            if tag in seq:
                #print('found tag {} in sequence {}'.format(tag, seq))
                sex_specific_forward_reads[tag].append((fqid, seq, plus, quality))
            else:
                pass

    # do we need to put back the sequence prefix and suffix? Yes.
    # do we need to write all four lines? Yes.
    # n counts the number of files written
    n = 0
    # forward_read_table looks like
    #   {
    #     <forward read file path>: (n, {dictionary of forward reads keyed by FASTQ id}),
    #     ... more like above
    #   }
    #
    # the key is a FASTQ id and the value is the file number
    # of the corresponding forward read file
    forward_read_table = {}
    for (tag, reads) in sex_specific_forward_reads.items():
        # read is a 4-ple of FASTQ lines
        n += 1
        output_file_name = output_file_name_template.format(n)
        with open(output_file_name, 'w') as confirmed_read_file:
            for read in reads:
                # output_file_name_template looks something like this
                #   'confirmed_male_specific_reads_{}_1.txt'
                confirmed_read_file.write('\n'.join(read))
                confirmed_read_file.write('\n')
                # the first element of read looks like this:
                #   @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633 1:N:0:ATTCCT
                # the id value will look like this:
                #   @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633
                fqid = read[0].split()[0]
                if output_file_name not in forward_read_table:
                    forward_read_table[output_file_name] = (n, {})
                forward_read_file_table = forward_read_table[output_file_name][1]
                # are the fqid unique?
                forward_read_file_table[fqid] = read
    return forward_read_table


# this function handles step 4
# the first argument is a dictionary like this:
#   {
#     <forward read file path>: (n, {dictionary of FASTQ 4-tuples keyed by FASTQ id}),
#     <forward read file path>: (n, {dictionary of FASTQ 4-tuples keyed by FASTQ id}),
#     <forward read file path>: (n, {dictionary of FASTQ 4-tuples keyed by FASTQ id}),
#   }
#
#  the second argument is the file path of the paired reads FASTQ file
#  the third argument is the template for the paired read output files
#  and looks something like this:
#    'female_specific_paired_read_{}_2.txt'
#
def write_sex_specific_paired_read_files(sex_specific_forward_read_table, paired_reads_file_path, output_file_name_template):
    # read the paired read FASTQ file
    # look in the sex_specific_forward_read_table for each paired id
    # if the paired id is found then write the corresponding paired read
    # to a file using the associated value in sex_specific_forward_read_table
    # as input to output_file_name_template.format()

    # accumulate all paired reads and then write them to files

    # sex_specific_paired_read_table accumulates lists of paired-read 4-ples that
    # will be written to the same paired reads output file; it is keyed with the
    # path of the output file with corresponding foward reads
    # sex_specific_paired_read_table looks like this
    #   {
    #      <forward reads file name>: (n, [list of paired read FASTQ 4-tuples])
    #      <forward reads file name>: (n, [list of paired read FASTQ 4-tuples]), ...
    #   }
    sex_specific_paired_read_table = {}
    for (fqidline, seq, plus, quality) in read_fastq(paired_reads_file_path):
        # the first element of read looks like this:
        #   @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633 1:N:0:ATTCCT
        # the fqid value will look like this:
        #   @HWI-ST1073:280:D18HTACXX:3:2105:19250:113633
        fqid = fqidline.split()[0]
        # look for the paired read fqid in the forward read dictionaries
        for (forward_reads_file_path, (n, forward_reads_table)) in sex_specific_forward_read_table.items():
            ##(n, forward_reads) = sex_specific_forward_read_table[fqid]
            if fqid in forward_reads_table:
                #print('found id {} from forward read file {} in paired read file {}'.format(
                #    fqid, forward_reads_file_path, paired_reads_file_path)
                #)
                if forward_reads_file_path not in sex_specific_paired_read_table:
                    sex_specific_paired_read_table[forward_reads_file_path] = (n, [])
                paired_reads = sex_specific_paired_read_table[forward_reads_file_path][1]
                paired_reads.append((fqidline, seq, plus, quality))
                break
            else:
                pass

    for (forward_reads_file_path, (n, paired_reads)) in sex_specific_paired_read_table.items():
        output_file_path = output_file_name_template.format(n)
        # only open this file once
        with open(output_file_path, 'w') as confirmed_paired_read_file:
            #print('writing {} paired reads to file {}'.format(len(paired_reads), output_file_path))
            for paired_read in paired_reads:
                confirmed_paired_read_file.write('\n'.join(paired_read))
                confirmed_paired_read_file.write('\n')


# the tags in the FASTQ files have prefix and suffix stuff around the tag
# as we know it from the RADseq output file
# this function returns groups of 4 lines from a FASTQ file as a 4-ple
def read_fastq(fastq_file_path):
    with open(fastq_file_path) as fastq_file:
        read_count = 0
        fastq_id = fastq_file.readline().strip()
        while len(fastq_id) > 0:
            read_count += 1
            seq = fastq_file.readline().strip()
            plus = fastq_file.readline().strip()
            quality = fastq_file.readline().strip()
            yield (fastq_id, seq, plus, quality)
            fastq_id = fastq_file.readline().strip()
    print('found {} reads in {}'.format(read_count, fastq_file_path))


def read_cmd_line():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        #'-wd',
        'working_dir_path',
        #default='.',
        help='working directory'
    )
    parser.add_argument(
        #'-mf',
        'seq_pattern_sex',
        #default='mmmmffff',
        help='indicate sex of each individual as a string, eg fmmf'
    )
    parser.add_argument(
        '-m',
        '--max_row_cluster_tag_count',
        type=int,
        default=1,
        help='maximum number of allowed cluster tags per row'
    )
    parser.add_argument(
        '-n',
        '--n_minus',
        type=int,
        default=1,
        help='positive integer to be subtracted from the number of individuals to give the minimum number required for a tag to be considered (putatively) sex-specific'
    )
    parser.add_argument(
        #'-rt',
        'rad_tools_output_file_name',
        #default='RADtools_output.txt',
        help='RAD tools output file'
    )
    parser.add_argument(
        'allele_counts_file_name',
        help='allele counts output file name'
    )
    parser.add_argument(
        #'-m1',
        'male_reads_1_file_name',
        #default='male_reads_1.fastq',
        help='FASTQ file with male forward reads'
    )
    parser.add_argument(
        #'-m2',
        'male_reads_2_file_name',
        #default='male_reads_2.fastq',
        help='FASTQ file with male paired-end reads'
    )
    parser.add_argument(
        #'-f1',
        'female_reads_1_file_name',
        #default='female_reads_1.fastq',
        help='FASTQ file with female forward reads'
    )
    parser.add_argument(
        #'-f2',
        'female_reads_2_file_name',
        #default='female_reads_2.fastq',
        help='FASTQ file with female paired-end reads'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = read_cmd_line()
    rsw(**vars(args))