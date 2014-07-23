__author__ = 'pcmarks'

# This version of the GEO soft data loader uses the SSDB server for the backend
# NB: CHANGE THE VARIABLE outer_directory TO POINT TO THE DIRECTORY THAT CONTAINS THE SOFT
# FILES.

import csv
from json import JSONEncoder
import os
from datetime import datetime
import pyssdb


outer_directory = 'data/'
data_directory = os.path.join(outer_directory, 'raw')

soft_file_attributes = [
    {
        'file_name': 'GSE35781_family.soft',
        'study': 'GSE35781',
        'profile': 'Expression-miRNA',
        'platform': 'Affymetrix_miRNA_Array',
        'gene_symbol_field': None},
    {
        'file_name': 'GSE26125_family.soft',
        'study': 'GSE26125',
        'profile': 'Expression-Genes',
        'platform': 'Human_Whole_Genome_Bioarray',
        'gene_symbol_field': 'Gene_Symbol'},
]

# A hard-coded dictionary of sample attributes:
#   the control attribute indicates whether the sample is control or diseased
#   the co_sample attribute points to the corresponding sample in an allied study.
sample_attributes_dict = {
    'GSM869376': {'control': True, 'co_sample': ''},
    'GSM869377': {'control': True, 'co_sample': ''},
    'GSM869378': {'control': True, 'co_sample': ''},
    'GSM869379': {'control': True, 'co_sample': 'GSM641583'},
    'GSM869380': {'control': True, 'co_sample': 'GSM641584'},
    'GSM869381': {'control': True, 'co_sample': 'GSM641585'},
    'GSM869382': {'control': True, 'co_sample': 'GSM641586'},
    'GSM869383': {'control': True, 'co_sample': 'GSM641587'},
    'GSM869400': {'control': True, 'co_sample': ''},
    'GSM869401': {'control': True, 'co_sample': ''},
    'GSM869402': {'control': True, 'co_sample': ''},

    'GSM641583': {'control': True, 'co_sample': 'GSM869379'},
    'GSM641584': {'control': True, 'co_sample': 'GSM869380'},
    'GSM641585': {'control': True, 'co_sample': 'GSM869381'},
    'GSM641586': {'control': True, 'co_sample': 'GSM869382'},
    'GSM641587': {'control': True, 'co_sample': 'GSM869383'},

    'GSM869399': {'control': False, 'co_sample': ''},
    'GSM869398': {'control': False, 'co_sample': ''},
    'GSM869397': {'control': False, 'co_sample': ''},
    'GSM869396': {'control': False, 'co_sample': ''},
    'GSM869387': {'control': False, 'co_sample': 'GSM641573'},
    'GSM869389': {'control': False, 'co_sample': 'GSM641568'},
    'GSM869384': {'control': False, 'co_sample': 'GSM641576'},
    'GSM869388': {'control': False, 'co_sample': 'GSM641570'},
    'GSM869386': {'control': False, 'co_sample': 'GSM641569'},
    'GSM869390': {'control': False, 'co_sample': 'GSM641575'},
    'GSM869391': {'control': False, 'co_sample': 'GSM641574'},
    'GSM869385': {'control': False, 'co_sample': ''},
    'GSM869394': {'control': False, 'co_sample': ''},
    'GSM869395': {'control': False, 'co_sample': ''},
    'GSM869392': {'control': False, 'co_sample': ''},
    'GSM869393': {'control': False, 'co_sample': ''},

    'GSM641573': {'control': False, 'co_sample': 'GSM869387'},
    'GSM641568': {'control': False, 'co_sample': 'GSM869389'},
    'GSM641576': {'control': False, 'co_sample': 'GSM869384'},
    'GSM641570': {'control': False, 'co_sample': 'GSM869388'},
    'GSM641569': {'control': False, 'co_sample': 'GSM869386'},
    'GSM641575': {'control': False, 'co_sample': 'GSM869390'},
    'GSM641574': {'control': False, 'co_sample': 'GSM869391'},

    'GSM641567': {'control': False, 'co_sample': ''},
    'GSM641571': {'control': False, 'co_sample': ''},
    'GSM641572': {'control': False, 'co_sample': ''},
    'GSM641577': {'control': False, 'co_sample': ''},
    'GSM641578': {'control': False, 'co_sample': ''},
    'GSM641579': {'control': False, 'co_sample': ''},
    'GSM641580': {'control': False, 'co_sample': ''},
    'GSM641581': {'control': False, 'co_sample': ''},
    'GSM641582': {'control': False, 'co_sample': ''},

}

## Database codes
source_code = '00'  # GEO, TCGA
cancer_code = '01'
profile_code = '02'
platform_code = '03'
sample_id_code = '04'
gene_code = '05'
vital_status_code = '06'
days_to_death_code = '07'
days_to_birth_code = '08'
days_to_last_followup_code = '09'
study_code = '10'  # For GEO only: GSEnnnnn number
sample_attributes_code = '11'  # dictionary of 'att': 'val' 's

cancers_code = '51'
profiles_code = '52'
platforms_code = '53'
sample_ids_code = '54'
genes_code = '55'
studies_code = '56'

db = pyssdb(host='localhost', port=8888)


def Put(key, value):
    """
     Just a convenience function
    :param key:
    :param value:
    """
    db.set(key, value)


"""
Begin the process:

"""
print "Start", datetime.now()
studies = []
for soft_file_attribute in soft_file_attributes:
    genes_list = []
    study = soft_file_attribute['study']
    studies.append(study)
    study_key = '|'.join([study_code, study])
    Put(study_key, '')

    profile = soft_file_attribute['profile']
    profiles = [profile]
    profile_key = '|'.join([study_code, study, profile_code, profile])
    Put(profile_key, '')

    platform = soft_file_attribute['platform']
    platforms = [platform]
    platform_key = '|'.join([profile_key, platform_code, platform])
    Put(platform_key, '')

    # Open the soft file for reading
    soft_file_name = soft_file_attribute['file_name']
    soft_file = open(os.path.join(data_directory, soft_file_name))
    line = soft_file.next()
    tokens = line.split(' = ')

    while '!Series_sample_id' != tokens[0]:
        line = soft_file.next()
        tokens = line.split(' = ')

    # Gather all of the Sample IDs only if they appear in the sample_attributes dictionary
    sample_ids = []
    while tokens[0] == '!Series_sample_id':
        sample_id = tokens[1].strip()
        if sample_id in sample_attributes_dict:
            sample_ids.append(sample_id)
        line = soft_file.next()
        tokens = line.split(' = ')
    print "number of samples", len(sample_ids)

    # Read until the first platform_table_begin line appears
    # Use this table to build the list of gene symbols for this platform
    while line.strip() != '!platform_table_begin':
        line = soft_file.next()
        tokens = line.split(' = ')
        if tokens[0] == '!Platform_data_row_count':
            number_of_symbols = int(tokens[1])

    print 'number of symbols', number_of_symbols

    # From this point open the soft file as a comma-separated value file
    # Read all of the gene symbols and store in a dictionary that will be used in
    # in the next phase.
    # Stop reading when the token '!platform_table_end' is read
    reader = csv.DictReader(soft_file, delimiter='\t')
    gene_symbol_dict = {}
    gene_symbol_field = soft_file_attribute['gene_symbol_field']
    if gene_symbol_field:
        for line in reader:
            if line['ID'] == '!platform_table_end':
                break
            gene_symbol = line[gene_symbol_field]
            if len(gene_symbol) > 0:
                gene_symbol_dict[line['ID']] = gene_symbol
    else:
        for line in reader:
            if line['ID'] == '!platform_table_end':
                break


    # Now begin reading the expression values for every sample, but only keep
    # those samples that appear in the sample_attributes dictionary
    line = soft_file.next()
    sample_count = 0
    while line:
        tokens = line.split(" = ")
        sample_id = tokens[1].strip('\n')
        if sample_id not in sample_attributes_dict:
            # Skip all the way to the end of this sample's expressions
            line = soft_file.next()
            while line:
                if line.startswith('!sample_table_end'):
                    try:
                        line = soft_file.next()  # Read past it
                    except StopIteration:
                        line = None
                        break
                    break
                line = soft_file.next()
        else:
            sample_id_key = '|'.join([platform_key, sample_id_code, sample_id])
            Put(sample_id_key, '')
            sample_attributes_key = '|'.join([platform_key, sample_id_code, sample_id, sample_attributes_code])
            sample_attributes = sample_attributes_dict[sample_id]
            Put(sample_attributes_key, JSONEncoder().encode(sample_attributes))
            print "values for", sample_id, "attributes:", sample_attributes
            # Read forward until the sample expression table begins
            line = soft_file.next()
            while line:
                if line.startswith('!sample_table_begin'):
                    break
                line = soft_file.next()

            # Read the embedded csv table containing the expression values
            # Use the gene symbol dictionary that we created above for matching symbol with
            # expression value.
            # Stop reading when we encounter the '!sample_table_end' token.
            reader = csv.DictReader(soft_file, delimiter='\t')
            value_count = 0
            for value_line in reader:
                if value_line['ID_REF'] == '!sample_table_end':
                    break
                value_count += 1
                if gene_symbol_field:
                    gene_index = value_line['ID_REF'].upper()
                    try:
                        gene = gene_symbol_dict[gene_index]
                    except KeyError:
                        #print "key error:", gene_index
                        #print value_line
                        continue
                    # Append the gene_index only if it is different than the gene symbol
                    if gene != gene_index:
                        gene_key = '|'.join([sample_id_key, gene_code, gene + '_' + gene_index])
                    else:
                        gene_key = '|'.join([sample_id_key, gene_code, gene])
                    expression = value_line['VALUE']
                else:
                    gene = value_line['ID_REF'].upper()
                    expression = value_line['VALUE']
                    gene_symbol_dict[gene] = gene
                    gene_key = '|'.join([sample_id_key, gene_code, gene])

                # Store the expression value only if there is one otherwise remove it from the gene dictionary so
                # that it doesn't show up in the gene symbol pick list.
                if len(expression) > 0:
                    Put(gene_key, expression)
                else:
                    try:
                        del gene_symbol_dict[value_line['ID_REF'.upper()]]
                    except KeyError:
                        pass

            # Don't read anymore soft files than are listed.
            sample_count += 1
            if sample_count == len(sample_ids):
                break
            else:
                try:
                    line = soft_file.next()
                except StopIteration:
                    line = None

    # No more soft files to process. Begin writing out lists of things such as
    # sample ids, gene symbols, profiles, etc.
    profiles_key = '|'.join([study_key, profiles_code])
    Put(profiles_key, JSONEncoder().encode(profiles))

    platforms_key = '|'.join([profile_key, platforms_code])
    Put(platforms_key, JSONEncoder().encode(platforms))

    sample_ids_key = '|'.join([platform_key, sample_ids_code])
    Put(sample_ids_key, JSONEncoder().encode(sample_ids))
    print "|sample_ids|", len(sample_ids)

    genes_list = []
    # Pair up the gene and gene_index if necessary
    if gene_symbol_field:
        for key, value in gene_symbol_dict.items():
            genes_list.append(value + '_' + key)
    else:
        for key, value in gene_symbol_dict.items():
            genes_list.append(value)

    genes_key = '|'.join([platform_key, genes_code])
    Put(genes_key, JSONEncoder().encode(genes_list))
    print "|gene_symbol_list|", len(genes_list)

    soft_file.close()

# Store the studies

studies_key = studies_code
Put(studies_key, JSONEncoder().encode(studies))

print "Stop", datetime.now()
