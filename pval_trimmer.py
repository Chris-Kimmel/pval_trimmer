'''
pval_trimmer.py
Trim the 5' and 3' ends of wide-form p-value CSV files.

Written for the Kim Lab (Aug 9, 2021)

The user should go to "constants for the user to customize" (a section in the
code below) and adjust the constants before running this script.

Chris Kimmel
chris.kimmel@live.com
'''

################################################################################
############################# imports and pragmas ##############################
################################################################################

# pylint: disable=invalid-name,line-too-long

import pandas as pd


################################################################################
##################### constants for the user to customize ######################
################################################################################

# full filepaths (including .csv extensions if applicable)
FILEPATH_TO_READ = '/fs/project/PAS1405/GabbyLee/project/m6A_modif/WT_cellular/23456_WT_cellular.csv'
FILEPATH_TO_WRITE = 'output.csv'

TRIM_FROM_5PRIME = 50 # Trim this many bases from the 5' end of each read
TRIM_FROM_3PRIME = 50 # Trim this many bases from the 3' end of each read

# Exempt reads from 5'-end trimming when they start upstream of
# EXEMPTION_BOUNDARY_5PRIME. Likewise, exempt a read from *3'-end* trimming if
# the read starts *downstream* of the 3' boundary.
EXEMPTION_BOUNDARY_5PRIME = 0
EXEMPTION_BOUNDARY_3PRIME = 9179

'''A (p-value, read_id, position) is kept if and only if this query is TRUE for
it. The query is listed here mainly to clarify what the script does. (If you
change the query and it doesn't work, it'll be hard to debug.)'''
QUERY = \
            "(dist_from_end_5prime >= @TRIM_FROM_5PRIME "\
            "or end_5prime_0b < @EXEMPTION_BOUNDARY_5PRIME) "\
            "and "\
            "(dist_from_end_3prime >= @TRIM_FROM_3PRIME "\
            "or end_3prime_0b > @EXEMPTION_BOUNDARY_3PRIME) "


################################################################################
################################# subroutines ##################################
################################################################################

def load_csv(filepath):
    '''Load per-read stats from a CSV file into a Pandas DataFrame'''
    retval = (
        pd.read_csv(filepath, header=0, index_col=0)
        .rename_axis('pos_0b', axis=1)
        .rename_axis('read_id')
    )
    retval.columns = retval.columns.astype(int)
    return retval

def longify(df):
    '''Convert dataframe output of load_csv to a long format'''
    return df.stack().rename('pval').dropna().reset_index()

def widify(df):
    '''Undo longify'''
    return df.set_index(['read_id', 'pos_0b']).unstack().droplevel(0, axis=1)


################################################################################
################################## load data ###################################
################################################################################

wide = load_csv(FILEPATH_TO_READ)

if not wide.index.is_unique:
    raise NotImplementedError("There is a duplicate read ID in the given file. "
        "This script requires unique read IDs.")

untrimmed = longify(wide)
del wide


################################################################################
############################## compute read ends ###############################
################################################################################

end = (
    untrimmed.groupby('read_id')
    .agg(
        end_5prime_0b=('pos_0b', min),
        end_3prime_0b=('pos_0b', max),
    )
    .reset_index()
) # should have only three columns: 'read_id', 'end_5prime_0b', 'end_3prime_0b'


################################################################################
################################## trim reads ##################################
################################################################################

trimmed = (
    untrimmed
    .merge(end, on='read_id', how='inner', validate='many_to_one')
    .assign(
        dist_from_end_5prime=lambda x: x['pos_0b'] - x['end_5prime_0b'],
        dist_from_end_3prime=lambda x: x['end_3prime_0b'] - x['pos_0b'],
    )
    .query(QUERY)
    .loc[:, ['read_id', 'pos_0b', 'pval']]
)


################################################################################
################################## write data ##################################
################################################################################

widify(trimmed).to_csv(FILEPATH_TO_WRITE)
