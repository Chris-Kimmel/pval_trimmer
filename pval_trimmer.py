'''
pval_trimmer.py
Trim the 5' and 3' ends of wide-form p-value CSV files.

Written for the Kim Lab (Aug 9, 2021)

The script takes six arguments from the command line. In order, they are:
    - trim_from_5prime (Trim this many bases from the 5' end of each read)
    - trim_from_3prime (Trim this many bases from the 3' end of each read)
    - exemption_boundary_5prime (Do not trim a read's 5' end if the 5' end of
        the untrimmed read is upstream of this position)
    - exemption_boundary_3prime (Do not trim a read's 3' end if the 3' end of
        the untrimmed read is downstream of this position)
    - filepath_to_read (include .csv extension if applicable)
    - filepath_to_write (include .csv extension if applicable)

Here is an example of correct command-line usage:
python3 pval_trimmer.py 50 50 0 9179 untrimmed.csv trimmed.csv

Chris Kimmel
chris.kimmel@live.com
'''

################################################################################
############################# imports and pragmas ##############################
################################################################################

# pylint: disable=invalid-name,line-too-long

import sys
import pandas as pd


################################################################################
################################## constants ###################################
################################################################################

'''A (p-value, read_id, position) is kept if and only if this query is TRUE for
it. The query is listed here mainly to clarify what the script does. (If you
change the query and it doesn't work, it'll be hard to debug.)'''
QUERY = \
            "(dist_from_end_5prime >= @trim_from_5prime "\
            "or end_5prime_0b < @exemption_boundary_5prime) "\
            "and "\
            "(dist_from_end_3prime >= @trim_from_3prime "\
            "or end_3prime_0b > @exemption_boundary_3prime) "


################################################################################
############################ command-line interface ############################
################################################################################

trim_from_5prime = int(sys.argv[1])
trim_from_3prime = int(sys.argv[2])
exemption_boundary_5prime = int(sys.argv[3])
exemption_boundary_3prime = int(sys.argv[4])
filepath_to_read = sys.argv[5]
filepath_to_write = sys.argv[6]


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

wide = load_csv(filepath_to_read)

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

widify(trimmed).to_csv(filepath_to_write)
