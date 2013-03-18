''' Kernelio is IO functions for working with the kernels output by STINGER
using the named results store. We have support for loading sparse vectors from a
file and then mmerging them into a big pandas dataframe.

'''

import pandas as pd

def load_data_histogram(pathfmt, batch):
    """Gets a batch into a series to be used in a the analysis

    Arguments:
    - `pathfmt`:
    - `batch`:
    """
    path = pathfmt % (batch)
    dframe = pd.read_csv(path, names=frozenset(['bin','count']))
    return dframe

def load_sparse_vec(pathfmt, batch, column=-1):
    """reads in a sparse vector that is represented as key value pairs.

    Arguments:
    - `pathfmt`:
    - `batch`:
    - `column`: which column to use if data is a table
    - `names`:
    - `'val']`:
    """
    path = pathfmt % (batch)
    dframe = pd.read_csv(path, index_col=0, header=None)
    dframe = pd.DataFrame(dframe[dframe.columns[column]],columns=[batch])
    nonzero = dframe[dframe!=0].dropna()
    return dframe

def load_batches(pathfmt, batches, column=-1):
    """Load a set of batches into a big data frame

    Arguments:
    - `pathfmt`:
    - `batches`:
    - `column`: which column to use if data is a table
    """
    series = [load_sparse_vec(pathfmt, b, column) for b in batches]
    frame = pd.DataFrame.join(series[0], series[1:], how='outer')
    frame.save
    return frame

def format_hdf_names(data_dir, kernel_name, init_sample, end_sample, stride):
    """ Use this to make the names for the frames within the HDF data store.
    The naming convention might change in the future and all consumers should be
    able to adapt.

    Arguments:
    - `data_dir`:
    - `kernel_name`:
    - `init_sample`:
    - `end_sample`:
    - `stride`:

    """
    store_name = data_dir+kernel_name+'.hf'
    frame_name = '%s.%d.%d.%d' % (kernel_name,
                                  init_sample, end_sample, stride)
    return store_name, frame_name

def load_hdf_table(store_name, frame_name):
    """

    Arguments:
    - `store_name`:
    - `frame_name`:
    """
    store = pd.HDFStore(store_name)
    df = store.select(frame_name)
    store.close()
    return df

def write_hdf_table(store_name, frame_name, frame):
    """ Wraps the storage into the HDFStore from pytables.

    Arguments:
    - `store_name`: the filename for the HDFStore constructor
    - `frame_name`: a string name for the data
    - `frame`: the data
    """
    store = pd.HDFStore(store_name)
    print('writing to hdfstore: %s/%s'% (store_name, frame_name))
    store.append(frame_name, frame)
    store.close()
    print('successfully wrote to hdfstore')

def load_csv_frame(DATA_DIR, KERNEL_NAME, FILENAME, TIMEINDEX):
    """ If HDF is not supported use this to load from a csv or fall back to
    the vectors each in one file.

    Arguments:
    - `DATA_DIR`:
    - `KERNEL_NAME`:
    - `FILENAME`:
    - `TIMEINDEX`:
    """

    try:
        df = pd.read_csv(FILENAME)
        print('sucessfully read file %s'% FILENAME)
        df = df.set_index('0')
        df.columns = df.columns.map(int)
    except:
        print('failed to read file %s\nperforming load from vectors'% FILENAME)
        df = load_batches(DATA_DIR+ KERNEL_NAME+".%d.csv",
                          TIMEINDEX, column=-1)
        print('serializing dataframe as %s' % FILENAME)
        df.to_csv(FILENAME)
    return df
