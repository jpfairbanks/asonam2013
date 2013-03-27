import argparse
import sklearn.svm as svm
import sklearn.covariance as sklcov
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

def skeptics(df, name_series, search_string, decision_func):
    """Looks for vertices based on an initial substring of their twitter handle
    Arguments:
    - `df`: 
    - `name_series`: the series from which to draw the names
    - `search_string`: the prefix that we are looking for
    - `decision_func`: the decision_function to use for prediction

    """
    skeps = name_series[name_series.str.startswith(search_string)].index
    skf = df.ix[skeps]
    print(skf)
    sk_scores = pd.Series(decision_func(skf[['mu','sigma']].values), index=skeps, name='score')
    skff = skf.join(sk_scores)
    print(skff)

def get_data(path):
    """Gets a frame containing 
    vtxid, stat1, stat2,..., names
    from disk

    :path: @todo
    :returns: a frame

    """
    df = pd.read_csv(path, index_col=0)
    if args.triangles:
        df = df[['mu','sigma','1']]
    return df

def get_args():
    """Use argparse to make this a command line tool
    :returns: @todo

    add new options here
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-b', '--betweenness',
                        help='use betweenness data for this',
                        action='store_true')
    parser.add_argument('-t', '--triangles',
                        help='use triangle data for this',
                        action='store_true')
    parser.add_argument('-e', '--ellipse',
                        help='use an elliptical envelop to identify the center',
                        action='store_true')

    parser.add_argument('-s', '--svm',
                        help='use a support vector machine with one class to find outliers',
                        action='store_true')
    parser.add_argument('-c', '--score',
                        help='Not Implemented! compute the score values as well as binary labels. These are not comparable across type',
                        action='store_true')
    parser.add_argument('-p', '--plot',
                        help='show a scatter plot with colors for the anomalies',
                        action='store_true')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    if args.betweenness:
        DATA_DIR = u'../data/post_process/'
        FILENAME = u'whitened_BC_stats.csv'
    if args.triangles:
        DATA_DIR = u'/scratch/jfairbanks/tri_bc_sandy/post_processing/'
        FILENAME = u'whitened_tri_stats.csv'
    df = get_data(DATA_DIR+FILENAME)
    smallf = df.ix[::1]
    data_mat = smallf[['mu','sigma']]
    names = df['1']#.str.lower()
    print(df.head(5))
    if args.ellipse:
        envelope = sklcov.EllipticEnvelope(assume_centered=True, contamination=.02)
        print('fitting')
        envelope.fit(data_mat.values)
        envelope.correct_covariance(data_mat.values)
        print('fitting done')

        print('predicting')
        scores = pd.Series(envelope.decision_function(data_mat), name='score')
        print('predicting done')
        flags = np.sign(scores)
        flags.name = 'pred'
        skl_out = pd.DataFrame({s.name:s for s in [scores, flags]})
    
    if args.svm:
        svdetector = svm.OneClassSVM(kernel='rbf', cache_size=2000)
        print('fitting')
        svdetector.fit(data_mat)
        print('fitting done for SVM based detector')
        print('predicting')
        scores = pd.Series(np.zeros((data_mat.index.shape)), index=data_mat.index, name='score')
        flags  = pd.Series(svdetector.predict(data_mat), index=data_mat.index, name='pred')
        skl_out = pd.DataFrame({s.name:s for s in [scores,flags]})
        print('predicting done')

    print('joining')
    labelf = smallf
    labelf = labelf.join(skl_out)
    print('joining')

    if args.plot:
        print('plotting')
        fig, ax = plt.subplots(1,1)
        ax.scatter(x=labelf.mu, y=labelf.sigma,
                     s=10, alpha=.7, c=-1*labelf.pred)
        print('plotting')
    anoms = labelf[labelf.pred == -1]
    anoms = anoms.sort(columns=['score', 'pred'])
    print(anoms.head(20))
    #skeptics(df, names, 'skep', envelope.decision_function)
