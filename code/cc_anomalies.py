import argparse
import sklearn.svm as svm
import sklearn.covariance as sklcov
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

def get_data(path):
    """Gets a frame containing 
    vtxid, stat1, stat2,..., names
    from disk

    :path: @todo
    :returns: a frame

    """
    df = pd.load(path)
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
def main():
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
        model = envelope

    if args.svm:
        #nu determines the amount of outliers that we find we want ten percent
        svdetector = svm.OneClassSVM(kernel='rbf',nu=.1,gamma=.3, cache_size=2000)
        print('fitting')
        svdetector.fit(data_mat)
        print('fitting done for SVM based detector')
        print('predicting')
        scores = pd.Series(np.zeros((data_mat.index.shape)), index=data_mat.index, name='score')
        flags  = pd.Series(svdetector.predict(data_mat), index=data_mat.index, name='pred')
        skl_out = pd.DataFrame({s.name:s for s in [scores,flags]})
        print('predicting done')
        model = svdetector 

    print('joining')
    labelf = smallf
    labelf = labelf.join(skl_out)
    print('joining')

    col1 = 'count'
    col2 = 'maxcc'
    col3 = 'mean'
    col4 = 'maxtri'
    if args.plot:
        print('plotting')
        fig, axes = plt.subplots(1,4)
        axes[0].scatter(x=labelf[col1], y=labelf[col2],
                     s=10, alpha=.3, c=-1*labelf.pred)
        axes[1].scatter(x=labelf[col3], y=labelf[col2],
                     s=10, alpha=.3, c=-1*labelf.pred)
        axes[2].scatter(x=labelf[col1], y=labelf[col3],
                     s=10, alpha=.3, c=-1*labelf.pred)
        axes[3].scatter(x=labelf[col3], y=labelf[col4],
                     s=10, alpha=.3, c=-1*labelf.pred)
        print('plotting')
    anoms = labelf[labelf.pred == -1]
    anoms = anoms.sort(columns=[col1, col2])
    print(anoms.head(20))
    return labelf, anoms, model
if __name__ == '__main__':
    args = get_args()
    DATA_DIR = u'./'
    FILENAME = u'ccfeatures.data'
    df = get_data(DATA_DIR+FILENAME).dropna()
    smallf = df.ix[::5]
    data_mat = smallf[['mean','count','maxcc','maxtri']]
    print(df.head(5))
    labelf, anoms, model = main()
    gf = labelf.groupby('pred')
    print(gf.describe())
    #pd.scatter_matrix(labelf, color=labelf.pred)
    pd.scatter_matrix(labelf[labelf.pred == 1], color='b')
    pd.scatter_matrix(labelf[labelf.pred == -1], color='r')
