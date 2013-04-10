import argparse
import sklearn.svm as svm
import sklearn.covariance as sklcov
import pandas as pd
import numpy as np
import scipy
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
FIGURE_PATH = "./figures/"
FIGURE_EXTENSION = 'png'
def get_data(path):
    """Gets a frame containing 
    vtxid, stat1, stat2,..., names
    from disk

    :path: @todo
    :returns: a frame

    """
    df = pd.load(path)
    return df

def save_description(gf):
    """Prints out the description of a frame and saves it to a file"""
    desctab = gf.describe().ix[[(-1.0,'mean'),(-1.0,'std'),(-1.0,'50%'),(1.0,'mean'),(1.0,'std'),(1.0,'50%')]]
    desctab = desctab.T.stack()
    print(desctab)
    tablestr = desctab.to_latex()
    #latexstr = "\begin{tabular}\n%s\n\end{tabular}"%tablestr
    fp = open(DATA_DIR+'svmtable.tex', 'w')
    fp.write(tablestr)
    fp.close()
    return tablestr

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
        svdetector = svm.OneClassSVM(kernel='rbf',nu=.10,gamma=.3, cache_size=2000)
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

    if args.plot:
        print('plotting')
        active_cols = labelf.columns[:-2]
        n = len(active_cols)
        #fig = plt.figure(figsize=(5,7))
        plots = plt.subplots(n,1,0)
        axes  = plots[1]
        fig   = plots[0]
        fig.subplots_adjust(.09,.06,.94,.94,.2,.5)
        for i in range(n):
            ax = axes[i]
            print(ax)
            col = labelf.columns[i]
            dcol = labelf[col]
            pos = dcol[labelf.pred ==1]
            neg = dcol[labelf.pred ==-1]
            plt.sca(ax)
            plt.hist([pos,neg], normed=True, color=['b','r'])
            ax.set_title(col)
        fig.savefig("%ssvm-outliers-hist.%s"% (FIGURE_PATH, FIGURE_EXTENSION))
        axes = pd.scatter_matrix(labelf[active_cols],marker='.', c=-1.5*labelf.pred,)# marker=labelf.pred+1)
        fig = plt.gcf()
        fig.subplots_adjust(.10,.12,.90,.90,.0,.0)
        print('saving')
        fig.savefig("%ssvm-outliers.%s"% (FIGURE_PATH, FIGURE_EXTENSION))
        #fig, axes = plt.subplots(1,4)
        #axes[0].scatter(x=labelf[col1], y=labelf[col2],
        #             s=10, alpha=.3, c=-1*labelf.pred)
        #axes[1].scatter(x=labelf[col3], y=labelf[col2],
        #             s=10, alpha=.3, c=-1*labelf.pred)
        #axes[2].scatter(x=labelf[col1], y=labelf[col3],
        #             s=10, alpha=.3, c=-1*labelf.pred)
        #axes[3].scatter(x=labelf[col3], y=labelf[col4],
        #             s=10, alpha=.3, c=-1*labelf.pred)
        #print('plotting')
    anoms = labelf[labelf.pred == -1]
    #anoms = anoms.sort(columns=[col1, col2])
    print(anoms.head(20))
    return labelf, anoms, model

if __name__ == '__main__':
    args = get_args()
    DATA_DIR = u'./'
    FILENAME = u'ccfeatures.data'
    df = get_data(DATA_DIR+FILENAME).dropna()
    print(df)
    df = (df - df.mean())/df.std()
    print(df)
    #smallf = df.ix[::5][['meanvalue','varvalue','meanderiv',]]
    smallf = df.ix[::5][['meanvalue','varvalue','meanderiv','varderiv',]]#'maxcc','maxtri']]
    data_mat = smallf  
    print(df.head(5))
    labelf, anoms, model = main()
    gf = labelf[labelf.columns[:-1]].groupby('pred')
    print(gf.describe())
    save_description(gf)
