import os, re
import configparser
import pandas as pd
from ROOT import TFile, TH1F, TF1, gStyle

# plot statistic box in fit canvas
gStyle.SetOptFit(0)
gStyle.SetOptStat(0)


def find(dtype, path):
    '''
    Returns the files pathes for .root files under the given path. 
    dtype can be anything, but for best results I recommend: 'inter' or 
    'physics' for intermediate and physics files.

    Parameters
    ----------
    dtype : string
        In principle it can be anything. Although, 'inter' and 'physics' will give
        the best results.

    path : string
        Path to the data.

    Returns
    -------
    list
        List with found file pathes + rootPath as the last entry.
    '''
    ofiles = []
    for root, dirs, files in os.walk(path):
        files.sort()
        if not root.find(dtype) > -1: continue
        for file in files:
            if not file.endswith('.root'): continue
            ofiles.append(os.path.join(root,file))
    print('{} files found.'.format(len(ofiles)))
    ofiles.append(path)
    return ofiles


def integrate(rootFilePathes, pedestalSubraction=True):
    '''
    Integrates the signal collected for MIP calibration. Integration happens the 
    following: I(L/2,L)-I(0,L/2) with L length of the waveform..
    This ensures pedestal correction. Offset correction is also applied.
    
    Parameters
    ----------
    rootFilePathes : list
        List with .root files pathes, e.g. output of mipcali.find function.
        The last entry is ignored.

    pedestalSubraction : bool
        Defaul it 'True'. The pedestal subtraction is done in a loop and thus 
        highly inefficient. Turn it off whenever possible.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one column per module, entries is the pedestal subracted
        integral.
    '''

    # read in the channel layout from the channels.info file
    ddict = {}
    try:
        ddict = _getChannelID(os.path.join(rootFilePathes[-1], 'channels.info'))
    except OSError as e:
        print(e)
        print('The channels.info file is mandatory for the channel layout!')
        return pd.DataFrame()

    # predefine a dict which holds the data for each module 
    ddata = {}
    for pico in ddict.keys():
        for channel in ddict[pico].keys():
            ddata[ddict[pico][channel]] = []

    # loop through the files and integrate
    for file in rootFilePathes[:-1]:
        tfile = TFile(file)
        subRunNumber = int(re.findall('\d+', file)[-1])
        for pico in ddict.keys():
            for channel in ddict[pico].keys():
                moduleNumber = ddict[pico][channel]
                th1 = 0
                try:
                    th1 = tfile.Get("{}_{}_{}".format(pico,channel,subRunNumber))
                    th1.GetBinContent(1)
                except AttributeError as excep:
                    try:
                        th1 = tfile.Get("{}_{}".format(pico,channel))
                        th1.GetBinContent(1)
                    except AttributeError as exc:
                        print(exc)
                # find the offset by averaging over the pedestal part
                offset = th1.Integral(1,int(round(th1.GetNbinsX()/2)))/(th1.GetNbinsX()/2.)
                # apply offset correction
                for i in range(1,th1.GetNbinsX()+1):
                    th1.AddBinContent(i, -offset)
                #integrate and add integral to data container
                ddata[moduleNumber].append(th1.Integral(int(round(th1.GetNbinsX()/2)), th1.GetNbinsX()) - th1.Integral(1,int(round(th1.GetNbinsX()/2))))

    return pd.DataFrame(ddata)
            

def pdStoTH1F(series, bins = 0, name = 'th1', lower=None, upper=None):
    '''
    Takes a pandas.Series and makes a histogram out of it.
    Example: th1f = pdStoTH1F(df.B)

    Parameters
    ----------
    series : pandas.Series
        Data

    bins : int
        Default = 0 -> Difference between highest and smallest value in integers.

    name : string
        Name the th1f ouput. This is ROOT specific stuff.

    lower : int
        Is the lower bound of the returning histogram.
        If not defined: lower bound is found by minimum of series.

    upper : int
        Is the upper bound of the returning histogram.
        If not defined: upper bound is found by maximum of series.

    Returns
    -------
    ROOT.th1f
    '''
    if lower == None:
        lower = round(series.min())
    if upper == None:
        upper = round(series.max())

    if bins == 0:
        bins = int(upper-lower)

    tf1 = TH1F(name, name, bins, lower, upper)
    for i in series.index:
        tf1.Fill(series[i])

    return tf1


def fitGausMax(th1f, sigmaRange = 1.):
    '''
    Returns the fitfunction of a gaus likelihood fit around the maximum with starting parameters:
    c       = th1f.GetMaximum()
    mean    = th1f.GetBinCenter(th1f.GetMaximumBin())
    sigma   = th1f.GetStdDev())


    Parameters
    ----------
    th1f : ROOT.TH1
        Data should be gaus like distributed.

    sigmaRange : float
        The fit is performed +- sigmaRange*sigma around the maximum of the 
        histogram.

    Returns
    -------
    ROOT.tf1
        Fitfunction.

    dict
        Fit results.

    '''

    fitf = TF1('f1', 'gaus',th1f.GetBinCenter(th1f.GetMaximumBin())- sigmaRange * th1f.GetStdDev(),th1f.GetBinCenter(th1f.GetMaximumBin())+th1f.GetStdDev()* sigmaRange)
    fitf.SetParameters(th1f.GetMaximum(),th1f.GetBinCenter(th1f.GetMaximumBin()), th1f.GetStdDev())
    th1f.Fit('f1', 'LR')
    
    fitResults = {  'amplitude' : fitf.GetParameter(0),
                    'mean' : fitf.GetParameter(1),
                    'mean_err' : fitf.GetParError(1),
                    'sigma' : fitf.GetParameter(2), 
                    'sigma_err' : fitf.GetParError(2),
                    'ndf' : fitf.GetNDF(),
                    'chisquare' : fitf.GetChisquare()}


    return fitf, fitResults


_dataChannel1 = 'B'
_dataChannel2 = 'C'

def defDataChannels( dataChannel1, dataChannel2 ):
    '''
    If you use different data channels than b & c, change that here globally.
    Options are two out of a,b,c,d.
    '''
    _dataChannel1 = dataChannel1
    _dataChannel2 = dataChannel2


def _getChannelID( path, dataChannel1 = _dataChannel1, dataChannel2 = _dataChannel2 ):
    '''
    Reads module numbers/names from the channel.ini file for correct naming 
    of the data.
    
    Returns
    -------
    dict
        It is a dict of a dict.

    Exceptions
    ----------
    OSError if file not found.
    '''
    if not os.path.isfile(path):
        raise OSError(path + ' not found.')
    parser = configparser.ConfigParser()
    parser.read(path)

    ddict = {}
    for section in parser.sections():
        sdict = {}
        for (key, val) in parser.items(section):
            if key[-1] == dataChannel1.lower():
                sdict[dataChannel1] = val
            elif key[-1] == dataChannel2.lower():
                sdict[dataChannel2] = val
        ddict[section] = sdict

    return ddict



