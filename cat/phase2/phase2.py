from . import _global as bg

import numbers as lnumbers
import math
import configparser, os, re, copy, platform

import pandas as pd
import pytz as tz
from datetime import datetime, timezone
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from tqdm.auto import tqdm as tqdm
from ROOT import TFile, TH1, TH1I, TH1F, TCanvas

from pathlib import Path
import os


def setPhase3(bool = False, direction=None):
    if bool:
        bg.picos = ['Z0', 'Z1', 'Z2', 'Z3']
        bg.phase3 = True

    if direction == None:
        raise IOError('No direction (fwd or bwd) specified!')
    elif isinstance(direction, str):
        if direction.find('fwd') > -1 or direction.find('bwd') > -1:
            bg.phase3_direction = direction
        else:
            raise IOError('Direction not known!')

def setDataPath(path=None):
    '''
    Sets the path to the data on the NAS.
    The default path is $HOME/NAS_futDet/claws/phase{2,3}/raw/

    If this does not apply to you, please define it yourself with the path variable.
    '''

    if path == None:
        if not bg.phase3:
            home = str(Path.home())
            bg.path_to_data = os.path.join(home,"NAS_futDet/claws/phase2/raw")
            if os.path.isdir(bg.path_to_data):
                print('Data path is set to {}'.format(bg.path_to_data))
                return
        else:
            home = str(Path.home())
            bg.path_to_data = os.path.join(home,"NAS_futDet/claws/phase3/data", bg.phase3_direction)
            if os.path.isdir(bg.path_to_data):
                print('Data path is set to {}'.format(bg.path_to_data))
                return

    else:
        if os.path.isdir(path):
            bg.path_to_data = path
            print('Data path is set to {}'.format(bg.path_to_data))
            return

    raise OSError('Data path could not be found! Please define!')


def getDataPath():
    '''
    Returns the current data path.
    '''
    return bg.path_to_data


def setIndexFilePath(path):
    bg.path_to_indexFile = path
    print("Index file path is set to {}".format(bg.path_to_indexFile))


def getIndexFilePath():
    return bg.path_to_indexFile



def getAllDatesOfData( pathToData=None):
    '''
    Lists all dates where data is avaiable.

    Parameters
    ----------
    pathToData : string
        path to the data, up to the level of the JJJJ-MM-DD dir structure.
        !Should be an absolute path!
    Returns
    -------
    dates : list of strings
    '''
    if pathToData == None:
        pathToData = getDataPath()

    listOfAllDates = os.listdir( pathToData )
    listOfAllDates = [ i for i in listOfAllDates if len(i)==len('JJJJ-MM-DD') ]
    listOfAllDates.sort()

    return listOfAllDates


def provideListOfDates(
    startDate,
    endDate,
    pathToData=None
    ):
    '''
    Creates a list of dates where data is available.

    Parameters
    ----------

    startDate : string
        should be given in the form of JJJJ-MM-DD

    endDate : string
        should be given in the form of JJJJ-MM-DD

    pathToData : string (optional)
        path to the data, up to the level of the JJJJ-MM-DD dir
        structure.

        !Should be an absolute path!

        Default is getDataPath()

    Returns
    -------
    dates : list of strings
    '''
    if pathToData == None:
        pathToData = getDataPath()

    dates = []
    files = []
    a = os.listdir( pathToData )
    b = [i for i in a if len(i)==len(startDate) ]
    b.sort()
    dates = b[ b.index( startDate ):b.index( endDate )+1 ]

    return dates


def findFiles(
        dates,
        pathToData=None,
        defRundNumRange=('',''),
        dType="online"
        ):
    '''
    Creates a list of paths to all .root-files of the specified time
    periode and data type.

    Parameters
    ----------

    dates : list of str of the form ['JJJJ-MM-DD', 'JJJJ-MM-DD', ...],
        can be accurate by using the function "provideListOfDates()"

    pathToData : string
        path to the data, up to the level of the JJJJ-MM-DD dir
        structure.

        !Should be an absolute path!

        The default is defined by setDefaultDataPath()

    defRundNumRange : Tuple
        By specifying a tuple of the form (start run #, end run #) one
        can have a look at a certain set of runs within the specified
        time range. Any entry which can be converted to an integer is
        acceptable

    dType : string

    Returns
    -------
    files : list of strings (file paths)
    '''
    if pathToData == None:
        pathToData = getDataPath()

    files = []

    for date in dates:
        pathToRuns = os.path.join(pathToData, date)
        if not os.path.isdir(pathToRuns):
            continue
        runsRaw = os.listdir( pathToRuns )
        runNumStart = int( runsRaw[0].split('-')[1] )
        runNumEnd = int(  runsRaw[-1].split('-')[1] )
        runs = []

        if not defRundNumRange == ('',''):
            if ( runNumStart <= int(defRundNumRange[0])
                or runNumEnd >= int(defRundNumRange[1]) ):
                for i in runsRaw:
                    if (int(defRundNumRange[0])
                        <= int( i.split('-')[1])
                        <= int(defRundNumRange[1]) ):
                        runs.append(i)
            else:
                runs = runsRaw
        else:
            runs = runsRaw

        if dType.find('inter') > -1:
            dType = 'raw/intermediate'
        elif dType.find('physics') > -1:
            dType = 'raw/physics'

        for i in runs:
            try:
                pathTof = os.path.join(pathToRuns, i, dType)
                os.path.isdir(pathTof)
            except Exception as e:
                print(e)
                continue
            try:
                filenames = os.listdir( pathTof )
            except Exception as exce:
                print(exce)
                continue
            for j in filenames:
                if not j.endswith('.root'):
                    continue
                files.append(os.path.join(pathTof, j))

    return files


def unixToHumanJSTTime(ts):
    '''
    Transforms the unix timestamp taken at KEK (JST) to human readable
    time format in JST timezone.

    !!! Bew ware: If this is used at KEK things might go wrong. !!!
                       !!! Please verify !!!

    Parameters
    ----------
    ts :

    Returns
    ----------
    The japanese standrad time in the format of YYYY-MM-DD HH:MM:SS
    '''
    utc_ts = datetime.fromtimestamp(ts)
    jst = tz.timezone('Asia/Tokyo')
    cet = tz.timezone('Europe/Berlin')
    jst_time = utc_ts + jst.utcoffset(utc_ts) - cet.utcoffset(utc_ts)
    return jst_time.strftime("%Y-%m-%d %H:%M:%S.%f")


def th1XScale( th1, factor ):
    '''
    Rescales the th1 xAxis by the given factor. Remember to change the
    xAxis legend!

    Parameters
    ----------
    th1 : ROOT.TH1
        TH1 to be changend.

    factor : float
        Scaling factor.

    Returns
    -------
    ROOT.TH1
    '''
    nbins = th1.GetNbinsX()
    xlow = th1.GetBinLowEdge(1)
    xup_old = th1.GetXaxis().GetBinUpEdge(th1.GetNbinsX()) +\
                th1.GetBinWidth(th1.GetNbinsX())
    xup_new = xup_old * factor
    th1.SetBins(nbins, xlow, xup_new)

    return th1


def makeAvg( dfPaths,
             picos,
             channels
            ):
    '''
    Averages the waveforms over the given events. Per event all given
    locations and channels are summed up.

    Parameters
    ----------
    dfPaths : pandas.DataFrame
        Data frame with index and paths for events only.

    picos : list of str
        List with all locations/picos which should be considered.
        Possible locations:
            - 'Top_Forward'
            - 'Top_Backward'
            - 'Bottom_Forward'
            - 'Bottom_Backward'
            e.g. picos = ['Top_Forward', 'Bottom_Backward']

    channels : list of str
        List with all channels which should be considered. Possible
        channels:
            - 'A'
            - 'B'
            - 'C'
            - 'D'
            e.g. channels = ['A', 'C']

    Returns
    -------
    th1 : ROOT.TH1I
        average waveform
    th2 : ROOT.TH1I
        signal heigth histogram.
    '''

    data = dfPaths.drop_duplicates()

    highestVal = -1
    th1 = 0
    th2 = 0
    thh = TH1F('Signal heigth histogram',
                'Signal height histogram',
                1000,
                0,
                1000
                )
    isFirst = True
    counter = 0
    for i, c1 in data.iterrows():
        firstEvent = True
        tfile = TFile(c1[0])
        #print('tfile: {0} \t {1}'.format(tfile, c1[0]))
        if not tfile: raise Exception( 'File not found: {}'.format(
                                                                c1[0]
                                                                )
                                     )

        for pico in picos:
            for channel in channels:
                dth1 = tfile.Get(pico
                                + '_'
                                + channel
                                + '_'
                                + bg.items_channel['reco']
                                )
                dmax = dth1.GetBinContent(dth1.GetMaximumBin())
                dmin = dth1.GetMinimum()
                if dmax > 1e6:
                    print('{0}_{1} max bin content too high: {2}'.format(pico,channel,dmax))
                    continue
                if dmin < -1e6:
                    print('{0}_{1} min bin content too low: {2}'.format(pico,channel,dmin))
                    continue

                if dmax > highestVal:
                    highestVal = dmax

                if firstEvent:
                    th2 = copy.deepcopy(dth1)
                    firstEvent = False
                else:
                    th2.Add(dth1)

                for i in range(dth1.GetNbinsX()):
                    dval = dth1.GetBinContent(i+1)
                    if dval > 1e-6:
                        thh.Fill(dval)

                counter += 1

        if isFirst:
            th1 = copy.deepcopy(th2)
            isFirst = False
        else:
            th1.Add(th2)
        tfile.Close()

    # create proper th1 title
    title = "AvgWf_{0}wfs_".format(counter)
    for i in range(len(picos)):
        if picos[i]=='Top_Forward':
            title += 'TF'
        elif picos[i] == 'Top_Backward':
            title += 'TB'
        elif picos[i] == 'Bottom_Forward':
            title += 'BF'
        elif picos[i] == 'Bottom_Backward':
            title += 'BB'
        if i >= len(picos)-1:
            title += '_'
            continue
        title += '-'

    for channel in channels:
        title += channel

    th1.SetTitle(title)
    print('# of overlaid events:'
         +' {0}\nmaxVal: {3}\nPicos {1}\nChannels {2}'.format(counter,
                                                            picos,
                                                            channels,
                                                            highestVal)
         )
    #th1.Scale(1./counter)
    # add waveforms over all channels, then average over runs
    th1.Scale(1./len(data))


    return th1, thh


def makeFFT(th1Wf):
    '''
    Simple Fourier Transformation approximation:
    This function determines the distance between all peaks and weigths
    that with the product of the heights of the two peaks. Therefore, the
    output unit on the y axis is MIPs x MIPs / bin width.

    Parameters
    ----------
    th1Wf : ROOT.TH1
        Waveform to analyse.

    Returns
    -------
    ROOT.TH1F
        Return the 'low cost' fast fourier transformation of the input
        waveform.
    '''
    # create dummy container filled with zeros
    ######
    ddata = [0. for bin in range(th1Wf.GetNbinsX())]

    # loop through th1Wf. Measure distance between two entries and fill
    # that into  the data container weigthed with the product of the two
    # entries
    ######
    for bin in range(1,th1Wf.GetNbinsX()-1):
        y0 = th1Wf.GetBinContent(bin)
        for sbin in range(bin+1,th1Wf.GetNbinsX()):
            y1 = th1Wf.GetBinContent(sbin)
            ddata[sbin-bin] += y0*y1

    ####################################################################
    # transform the data from the dummy container to a th1Wf
    ####################################################################
    # create new th1Wff for fft content
    ######
    th2 = TH1F('tff',
            'tff',
            th1Wf.GetNbinsX(),
            0,
            th1Wf.GetBinCenter(th1Wf.GetNbinsX())
            )
    th2.SetXTitle('Time [ms]')
    th2.SetYTitle('MIPs x MIPs / {0}us [x1000]'.format('8'))
    th2.SetTitle('{0}_{1}'.format(th1Wf.GetTitle(),'FFT'))
    for bin in range(len(ddata)):
        th2.SetBinContent(bin,ddata[bin])
    th2.Scale(1./1000.)

    return th2


def saveCanvas( c1, name, info = ""):
    '''
    Saves the TCanvas in several formats including the info file;
    !!! 2018-06-01 jpg does not work because of missing codecs !!!

    '''
    c1.SaveAs(name + ".root")
    c1.SaveAs(name + ".pdf")
    c1.SaveAs(name + ".C")
    c1.SaveAs(name + ".png")
#    c1.SaveAs(name + ".jpg")
    info = name + "\n" + info
    with open(name+".info",'w') as f:
        f.write(info)


def saveFig(fig,
            ax,
            title="",
            pdf = True,
            jpg = True,
            png = True,
            dpi_fix = 150
            ):
    '''
    Saves a Figure as PDF, JPG and PNG with a specified quality

    Parameters
    ----------
    fig :
        Path to info*.ini file
    ax :
        matplotlip output ....
    optional - title : string
                    string with no white spaces or other path forbidden
                    asci symboles
    optional - pdf : bool
    optional - jpg : bool
    optional - png : bool

    optional - dpi_fix : integer

    Returns
    -------
    creates a new file for each data type
    '''
    dtitle = 0
    if title == "":
        dtitle = ax.get_title()
    else:
        dtitle = title
    if pdf:
        fig.savefig('{0}.pdf'.format(dtitle), dpi=dpi_fix)
    if jpg:
        fig.savefig('{0}.jpg'.format(dtitle), dpi=dpi_fix)
    if png:
        fig.savefig('{0}.png'.format(dtitle), dpi=dpi_fix)


def TH1ToDf(th1, name=None):
    '''
    Converts a ROOT.TH1 to a pandas.DataFrame.

    Parameters
    ----------
    th1 : ROOT.TH1
        Root histogram.

    name : str
        Give the y axis a proper naming.

    Returns
    -------
    pandas.DataFrame
        The index is the x axis, bin content is column 'y' or 'name'.
    '''
    content_name = 'y' if name is None else name
    d = { 'x' : [], content_name : []}
    for i in range(1,th1.GetNbinsX()):
        d['x'].append(th1.GetBinCenter(i))
        d[content_name].append(th1.GetBinContent(i))

    return pd.DataFrame(d).set_index('x')


def makeAveragedWf( rawWFdatafram ):
    '''
    Calculate average CLAWS waveform

    Parameters
    ----------
    rawWFdatafram : pandas data frame
        data frame consisting of waveforms for one event each.

    Returns
    -------
    avWFdatafram : pandas data frame
        returns the averaged claws wf and it's cumulatative sum.
    '''
    avWFdatafram = rawWFdatafram
    avWFdatafram['cumsum'] = avWFdatafram.cumsum()/avWFdatafram.y.sum()*100
    avWFdatafram.index = avWFdatafram.index/1000.

    return avWFdatafram


def makePathToInfo( filepath ):
    '''
    Helper function for readEvent function, only for early phase 2 data:
    Find info*.ini to corresponding online file

    Parameters
    ----------
    filepath : str
        Path to online*.root file

    Returns
    -------
    str
        Returns the info*.ini file to the given online*.root file.
    '''
    filepath = filepath.replace('online', 'raw/physics',1)
    filepath = filepath.replace('online', 'info')
    filepath = filepath.replace('root','ini')

    return filepath


def makePathToRaw( filepath ):
    '''
    Helper function for readEvent function, only for early phase 2 data:
    Find physics*.root to corresponding online file

    Parameters
    ----------
    filepath  : str
        Path to online*.root file

    Returns
    -------
    str
        Returns the physics*.root file to the given online*.root file.
        Note: There is no guarantee that the physics file exists! The returned
        path is a hypothetical one.
    '''
    filepath = filepath.replace('online', 'raw/physics',1)
    filepath = filepath.replace('online', 'physics')

    return filepath


def getTimestamp( pIni ):
    '''
    Helper function for readEvent function, only for early phase 2 data:
    Optain timestamp for early phase 2 data
    Necessary since the timestamp is not available in some root files at
    the beginning of phase2.

    Parameters
    ----------
    pIni : str
        Path to info*.ini file

    Returns
    -------
    timestamp : str
        Returns the timestamp written in the info*.ini file.
    '''
    if not os.path.isfile(pIni):
        raise Exception("Ini-file not found: " + pIni)
    config = configparser.ConfigParser()
    config.read(pIni)
    timestamp = 0
    counter = 0
    for key in config[config.sections()[0]]:
        if key.endswith('timestamp'):
            timestamp += config.getfloat(config.sections()[0], key)
            counter += 1
    timestamp /= counter

    return timestamp

def readEvent( pathToRootFile, wf_maxLength=None, wf_minLength=None , correctCirculatingBackground=False):
    '''
    Reads one event/file
    Is dependent on the getTimestamp and makePathToInfo helper functions
    for early phase 2 data

    Parameters
    ----------
    pathToRootFile : string
        Path to the phase2 root file.

    wf_maxLength : float
        Cut all wfs longer at wf_maxLength, if they are longer. wf_maxLength is
        the waveform length calculated as wf.GetNbinsX()*wf.GetBinWidth(1)+wf.GetBinWidth(1)/2.

    wf_minLength : int
        Reject all wfs shorter than wf_minLength. The unit is similar to wf_maxLength above.

    correctCirculatingBackground : bool
        If true, calculates the average bin content of the first 100us and subtracts this value from each bin.
        Negative entries are set to zero. Default: False
        

    Returns
    -------
    datas : array
        array of dictonaries

    '''

    tfile = TFile(pathToRootFile)
    if not tfile:
        raise IOError("TFile is empty or does not exist: " + tfile)

    # fix for phase3 data
    location = None
    if bg.phase3:
        if pathToRootFile.find('fwd') > -1:
            location = 'FWD_'
        elif pathToRootFile.find('bwd') > -1:
            location = 'BWD_'
        else: location = ''

    datas = []
    # extract the runnumber and subrunnumber from the filename
    numbers = re.findall('\d+',pathToRootFile)
    runNumber = numbers[-2]
    subRunNumber = numbers[-1]

    # predefine all variables
    tsV = 0
    lerCV = 0
    herCV = 0
    lerGV = 0
    herGV = 0
    lerIDV = 0
    herIDV = 0

    for pico in bg.picos:
        # fix for phase3 data
        if bg.phase3:
            pico = location + pico

        # read in pico specific data
        # check each time if branch exists
        ts = tfile.Get(pico + '_' + bg.items_pico['ts'])
        if not ts:
            tsV = getTimestamp(makePathToInfo(pathToRootFile))
        else:
            tsV = ts.GetBinContent(1)

        lerC = tfile.Get(pico + '_' + bg.items_pico['lerCur'])
        herC = tfile.Get(pico + '_' + bg.items_pico['herCur'])
        if not (lerC or herC):
            lerCV = math.nan
            herCV = math.nan
        else:
            lerCV = lerC.GetBinContent(1)
            herCV = herC.GetBinContent(1)

        lerG = tfile.Get(pico + '_' + bg.items_pico['lerBGate'])
        herG = tfile.Get(pico + '_' + bg.items_pico['herBGate'])
        if not (lerG or herG):
            lerGV = math.nan
            herGV = math.nan
        else:
            lerGV = lerG.GetBinContent(1)
            herGV = herG.GetBinContent(1)

        lerID = tfile.Get(pico + '_' + bg.items_pico['lerID'])
        herID = tfile.Get(pico + '_' + bg.items_pico['herID'])
        if not (lerID or herID):
            lerIDV = math.nan
            herIDV = math.nan
        else:
            lerIDV = lerID.GetBinContent(1)
            herIDV = herID.GetBinContent(1)

        for channel in bg.channels:
            data = {}
            tth1 = tfile.Get(pico
                            + '_'
                            + channel
                            + '_'
                            + bg.items_channel['reco'])

            if not isinstance(tth1, TH1):
                raise Exception(
                        "{0}: Reco {1} {2} not available".format(
                            pathToRootFile, pico,channel)
                            )

            # in some cases it makes sense to cut all wfs at a certain length
            if isinstance(wf_maxLength, lnumbers.Number) and (tth1.GetNbinsX()*tth1.GetBinWidth(1)+tth1.GetBinWidth(1)/2.) > wf_maxLength:
                xlow = tth1.GetBinCenter(1)-tth1.GetBinWidth(1)/2.
                maxBins = int(round(wf_maxLength/tth1.GetBinWidth(1)))
                xhigh = tth1.GetBinCenter(maxBins)+tth1.GetBinWidth(1)/2.

                #xhigh = tth1.GetBinCenter(wf_maxLength)+tth1.GetBinWidth(wf_maxLength)/2
                _tth1 = TH1F(f'{tth1.GetName()}_shrinked', f'{tth1.GetName()}_shrinked', maxBins, xlow, xhigh)

                for _i in range(tth1.GetNbinsX()):
                    _tth1.SetBinContent(_i+1, tth1.GetBinContent(_i+1))
                tth1 = _tth1
                
            # ....or to reject too small wfs
            if isinstance(wf_minLength, lnumbers.Number) and (tth1.GetNbinsX()*tth1.GetBinWidth(1)+tth1.GetBinWidth(1)/2.) < wf_minLength:
                continue


            _binNo = int(round(100/tth1.GetBinWidth(1)))
            meanContent100us = tth1.Integral(1,_binNo)/_binNo
            data['circB_mean'] = meanContent100us
            if correctCirculatingBackground:
                # Suppress normal background to highlight injection background.
                # Find average bin content of the first 100 us and subtract this value from each bin. 
                for _i in range(tth1.GetNbinsX()):
                    tth1.SetBinContent(_i+1, tth1.GetBinContent(_i+1)-meanContent100us if tth1.GetBinContent(_i+1)-meanContent100us > 0 else 0)


            onePeV = tfile.Get(pico
                                + '_'
                                + channel
                                + '_'
                                + bg.items_channel['1pe']
                                ).GetBinContent(1)
            onePeRatioV = tfile.Get(pico
                                    + '_'
                                    + channel
                                    + '_'
                                    + bg.items_channel['1peVStotal']
                                    ).GetBinContent(1)

            data['maxVal'] = tth1.GetMaximum()
            data['minVal'] = tth1.GetMinimum()

            # get the time until 50/90% of the MIPs can be found in the wf
            timeTo33 = 0
            timeTo50 = 0
            timeTo66 = 0
            timeTo90 = 0
            sum = 0
            totalMIPs = tfile.Get( pico
                            + '_'
                            + channel
                            + '_'
                            + bg.items_channel['mips']
                            ).GetBinContent(1)
            if totalMIPs > 0:
                for i in range(tth1.GetNbinsX()):
                    sum += tth1.GetBinContent(i+1)
                    if (timeTo33 == 0) and (sum / totalMIPs > 0.33):
                        timeTo33 = tth1.GetBinCenter(i+1)
                    elif (timeTo50 == 0) and (sum / totalMIPs > 0.5):
                        timeTo50 = tth1.GetBinCenter(i+1)
                    elif (timeTo66 == 0) and (sum / totalMIPs > 0.66):
                        timeTo66 = tth1.GetBinCenter(i+1)
                    elif (timeTo90 == 0) and (sum / totalMIPs > 0.9):
                        timeTo90 = tth1.GetBinCenter(i+1)
                        break
            data['timeTo33'] = timeTo33
            data['timeTo50'] = timeTo50
            data['timeTo66'] = timeTo66
            data['timeTo90'] = timeTo90


            data['maxValBin'] = tth1.GetMaximumBin()
            data['rate'] = tfile.Get( pico
                            + '_'
                            + channel
                            + '_'
                            + bg.items_channel['mips']
                            ).GetBinContent(1) \
                            / ( (tth1.GetBinCenter( tth1.GetNbinsX() )
                                +tth1.GetBinWidth(1)/2.0)
                                *1.e-6) \
                                / 1.e3\
                                / 1
            # for the three lines above, top to bottom
                # time per bin, phase2/phase3 = 1e-6/10e-6
                # normalization to kHz
                # area of scintillator
            data['timestamp'] = tsV
            data['time'] = unixToHumanJSTTime(tsV)
            data['runNumber'] = runNumber
            data['subRunNumber'] = subRunNumber
            data['location'] = pico
            data['channel'] = channel
            data['lerCurr'] = lerCV
            data['herCurr'] = herCV
            data['lerGate'] = lerGV
            data['herGate'] = herGV
            data['lerID'] = lerIDV
            data['herID'] = herIDV
            data['NbinsX'] = tth1.GetNbinsX()
            data['xBinWidth'] = tth1.GetBinWidth(1)
            data['wf_len'] = tth1.GetBinWidth(1)*data['NbinsX']+tth1.GetBinWidth(1)/2.
            data['path'] = pathToRootFile
            data['OnePe'] = onePeV
            data['OnePeRatio'] = onePeRatioV
            data['rawExists'] = os.path.isfile(makePathToRaw(pathToRootFile))


            datas.append(data)

    # sum the rate over all channels
    sum_rate = 0
    for tmp_data in datas:
        # exclude outliers. This value is arbitrarily chosen.
        if tmp_data['rate'] > 1e3:
            continue
        sum_rate += tmp_data['rate']

    if sum_rate > 1e6 or sum_rate < 1:
        sum_rate = np.NaN

    for tmp_data in datas:
        tmp_data['rate_sum'] = sum_rate

    tfile.Close()
    return datas



def getPandasDF(files=None, events=None, onePerRunNumber=False, dtype = 'online', hideProgress=False, **kwargs):
    '''
    Creates a pandas.DataFrame from the data given in the range defined by
    first- and lastEvent.

    Parameters
    ----------
    files : list of strings
        list of files paths of the files to be analyzed

        If files is not give, events must be given!

    events : tuple of two ints, e.g. (int, int)
        10 digit number for each the first and last event which should be considered.
        The first 6 digits are the runNumber, the last four the subRunNumber,
        e.g. (5001080000, 5001080199)

        If no events are given, a list fo files must be given!

    onePerRunNumber : bool
        Default: False
        If set to True, only one file per runNumber is considered. This can be
        useful when working with eg. 1pe integral data since it is the same
        value for one run.

    dtype : str
        Options are 'physics' or 'online', depending on the data you want.

    hideProgess : bool
        Turn on/off the progress bar by tqdm. Default: False

    **kwargs
        ...are given to the readEvent method. See cat.readEvent for further details.
    '''

    # predefine variables
    tmin = 0
    tmax = 0
    data = []
    datas = []
    counter = 0
    lrunNumber = []

    if not (isinstance(files, list) or isinstance(files, np.ndarray)):
        if len(events) != 2:
            raise Exception('events has to few/many entries!')
        files = _findByIndexFile(events[0],
                                events[1],
                                dtype = dtype
                                )

    for file in tqdm(files, disable=hideProgress, desc='Files'):
        # if 'oneperrunnumber' is true, check if run was processed already
        cRunNumber = int(re.findall('\d+',file)[-2])
        if onePerRunNumber:
            if cRunNumber in lrunNumber:
                continue
        # keep a list of which runnumbers are processed for
        # 1. if onePerRunNumber is true
        # 2. to print out how many files were found
        lrunNumber.append(cRunNumber)

        try:
            # read event and find correct timerange the data was recorded in
            data = readEvent(file, **kwargs)
            if counter == 0:
                tmin = data[0]['timestamp']
                tmax = data[0]['timestamp']
            else:
                for channel in data:
                    dt = channel['timestamp']
                    if tmin > dt:
                        tmin = dt
                    elif tmax < dt:
                        tmax = dt
            datas.append(data)
        except Exception as exce:
            pass
#            print('{0} not available: {1}'.format(file, exce.args))

        counter += 1

    print('Found {} files.'.format(len(lrunNumber)))
    print('TMin: {0} \t TMax: {1}'.format(tmin,tmax))
    print('TMin: {0} JST\t TMax: {1} JST'.format(unixToHumanJSTTime(tmin), unixToHumanJSTTime(tmax)))

    try:
        # create pandas.DataFrame and define some types
        output = pd.DataFrame([pair for data in datas for pair in data])
        output['time'] = pd.to_datetime(output['time'],format="%Y-%m-%d %H:%M:%S.%f")
        output['runNumber'] = output.runNumber.astype('int')
        output['subRunNumber'] = output.subRunNumber.astype('int')
    except KeyError as keyerr:
        print('Key error was caught. Key {0} not available!'.format(keyerr.args))
        return None
    return output


def _findByIndexFile(   firstEvent,
                        lastEvent,
                        pathToIndexFile=None,
                        pathToData=None,
                        dType="online"):
    '''
    Finds the events in the given range thanks to a index file.
    This is much faster when in the MPP since the data is on a nas and network is slow.

    Parameters
    ----------
    firstEvent : int
    lastEvent : int

    pathToIndexFile : str
        Absolute path to the index file.

    pathToData : str
        The index file only contains the path information in the data folder.
        Therefore the path to the data folder needs to be specified.

    dType : str
        This can be either 'physics' or 'online'.

    Returns
    -------
    list of str
        Absolute paths to the files.
    '''
    if pathToIndexFile == None:
        pathToIndexFile = getIndexFilePath()
    if pathToData == None:
        pathToData = getDataPath()

    # extract the run and subrun numbers
    fRunN = str(firstEvent)[:6]
    fSubRunN = str(firstEvent)[-4:]
    lRunN = str(lastEvent)[:6]
    lSubRunN = str(lastEvent)[-4:]

    # make a list with all runnumbers in the corresponding range
    runList = [run for run in range(int(fRunN),int(lRunN)+1)]
    subRunList = [run for run in range(int(fSubRunN),int(lSubRunN)+1)]

    with open(pathToIndexFile, 'r') as f:
        lines = f.readlines()
        lines = [line.rstrip('\n') for line in lines if line.find('.root') > -1]
        lines = [os.path.join(pathToData,line) for line in lines if __searchByIndexFile_helper__(dType,line,fRunN,fSubRunN,lRunN,lSubRunN)]
        print('Found {0} possible files.'.format(len(lines)))
        return lines


def __searchByIndexFile_helper__(   dType,
                                    x,
                                    fRunN,
                                    fSubRunN,
                                    lRunN,
                                    lSubRunN
                                    ):
    if not x.find(dType) > -1:
        return False
    if not (int(fRunN) <= int(re.findall('\d+',x)[-2]) <= int(lRunN)):
        return False
    if (int(fRunN) == int(re.findall('\d+',x)[-2]) and int(fSubRunN) > int(re.findall('\d+',x)[-1])):
        return False
    if ( int(lRunN) == int(re.findall('\d+',x)[-2]) and int(lSubRunN) < int(re.findall('\d+',x)[-1]) ):
        return False
    return True
