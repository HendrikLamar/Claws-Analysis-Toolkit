from . import _global as bg

import math
import configparser, os, re, copy, platform

import pandas as pd
import pytz as tz
from datetime import datetime, timezone
import matplotlib
import matplotlib.pyplot as plt

from ROOT import TFile, TH1, TH1I, TH1F, TCanvas


def getMachineDependendPathToData():
    ''' Returns the root path to the data depending on which machine the code
    is executed.
    Implemented machines: pcbelle20, sapporo, claws2-daq
    '''
    machine = platform.node()

    if machine == 'sapporo':
        return "/home/hwindel/NAS_futDet/claws/data/phase2/raw"
    elif machine == 'pcbelle20' or 'pcilc11':
        return "/home/iwsatlas1/hwindel/NAS_futDet/claws/data/phase2/raw"
    elif machine == 'claws2-daq':
        return "/mnt/claws_data/data"


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


def findByIndexFile(firstEvent, 
                    lastEvent,                     
                    fidx,
                    rpath,
                    dType="online"
                    ):
    '''
    Finds the files you need by searching in the index file. This works only
    when using the setup at MPP. There is no index file at the claws computer
    at KEK.
    '''

    # extract the run and subrun numbers
    fRunN = str(firstEvent)[:6]
    fSubRunN = str(firstEvent)[-4:]
    lRunN = str(lastEvent)[:6]
    lSubRunN = str(lastEvent)[-4:]

    # make a list with all runnumbers in the corresponding range
    runList = [run for run in range(int(fRunN),int(lRunN)+1)]
    print('runList: ', runList)
    subRunList = [run for run in range(int(fSubRunN),int(lSubRunN)+1)]
    print( 'subRunList: ', subRunList)

#    fidx = '/home/iwsatlas1/hwindel/workspace/beast_daq/phase2/claws/analysis/jupyter/index.idx'
#    rpath = getMachineDependendPathToData()
    print( 'rpath: ', rpath)

    with open(fidx, 'r') as f:
        lines = f.readlines()
        print( 'lines: ', lines)
        lines = [line.rstrip('\n') for line in lines if line.find('.root') > -1]
        print( 'lines: ', lines)
        lines = [os.path.join(rpath,line) for line in lines if __searchByIndexFile_helper__(dType,line,fRunN,fSubRunN,lRunN,lSubRunN)]
        print('Found {0} possible files.'.format(len(lines)))
        return lines


def unixToHumanJSTTime(ts):
    '''
    Transforms the unix timestamp taken at KEK (JST) to human readable time
    format in JST timezone.

    If this is used at KEK things might go wrong. Please verify!
    '''
    utc_ts = datetime.fromtimestamp(ts)
    jst = tz.timezone('Asia/Tokyo')
    cet = tz.timezone('Europe/Berlin')
    jst_time = utc_ts + jst.utcoffset(utc_ts) - cet.utcoffset(utc_ts)
    return jst_time.strftime("%Y-%m-%d %H:%M:%S.%f")


def readEvent( pathToRootFile ):
    '''
    Reads one event/file.

    Parameters
    ----------
    pathToRootFile : string
        Path to the phase2 root file.

    '''

    tfile = TFile(pathToRootFile)
    if not tfile:
        raise Exception("TFile is empty or does not exist: " + tfile)

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
            tth1 = tfile.Get(pico + '_' + channel + '_' + bg.items_channel['reco'])
            if not tth1:
                raise Exception("{0}: Reco {1} {2} not available".format(pathToRootFile,pico,channel))


            onePeV = tfile.Get(pico + '_' + channel + '_' + bg.items_channel['1pe']).GetBinContent(1)
            onePeRatioV = tfile.Get(pico + '_' + channel + '_' + bg.items_channel['1peVStotal']).GetBinContent(1)
            

            data['maxVal'] = tth1.GetMaximum()
            data['minVal'] = tth1.GetMinimum()
            data['maxValBin'] = tth1.GetMaximumBin()
            data['rate'] = tfile.Get(pico + '_' + channel + '_' + bg.items_channel['mips']).GetBinContent(1) / ((tth1.GetBinCenter(tth1.GetNbinsX())+tth1.GetBinWidth(1)/2.)*1.e-6)/1.e6/4.
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
            data['xbins'] = tth1.GetNbinsX()
            data['path'] = pathToRootFile
            data['1pe'] = onePeV
            data['1peRatio'] = onePeRatioV

            datas.append(data)

    tfile.Close()
    return datas


def getPandasDF(firstEvent, 
                lastEvent,
                rpath, 
                onePerRunNumber=False,
                pathToIndex='/home/iwsatlas1/hwindel/workspace/beast_daq/phase2/claws/analysis/jupyter/index.idx'
                ):
    '''
    Creates a pandas.DataFrame from the data given in the range defined by
    first- and lastEvent.

    Parameters
    ----------
    firstEvent : int
        10 digit number with the first event which should be considered.
        The first 6 digits are the runNumber, the last four the subRunNumber,
        e.g. 5001080000

    lastEvent : int
        10 digit number with the last event which should be considered.
        The first 6 digits are the runNumber, the last four the subRunNumber,
        e.g. 5001080199

    onePerRunNumber : bool
        Default: False
        If set to True, only one file per runNumber is considered. This can be
        useful when working with eg. 1pe integral data since it is the same
        value for one run.
    '''

    # predefine variables
    tmin = 0
    tmax = 0
    data = []
    datas = []
    counter = 0
    lrunNumber = []
    files = findByIndexFile(firstEvent,
                            lastEvent, 
                            pathToIndex,
                            rpath
                            )
#    files = findRootFiles_online(firstEvent, lastEvent)

    for file in files:
        # if 'oneperrunnumber' is true, check if run was processed already
        cRunNumber = int(re.findall('\d+',file)[-2])
        if onePerRunNumber:
            if cRunNumber in lrunNumber:
                continue
        # keep a list of which runnumbers are processed for
        # 1. if onePerRunNumber is true
        # 2. to print out how many files were found
        lrunNumber.append(cRunNumber)

        if counter %10000 == 0:
            print('#{0} {1}'.format(counter,file))
        
        try:
            # read event and find correct timerange the data was recorded in
            data = readEvent(file)
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
            print('{0} not available: {1}'.format(file, exce.args))
        
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

def th1XScale( th1, factor ):
    '''
    Rescales the th1 xAxis by the given factor. Remember to change the xAxis legend!

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
    xup_old = th1.GetXaxis().GetBinUpEdge(th1.GetNbinsX())+th1.GetBinWidth(th1.GetNbinsX())
    xup_new = xup_old * factor
    th1.SetBins(nbins, xlow, xup_new)

    return th1


def makeAvg( dfPaths, picos, channels):
    '''
    Averages the waveforms over the given events. Per event all given locations and channels are summed up.

    Parameters
    ----------
    dfPaths : pandas.DataFrame
        Data frame with index and paths for events only.

    picos : list of str
        List with all locations/picos which should be considered. Possible locations:
            - 'Top_Forward'
            - 'Top_Backward'
            - 'Bottom_Forward'
            - 'Bottom_Backward'
            e.g. picos = ['Top_Forward', 'Bottom_Backward']

    channels : list of str
        List with all channels which should be considered. Possible channels:
            - 'A'
            - 'B'
            - 'C'
            - 'D'
            e.g. channels = ['A', 'C']

    Returns
    -------
    Tuple of ROOT.TH1I and ROOT.TH1I
    The first one is the average waveform. The second is the signal heigth histogram.
    '''

    data = dfPaths.drop_duplicates()

    highestVal = -1
    th1 = 0
    th2 = 0
    thh = TH1F('Signal heigth histogram','Signal height histogram', 1000,0,1000)
    isFirst = True
    counter = 0
    for i, c1 in data.iterrows():
        firstEvent = True
        tfile = TFile(c1[0])
        #print('tfile: {0} \t {1}'.format(tfile, c1[0]))
        if not tfile: raise Exception('File not found: {}'.format(c1[0]))

        for pico in picos:
            for channel in channels:
                dth1 = tfile.Get(pico + '_' + channel + '_' + bg.items_channel['reco'])
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
    print('# of overlaid events: {0}\nmaxVal: {3}\nPicos {1}\nChannels {2}'.format(counter, picos, channels, highestVal))
    th1.Scale(1./counter)


    return th1, thh



def makeFFT(th1Wf):
    '''
    Determines the distance between all peaks and weigths that with the product of
    the heights of the two peaks. Therefore, the output unit on the y axis is 
    MIPs x MIPs / bin width.

    Parameters
    ----------
    th1Wf : ROOT.TH1
        Waveform to analyse.

    Returns
    -------
    ROOT.TH1F
        Return the 'low cost' fast fourier transformation of the input waveform.
    '''

    ddata = [0. for bin in range(th1Wf.GetNbinsX())] # create dummy container filled with zeros

    # loop through th1Wf. Measure distance between two entries and fill that into 
    # the data container weigthed with the product of the two entries
    for bin in range(1,th1Wf.GetNbinsX()-1):
        y0 = th1Wf.GetBinContent(bin)
        for sbin in range(bin+1,th1Wf.GetNbinsX()):
            y1 = th1Wf.GetBinContent(sbin)
            ddata[sbin-bin] += y0*y1
    
    # transform the data from the dummy container to a th1Wf
    th2 = TH1F('tff','tff',th1Wf.GetNbinsX(),0,th1Wf.GetBinCenter(th1Wf.GetNbinsX())) # create new th1Wff for fft content
    th2.SetXTitle('Time [ms]')
    th2.SetYTitle('MIPs x MIPs / {0}us [x1000]'.format('8'))
    th2.SetTitle('{0}_{1}'.format(th1Wf.GetTitle(),'FFT'))
    for bin in range(len(ddata)):
        th2.SetBinContent(bin,ddata[bin])
    th2.Scale(1./1000.)
    
    return th2



def saveCanvas( c1, name, info = ""):
    ''' Saves the TCanvas in several formats + info file.
    2018-06-01 jpg does not work because of missing codecs.
    '''
    c1.SaveAs(name + ".root")
    c1.SaveAs(name + ".pdf")
    c1.SaveAs(name + ".C")
    c1.SaveAs(name + ".png")
#    c1.SaveAs(name + ".jpg")
    info = name + "\n" + info
    with open(name+".info",'w') as f:
        f.write(info)

def saveFig( fig, ax, title=""):
    dtitle = 0
    if title == "":
        dtitle = ax.get_title()
    else:
        dtitle = title
    fig.savefig('{0}.pdf'.format(dtitle), dpi=150)
    fig.savefig('{0}.jpg'.format(dtitle), dpi=150)
    fig.savefig('{0}.png'.format(dtitle), dpi=150)


def makePathToInfo( filepath ):
    '''
    Returns the info*.ini file to the given .root file.
    '''
    filepath = filepath.replace('online', 'raw/physics',1)
    filepath = filepath.replace('online', 'info')
    filepath = filepath.replace('root','ini')

    return filepath


def getTimestamp( pIni ):
    '''
    Returns the timestamp written in the info*.ini file.
    This method is needed since the timestamp is not available in some root files
    at the beginning of phase2.
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


def TH1ToDf(th1):
    '''
    Converts a ROOT.TH1 to a pandas.DataFrame.

    Parameters
    ----------
    th1 : ROOT.TH1
        Root histogram.

    Returns
    -------
    pandas.DataFrame
        The index is the x axis, bin content is column 'y'.
    '''
    d = { 'x' : [], 'y' : []}
    for i in range(1,th1.GetNbinsX()):
        d['x'].append(th1.GetBinCenter(i))
        d['y'].append(th1.GetBinContent(i))

    return pd.DataFrame(d).set_index('x')




# This class was never used...if you use it, use it with care!
class Peast:

    def __init__(self, firstEvent=0, lastEvent=0, onePerRunNumber=False, df=None):
        if firstEvent == 0:
            self.df = df
        else:
            self.df = getPandasDF( firstEvent, lastEvent, onePerRunNumber)     

        self.fft = 0
        self.avg = 0
        self.peak = 0
        self.__makeWhat()
    
    def makeAvg(self, picos=['Top_Forward', 'Top_Backward', 'Bottom_Forward', 'Bottom_Backward'], channels=['A', 'B' ,'C', 'D']):
        '''
        Averages the waveforms over the given events. Per event all given locations and channels are summed up.

        Parameters
        ----------
        picos : list of str
            List with all locations/picos which should be considered. Possible locations:
                - 'Top_Forward'
                - 'Top_Backward'
                - 'Bottom_Forward'
                - 'Bottom_Backward'
                e.g. picos = ['Top_Forward', 'Bottom_Backward']

        channels : list of str
            List with all channels which should be considered. Possible channels:
                - 'A'
                - 'B'
                - 'C'
                - 'D'
                e.g. channels = ['A', 'C']
        '''

        self.avg, self.peak = makeAvg(self.df.path.to_frame(), picos, channels)
        self.__makeWhat()

    def makeFft(self):
        self.fft = makeFFT(self.avg)
        self.__makeWhat()

    def plot(self, what):
        '''
        Plots the graph determined in the 'what' parameter.
        
        Parameter
        ---------
        what : string
            Options are:
                - fft
                - peak
                - avg
                - overview
        
        '''
        dth1 = 0 
        dkey = 0
        for key, value in self.what.items():
            if key == what:
                dkey = key
                dth1 = value
       
        if dkey == 'overview':
            self.plotOverview(self.df)
        else:
            self.c1 = TCanvas('sc1','sc1',1066,600)
            dth1.Draw('hist')

            self.c1.Draw()


    def save(self, fname, info=''):
        '''
        Saves the last plotted ROOT graph with the given file name. The info parameter
        can be filled additional information and is stored in a text file.

        Parameters
        ----------
        fname : string
            Filename only. The ending (e.g. .pdf or .png) is set automatically.
        info : string
            Info text for more specific explanations for the plot. The text is saved
            in a text file.
        '''

#    def df(self):
#        '''
#        Returns
#        -------
#        pandas.DataFrame
#            Returns the data.
#        '''
#        return self.df

    @staticmethod
    def plotOverview(df = None):
        try:
            df.plot(x='timestamp', y=['herCurr','lerCurr','rate'], secondary_y=['rate'], figsize=(18,9), ls='None', marker='.') 
        except (NameError, AttributeError, KeyError) as excep:
            if type(excep).__name__ == 'KeyError':
                print('KeyError\nOne of the following keys does not exist: timestamp, herCurr, lerCurr, rate')
            else:
                print('The parameter you entered is no pandas.DataFrame')

    def __makeWhat(self):
        self.what = {'fft' : self.fft, 'avg' : self.avg, 'peak' : self.peak, 'overview' : True}

    def xScale(self, factor, what):
        th1XScale(self.what[what], factor)































