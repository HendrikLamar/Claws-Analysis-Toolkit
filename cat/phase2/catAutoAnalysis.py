import matplotlib as mpl
import matplotlib.style as style
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import time as t
import datetime
import os
import pylatex

from cat.phase2 import phase2
from cat.phase2 import _global as globalvar

from ROOT import TFile, TH1, TCanvas, TH1F, ROOT, TF1, gStyle
from copy import deepcopy

#from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Matrix, Alignat

########################################################################                       
# plotting Current of Collider ring and the measured rate
########################################################################                       
def plotColRingCurrWithRate(clawsDataframe,
                            colliderRing,
                            sampleName,
                            minuteLocator = True,
                            hourLocator = False,
                            dayLocator = False,
                            mountLocator = False,
                            minorFormatter = False,
                            majorFormatter = True
                            ):
    '''
    Parameters                                                                                 
    ----------
    
    Returns
    -------
    '''
    title = colliderRing + ' Injection Study ' + sampleName
    fig, ax = plt.subplots(figsize=(16,9))
    ax2 = ax.twinx()
    line1, = ax.plot(clawsDataframe.time, 
                     clawsDataframe.rate, 
                     ls='None', 
                     markersize=5, 
                     marker='.', 
                     color='g'
                    )
    line2, = ax2.plot(clawsDataframe.time, 
                      clawsDataframe.herCurr, 
                      ls='None', 
                      markersize=5, 
                      marker='.', 
                      color='b'
                     )
    ax.set_ylabel('rate [MHz/cm²]')
    ax.set_yscale("log", nonposy='clip')
    ax2.set_ylabel(colliderRing + ' current [mA]')
    ax.set_xlabel('time')
    
    if minuteLocator:
        mmin = mdates.MinuteLocator(byminute=range(0,60,15))
        ax.xaxis.set_minor_locator(mmin)
    if hourLocator:
        mhours = mdates.HourLocator()
        ax.xaxis.set_major_locator(mhours)
    if dayLocator:
        mdays = mdates.DayLocator(interval=5)
        ax.xaxis.set_major_locator(dayLocator)
    if mountLocator:
        mmonth = mdates.MonthLocator(bymonthday=1)
        ax.xaxis.set_major_locator(mountLocator)
    if minorFormatter:
        miformatter = mdates.DateFormatter('%M')
        ax.xaxis.set_minor_formatter(miformatter)
    if majorFormatter:
        maformatter = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(maformatter)

    ax.xaxis.set_tick_params(pad=18, which='major')
    plt.title(title)
    leg = plt.legend((line1, line2), 
                     ('Rate', 'HER Current'), 
                     markerscale=5, loc='best'
                    )
    figTitle = title.replace(" ", "_")
    
    return fig, ax, figTitle

########################################################################                       
# plot the averaged waveform of claws with it's cumulative sum
########################################################################  
def plotAveragedWf( awfDataframe, additionalInfo=''):
    if not additionalInfo == '':
        additionalInfo = ': ' + additionalInfo
    title = 'Averaged Waveform' + additionalInfo

    ax = awfDataframe.plot(secondary_y='cumsum', 
                 drawstyle='steps-mid', 
                 linewidth=0.5, 
                 legend=None, 
                 title=title, 
                 ylim=(0,260))#, figsize=(6.8,4.8))
    ax.legend(['Injection Wf', 'second'],loc='lower right')
    ax.set_ylabel('MIPs / 8 $\\mu$s')
    #ax.set_ybound(0,320)
    ax.right_ax.set_ylabel('Cumulative [%]')
    ax.set_xlabel('Time [ms]')
    ax.set_xticks(np.arange(0,10,0.5), minor=True)
    ax.get_figure().text(0.75,0.85,
                         'Claws preliminary', 
                         color='gray', 
                         alpha=0.5, 
                         fontsize=15)
    ax.get_figure().text(0.48,0.85, 
                         '{} kHz/cm²'.format(int(round(
                             awfDataframe.y.sum()
                                 /5e-3
                                 /4.
                                 /1e3,
                             0
                         ))), 
                         fontsize=15)
    fig = ax.get_figure()
    
    return fig, ax, title

########################################################################                       
# plot a histogram how often the peaks distance, identified in the calws 
# waveform, occurs (slang: simple Fast Fourier Transform (FFT) )
########################################################################  
def plotPeakAna( paDataframe, additionalInfo=''):
    if not additionalInfo == '':
        additionalInfo = ': ' + additionalInfo
    title = 'Peak Analysis' + additionalInfo
    ax = paDataframe.plot(drawstyle='steps-mid', 
                          linewidth=0.5, 
                          legend=None, 
                          title=title, 
                          logy=False, 
                          logx=True, 
                          xlim=(1,1e4), 
                          ylim=(1e-1,70)
                         )
    ax.legend(['Injection Wf', 'second'],loc='upper left')
    ax.set_ylabel('Weighted Entries [MIP² / (8 $\\mu$s)³]')
    #ax.set_ybound(0,320)
    ax.set_xlabel('Peak to Peak Distance [ms]')
    #ax.set_xticks(np.arange(0,5,0.5), minor=True)
    ax.get_figure().text(0.75,0.85,
                         'Claws preliminary', 
                         color='gray', 
                         alpha=0.5, 
                         fontsize=15
                        )
    fig = ax.get_figure()
    
    return fig, ax, title

########################################################################                       
# Fast Fourier Transform (FFT) waveform 
########################################################################  
def autoAnalysis(startDate,
                 endDate,
                 pathToData,
                 dataType = 'online',
                 defRundNum = ('','')
                ):
    '''
    Parameters                                                                                 
    ----------
    day : str
        a date given in the format JJJJ-MM-DD
        
    pathToData : str
        absolute path, can be optained by:
        os.path.abspath( relative_path )
        
    dataType : str
        flog wich type of recorded data should be used
    
    Returns
    -------
    figures
        writes out a new file for each figure created
        
    figureFileNames : dictonary
        contains the filename ( string ), of each of the craeted 
        figures
        
    variablesOfInterest : dictonary
        contains all variables of interest 
    '''
    #---------------------------------------------------------------
    # Check if dates are set, prepare time slot given to plotting 
    # functions, and optain all filepaths for the days data
    #---------------------------------------------------------------
    if startDate == endDate:
        timeSlot = startDate
    else:
        timeSlot = startDate + ' until ' + endDate
    dates = phase2.provideListOfDates( startDate,
                                     endDate,
                                     pathToData
                                     )
    
    files = phase2.findFiles( dates, 
                            pathToData,
                            defRundNumRange=defRundNum,
                            dType=dataType 
                            )
    #---------------------------------------------------------------
    # create a pandas dataframe from the raw data
    #---------------------------------------------------------------
    dataframe = phase2.getPandasDF( files ) 
   
    #---------------------------------------------------------------
    # Plot the rate and the Current for the corresponding collider ring
    #---------------------------------------------------------------
    fig_raHcurr, ax_raHcurr, figTitle_raHcurr = plotColRingCurrWithRate( 
                                                    dataframe,
                                                    'HER',
                                                    timeSlot
                                                    )
    phase2.saveFig(fig_raHcurr, 
                   ax_raHcurr,  
                   title=figTitle_raHcurr,
                   jpg=False,
                   pdf=False
                  )
    fig_raLcurr, ax_raLcurr, figTitle_raLcurr = plotColRingCurrWithRate( 
                                                    dataframe,
                                                    'LER',
                                                    timeSlot
                                                    )
    phase2.saveFig(fig_raLcurr, 
                   ax_raLcurr,  
                   title=figTitle_raLcurr,
                   jpg=False,
                   pdf=False
                   )
    #---------------------------------------------------------------
    # Calculate the Average Waveform and plot it
    #---------------------------------------------------------------
    tref, sref = phase2.makeAvg( dataframe.path.to_frame(), 
                                globalvar.picos, 
                                globalvar.channels
                               )
    DataframeAvg = phase2.TH1ToDf(tref)
    averagedWaveform = makeAveragedWf( DataframeAvg )
    ## Plot
    fig_avWf, ax_avWf, figTitle_avwf = plotAveragedWf(
                                            averagedWaveform, 
                                            additionalInfo='HER'+startDate
                                            )
    print( figTitle_avwf )
    phase2.saveFig(fig_avWf, 
                   ax_avWf,  
                   title=figTitle_avwf,
                   jpg=False,
                   pdf=False
                   )
    #---------------------------------------------------------------
    # create pandas dataframe of peak distance analysied data and 
    # plot in histogram
    #---------------------------------------------------------------
    dataframe_sFFT = phase2.TH1ToDf( phase2.makeFFT( tref ) )
    ## Plot
    fig_PAna, ax_PAna, figTitle_PAna = plotPeakAna( dataframe_sFFT )
    phase2.saveFig(fig_PAna, 
                   ax_PAna,  
                   title=figTitle_PAna,
                   jpg=False,
                   pdf=False
                   )
    
    #---------------------------------------------------------------
    # Fourier Transform Averaged Waveform and plot in histogram
    #---------------------------------------------------------------
    numPyarrAvg = DataframeAvg.y.as_matrix()
    FFTnumPyarrAvg = fft(numPyarrAvg)
    binNumberForFFTfunction = len(DataframeAvg)/2
    ## Plot
    fig_FFT, ax_FFT, figTitle_FFT = plotFFT( FFTnumPyarrAvg, 
                                            binNumberForFFTfunction, 
                                            additionalInfo=''
                                            )
    phase2.saveFig(fig_FFT, 
               ax_FFT,  
               title=figTitle_FFT,
               jpg=False,
               pdf=False
               )
    #---------------------------------------------------------------
    # give the number of claws-runs taken
    # BEWARE: these runs do not correspond to the Belle II run 
    # numbers
    #---------------------------------------------------------------
    runs = []
    for i in dates:
        runs += os.listdir(pathToData + '/' + i)
        
    #---------------------------------------------------------------
    # write out variables of interest
    #---------------------------------------------------------------   
    variablesOfInterest = {'firstrun' : runs[0],
                            'lastrun' : runs[-1],
                            'meanHERcurrent' : dataframe.herCurr.mean(),
                            'maxHERcurrent' : dataframe.herCurr.max(),
                            'meanLERcurrent' : dataframe.lerCurr.mean(),
                            'maxLERcurrent' : dataframe.lerCurr.max(),
                            'meanRate' : dataframe.rate.mean(),
                            'maxRate' : dataframe.rate.max(),
                            'meanMxVal' : dataframe.maxVal.mean(),
                            'meanMaxValBin' : dataframe.maxValBin.mean(),
                            'maxMaxValBin' : dataframe.maxValBin.max(),
                            'absMaxVal' : dataframe.maxVal.max(),
                            'meanMinVal' : dataframe.minVal.mean(),
                            'minMinVal' : dataframe.minVal.min()
                            }
    
    figureFileNames = { 'raHcurr' : figTitle_raHcurr,
                         'raLcurr' : figTitle_raLcurr,
                         'averageWf' : figTitle_avwf,
                         'peakDist' : figTitle_PAna,
                         'fft' : figTitle_FFT
                        }
                 
                 
    return figureFileNames, variablesOfInterest


########################################################################                       
# Automated Analysis
########################################################################  
def autoAnalysis(dates,
                pathToData,
                dataType = 'online'
                ):
    '''
    Parameters                                                                                 
    ----------
    day : str
        a date given in the format JJJJ-MM-DD
        
    pathToData : str
        absolute path, can be optained by:
        os.path.abspath( relative_path )
        
    dataType : str
        flog wich type of recorded data should be used
    
    Returns
    -------
    figures
        writes out a new file for each figure created
        
    figureFileNames : dictonary
        contains the filename ( string ), of each of the craeted 
        figures
        
    variablesOfInterest : dictonary
        contains all variables of interest 
    '''
    #---------------------------------------------------------------
    # Check if dates are set, prepare time slot given to plotting 
    # functions, and optain all filepaths for the days data
    ################################################################
    if len(dates) > 1:
        timeSlot = dates[0] + 'until' + dates[-1]
    elif len(dates)<1:
        print( 'Please enter dates to be processed!' )
    else:
        timeSlot = dates[0]
    
    files = phase2.findFiles( dates, 
                            pathToData,
                            defRundNumRange=defRundNum,
                            dType=dataType 
                            )
    ################################################################
    # create a pandas dataframe from the raw data
    ################################################################
    dataframe = phase2.getPandasDF( files ) 
   
    ################################################################
    # Plot the rate and the Current for the corresponding collider ring
    ################################################################
    fig_raHcurr, ax_raHcurr, figTitle_raHcurr = plotColRingCurrWithRate( 
                                                    dataframe,
                                                    'HER',
                                                    timeSlot
                                                    )
    phase2.saveFig(fig_raHcurr, 
                   ax_raHcurr,  
                   title=figTitle_raHcurr,
                   jpg=False,
                   pdf=False
                  )
    fig_raLcurr, ax_raLcurr, figTitle_raLcurr = plotColRingCurrWithRate( 
                                                    dataframe,
                                                    'LER',
                                                    timeSlot
                                                    )
    phase2.saveFig(fig_raLcurr, 
                   ax_raLcurr,  
                   title=figTitle_raLcurr,
                   jpg=False,
                   pdf=False
                   )
    ################################################################
    # Calculate the Average Waveform and plot it
    ################################################################
    tref, sref = phase2.makeAvg( dataframe.path.to_frame(), 
                                globalvar.picos, 
                                globalvar.channels
                               )
    DataframeAvg = phase2.TH1ToDf(tref)
    averagedWaveform = makeAveragedWf( DataframeAvg )
    ## Plot
    fig_avWf, ax_avWf, figTitle_avwf = phase2.plotAveragedWf(
                                            averagedWaveform, 
                                            additionalInfo='HER'+startDate
                                            )
    phase2.saveFig(fig_avWf, 
                   ax_avWf,  
                   title=figTitle_avWf,
                   jpg=False,
                   pdf=False
                   )
    
    ################################################################
    # create pandas dataframe of peak distance analysied data and 
    # plot in histogram
    ################################################################
    dataframe_sFFT = phase2.TH1ToDf( phase2.makeFFT( tref ) )
    ## Plot
    fig_PAna, ax_PAna, figTitle_PAna = plotPeakAna( dataframe_sFFT )
    phase2.saveFig(fig_PAna, 
                   ax_PAna,  
                   title=figTitle_PAna,
                   jpg=False,
                   pdf=False
                   )
    
    ################################################################
    # Fourier Transform Averaged Waveform and plot in histogram
    ################################################################
    numPyarrAvg = DataframeAvg.y.as_matrix()
    FFTnumPyarrAvg = fft(numPyarrAvg)
    binNumberForFFTfunction = len(DataframeAvg)/2
    ## Plot
    fig_FFT, ax_FFT, figTitle_FFT = plotFFT( FFTnumPyarrAvg, 
                                            binNumberForFFTfunction, 
                                            additionalInfo=''
                                            )
    
    ################################################################
    # give the number of claws-runs taken
    # BEWARE: these runs do not correspond to the Belle II run 
    # numbers
    ################################################################
    runs = []
    for i in dates:
        runs += os.listdir(pathToData + '/' + i)
        
    ################################################################
    # write out variables of interest
    ################################################################   
    variablesOfInterest = {'firstrun' : runs[0],
                            'lastrun' : runs[-1],
                            'meanHERcurrent' : dataframe.herCurr.mean(),
                            'maxHERcurrent' : dataframe.herCurr.max(),
                            'meanLERcurrent' : dataframe.lerCurr.mean(),
                            'maxLERcurrent' : dataframe.lerCurr.max(),
                            'meanRate' : dataframe.rate.mean(),
                            'maxRate' : dataframe.rate.max(),
                            'meanMxVal' : dataframe.maxVal.mean(),
                            'meanMaxValBin' : dataframe.maxValBin.mean(),
                            'maxMaxValBin' : dataframe.maxValBin.max(),
                            'absMaxVal' : dataframe.maxVal.max(),
                            'meanMinVal' : dataframe.minVal.mean(),
                            'minMinVal' : dataframe.minVal.min()
                            }
    
    figureFileNames = { 'raHcurr' : figTitle_raHcurr,
                         'raLcurr' : figTitle_raLcurr,
                         'averageWf' : figTitle_avwf,
                         'peakDist' : figTitle_PAna,
                         'fft' : figTitle_FFT
                        }
                 
                 
    return figureFileNames, variablesOfInterest

########################################################################                       
# Automated TeX and PDf file creation
########################################################################  
def writeAnalysisToFile(figureFileNames, 
                        pathTofigures,
                        varOfI,
                        fileName):
    '''
    Parameters                                                                                 
    ----------
    
    Returns
    -------
    '''
    fileName = os.path.abspath('./') +'/'+ fileName
    clawslogo_filename = os.path.abspath('claws3.png')
    geometry_options = {"tmargin": "1cm"}#, "lmargin": "3cm"}
    doc = pytex.Document(geometry_options=geometry_options)

    #---------------------------------------------------------------
    # Create Header of PDF File
    #---------------------------------------------------------------
    doc.preamble.append(pytex.Command('title', 'CLAWS ++ Background Examination'))
    doc.preamble.append(pytex.Command(
        'author','Hendrik Windel and Thomas Kraetzschmar'))
    doc.preamble.append(pytex.Command('date', pytex.NoEscape(r'\today')))
    doc.append(pytex.NoEscape(r'\maketitle'))

    with doc.create(pytex.Figure(position='h!')) as clawslogo:
        clawslogo.add_image(clawslogo_filename, width='120px')

    
    #---------------------------------------------------------------
    # Start document input
    #---------------------------------------------------------------
    with doc.create( pytex.Section('Overview of Injection Study' + fileName[23:-1]) ):
        doc.append('In this time periode, the runs '
                   + varOfI['firstrun']
                   + ' up to '
                   + varOfI['lastrun']
                   + 'were performed.')

        with doc.create(pytex.Itemize()) as itemize:
            itemize.add_item('Mean HER Beam Current: ' 
                             + str( round(varOfI['meanHERcurrent'], 2) ) 
                             + ' [mA]')
            itemize.add_item('Max HER Beam Current: ' 
                             + str( round(varOfI['maxHERcurrent'], 2) ) 
                             + ' [mA]')
            itemize.add_item('Mean LER Beam Current: ' 
                             + str( round(varOfI['meanLERcurrent'], 2) ) 
                             + ' [mA]')
            itemize.add_item('Max LER Beam Current: ' 
                             + str( round(varOfI['maxLERcurrent'], 2) ) 
                             + ' [mA]')
            itemize.add_item('Mean Rate: ' 
                             + str( round(varOfI['meanRate'], 2) ) 
                             + ' [MHz/cm²]')
            itemize.add_item('Maximum Rate: ' 
                             + str( round(varOfI['maxRate'], 2) ) 
                             + ' [MHz/cm²]')
            
        with doc.create(pytex.Subsection('Variables of Interest')):
            with doc.create(pytex.Figure(position='h!')) as fig1:
                fig1.add_image( figureFileNames['raHcurr'], width='480px' )
                #fig1.add_caption('Look it\'s on its back')
            with doc.create(Figure(position='h!')) as fig1:
                fig1.add_image( figureFileNames['raLcurr'], width='480px')
                #fig1.add_caption('Look it\'s on its back')
            with doc.create(Figure(position='h!')) as fig1:
                fig1.add_image( figureFileNames['peakDist'], width='480px')
                #fig1.add_caption('Look it\'s on its back')
            with doc.create(Figure(position='h!')) as fig1:
                fig1.add_image( figureFileNames['fft'], width='480px')
                #fig1.add_caption('Look it\'s on its back')

    print( fileName )
    doc.generate_pdf( fileName, clean_tex=False)
    
########################################################################                       
# This routine starts and handles the automated analysis
########################################################################  
#####################################################---------------------------
ca = c.Calendar()
sundays = []
calender = []
processedDay = t.strftime("%Y-%m-%d")
calenderFoulded = ca.yeardatescalendar(2019)
for quaterYear in calenderFoulded:
    for month in quaterYear:
        for week in month:
            sundays.append(week[6])
            calender += week
sundays = [i.strftime('%Y-%m-%d') for i in sundays ]
calender = [i.strftime('%Y-%m-%d') for i in calender ]
pathToData = os.path.abspath('../raw')
dataType = 'online'
processedDay = ''

while True:
    t.localtime(ts)
    localDate = t.strftime("%Y-%m-%d")# %w %H:%M:%S.%f")
    listOfAllDates = getAllDatesOfData( pathToData )
    listOfAllDates.append(processedDay)
    listOfAllDates.append(localDate)
    print(listOfAllDates[-4:-1])
    

    if ((listOfAllDates[-1] == localDate and not localDate == processedDay) or
        (not listOfAllDates[-1] == localDate 
             and not listOfAllDates[-1] == processedDay) or
        (not listOfAllDates[-1] == localDate and processedDay == '')
       ):
        print( 'performing day analysis')
        try:
            processedDayDataIndex = listOfAllDates.index( processedDay )
        except ValueError:
            if processedDay == '':
                processedDayDataIndex = 0
            else:
                print('Data for last processed Day does not exist')
                break
        else:
            day = listOfAllDates[processedDayDataIndex]
            fileName = 'CLAWS_analysis_output_' + day
            if os.path.isfile( fileName ):
                figureFileNames, variablesOfInterest = autoAnalysis(
                    day,
                    day,
                    pathToData,
                    dataType='online',
                    defRundNum=('','') )
                writeAnalysisToFile(figureFileNames, 
                                    pathTofigures,
                                    variablesOfInterest,
                                    fileName)
            processedDay = listOfAllDates[ processedDayDataIndex + 1 ]
            print( processedDay, ' day analysis processed')
    try:
        a=sundays.index(localDate)
    except ValueError:
        pass
    else:
        print("make week analysis")
        sundayOfTheYear = calender.index(localDate)
        mondayOfTheYear = sundayOfTheYear - 7
        weekOfData = []
        for i in range(7):
            try:
                d = listOfAllDates.index(calender[mondayOfTheYear + 0] )
            except ValueError:
                pass
            else:
                weekOfData.append(d)
        if not weekOfData == []:
            fileName = ('CLAWS_analysis_output_' 
                        + weekOfData[0] 
                        + '--'
                        + weekOfData[-1])
            
            if os.path.isfile( fileName ):
                figureFileNames, variablesOfInterest = autoAnalysis(
                    weekOfData[0],
                    weekOfData[-1],
                    pathToData,
                    dataType='online',
                    defRundNum=('','') )
                writeAnalysisToFile(figureFileNames, 
                                    pathTofigures,
                                    variablesOfInterest,
                                    fileName)
        else:
            print('No data this week')
