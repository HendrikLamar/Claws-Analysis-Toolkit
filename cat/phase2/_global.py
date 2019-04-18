#picos = ['Top_Forward', 'Top_Backward', 'Bottom_Forward', 'Bottom_Backward']
picos = ['FWD_Z0', 'FWD_Z1', 'FWD_Z2', 'FWD_Z3', 'BWD_Z0', 'BWD_Z1', 'BWD_Z2', 'BWD_Z3']
channels = ['A', 'B', 'C','D']
items_channel = {'1pe' : '1peNormalized',
        '1peVStotal' : '1peVStotal', 'mips' : 'mipAbsolute', 'reco' : 'wfReco', 'cdf' : 'cdf'}
items_pico = {'run' : 'beast_runNumber', 'subRun' : 'beast_subRunNumber', 'ts' : 'timestamp',
        'lerID' : 'LER_injID', 'herID' : 'HER_injID', 'lerCur' : 'LER_curr', 'herCur' : 'HER_curr',
        'lerBGate' : 'LER_bGate', 'herBGate' : 'HER_bGate'}

path_to_data = ""
