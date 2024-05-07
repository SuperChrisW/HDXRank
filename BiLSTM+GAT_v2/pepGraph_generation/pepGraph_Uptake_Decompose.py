# pyhdx environment
import pandas as pd
import numpy as np
from scipy.optimize import lsq_linear

from pyhdx import read_dynamx, HDXMeasurement
from pyhdx.process import filter_peptides, apply_control, correct_d_uptake
from pyhdx.models import Coverage
from pyhdx.batch_processing import StateParser
from pyhdx.fitting import fit_d_uptake
from pyhdx.config import cfg

root_dir = ''
summary_df = pd.read_excel('/home/lwang/models/HDX_LSTM/data/RTT_BCD/hdock/rtt+bcd_hdock.xlsx', sheet_name='Sheet1') 

fpath = '/home/lwang/models/HDX_LSTM/data/RTT_BCD/hdock/HDX_files/PXD023434.csv'
data = read_dynamx(fpath)
data['rfu'] = data['%d'] / 100

peptides = data[(data['protein'] == 'bcdFL') & (data['state'] == 'Bcd')]

peptides_corrected = correct_d_uptake(peptides, drop_first=1, d_percentage=90.0)

sequence = "MAVLCGVCGIKEFKYKCPRCLVQTCSLECSKKHKTRDNCSGQTHDPKEYISSEALKQADDDKHERNAYVQRDYNYLTQLKRMVHVQKMDARMKNKRVLGPVGGHNSNFKKRRYDIDEDDRDSTECQRIIRRGVNCLMLPKGMQRSSQNRSKWDKTMDLFVWSVEWILCPMQEKGEKKELFKHVSHRIKETDFLVQGMGKNVFQKCCEFYRLAGTSSCIEGEDGSETKEERTQILQKSGLKFYTKTFPYNTTHIMDSKKLVELAIHEKCIGELLKNTTVIEFPTIFVAMTEADLPEGYEVLHQEPRPLEHTSTLNKFIDNAREEEDAEEDSQPTEEPVQKETQDASDSDSDSDDDYNPGLSMDFLTA"
temperature, pH = 273.15 + 25, 7.5

hdxm = HDXMeasurement(
    peptides_corrected, sequence=sequence, pH=pH, temperature=temperature, name="My HDX measurement"
)

hdx_t = hdxm[-1]
fit_result_1 = fit_d_uptake(hdx_t, r1=1.0, repeats=20)