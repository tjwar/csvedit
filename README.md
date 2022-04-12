# csvedit
for manipulating ecommerce inventory csv files for different import requirements. e.g. generating explicit parent/child relationships, transforming, collating values 

reqs:

from csv import *
import copy
import sys
import re
import pandas as pd
from pandas import *
from Bio.Align import PairwiseAligner
from colour import Color
from pprint import pprint
import datetime
