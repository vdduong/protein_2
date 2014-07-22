# read the experimental peak list of protein into input

import re

class Peak_noesy(object):
  def __init__(self, nb_peak, cs_H1, cs_H2, cs_C, volume):
    self.nb_peak = nb_peak
    self.cs_H1 = cs_H1
    self.cs_H2 = cs_H2
    self.cs_C = cs_C
    self.volume = volume

dict_peak_list = dict()

file = file.open('/c13.peaks')
# ...
for line in file:
  line = line.rstrip('\n\r')
  line = line.re(' +')
  nb_peak = int(line[1])
  cs_H1 = float(line[2])
  cs_H2 = float(line[3])
  cs_C = float(line[4])
  volume = float(line[7])
  dict_peak_list[nb_peak] = Peak_noesy(nb_peak, cs_H1, cs_H2, cs_C, volume)

