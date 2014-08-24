# catalog of amino acid protons
class Amino_acid:
  def __init__(name, dict_protons):
    self.name = name
    self.dict_protons = dict_protons

name = 'Alanine' # symbol ?
dict_protons = dict()
dict_protons[name] = dict()
dict_protons['Alanine']['H'] = 0.0 # chemical shift of HN in Alanine ?
dict_protons['Alanine']['HA'] = 0.0 # chemical shift of HA in Alanine ?
dict_protons['Alanine']['MB'] = 0.0

dict_protons['Arginine']['H'] = 0.0
dict_protons['Arginine']['HA'] = 0.0
dict_protons['Arginine']['HB_2'] = 0.0
dict_protons['Arginine']['HB_3'] = 0.0
dict_protons['Arginine']['QG'] = 0.0
dict_protons['Arginine']['QD'] = 0.0
dict_protons['Arginine']['HE'] = 0.0

dict_protons['Asparagine']['H'] = 0.0
dict_protons['Asparagine']['HA'] = 0.0
dict_protons['Asparagine']['HB_2'] = 0.0
dict_protons['Asparagine']['HB_3'] = 0.0
dict_protons['Asparagine']['HD21'] = 0.0
dict_protons['Asparagine']['HD22'] = 0.0

dict_protons['Aspartate']['H'] = 0.0
dict_protons['Aspartate']['HA'] = 0.0
dict_protons['Aspartate']['HB_2'] = 0.0
dict_protons['Aspartate']['HB_3'] = 0.0

dict_protons['Cysteine']['H'] = 0.0
dict_protons['Cysteine']['HA'] = 0.0
dict_protons['Cysteine']['HB_2'] = 0.0
dict_protons['Cysteine']['HB_3'] = 0.0

dict_protons['Glutamine']['H'] = 0.0
dict_protons['Glutamine']['HA'] = 0.0
dict_protons['Glutamine']['HB_2'] = 0.0
dict_protons['Glutamine']['HB_3'] = 0.0
dict_protons['Glutamine']['QG'] = 0.0
dict_protons['Glutamine']['HE21'] = 0.0
dict_protons['Glutamine']['HE22'] = 0.0

