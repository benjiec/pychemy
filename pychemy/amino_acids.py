from pychemy.molmass import Formula
from pychemy.elements import ELEMENTS, Isotope


class Amino_Acid(object):
  def __init__(self, residue = '', mod = '', formula = '', losses = []):
    self.residue = residue
    self.mod = mod
    self.formula = formula
    if formula:
      self.mass = Formula(formula).isotope.mass
    else:
      self.mass = 0.0
    self.losses = losses
    AMINO_ACIDS.add(self)

  def __str__(self):
    return self.mod


class Loss(object):
  def __init__(self, name, formula):
    self.display_name = name
    self.formula = formula
    self.mass = Formula(formula).isotope.mass
    LOSSES.add(self)

    def __str__(self):
      return self.formula


class Amino_Acid_Registry(object):
  def __init__(self):
    self._list = []
    self._dict = {}
    
  def __add__(self, amino_acid):
    if amino_acid.mod in self._dict.keys():
      raise KeyError
    else:
      self._list.append(amino_acid)
      self._dict[amino_acid.mod] = amino_acid
  
  def __str__(self):
    return "[%s]" % ", ".join(aa.mod for aa in self._list)
    
  def __iter__(self):
    return iter(self._list)
    
  def __len__(self):
    return len(self._list)

  def __getitem__(self, key):
    return self._dict[key]
  
  def add(self, amino_acid):
    self.__add__(amino_acid)  


class Loss_Registry(object):
  def __init__(self):
    self._list = []
    self._dict = {}
  
  def __add__(self, loss):
    if loss.formula in self._dict.keys():
      raise KeyError
    else:
      self._list.append(loss)
      self._dict[loss.formula] = loss
  
  def __str__(self):
    return "[%s]" % ", ".join(loss.formula for loss in self._list)

  def __iter__(self):
    return iter(self._list)

  def __len__(self):
    return len(self._list)

  def __getitem__(self, key):
    try:
      return self._dict[key]
    except:
      raise KeyError

  def add(self, loss):
    self.__add__(loss)


AMINO_ACIDS = Amino_Acid_Registry()
LOSSES = Loss_Registry()


Loss(name = 'H_2O', formula = 'OH2')
Loss(name = 'NH_3', formula = 'NH3')
Loss(name = 'CO', formula = 'CO')


Amino_Acid(residue = 'A', mod = 'A', formula = 'C3H5N1O1')
Amino_Acid(residue = 'R', mod = 'R', formula = 'C6H12N4O1')
Amino_Acid(residue = 'N', mod = 'N', formula = 'C4H6N2O2')
Amino_Acid(residue = 'D', mod = 'D', formula = 'C4H5N1O3')
Amino_Acid(residue = 'C', mod = 'C', formula = 'C3H5N1O1S1')
Amino_Acid(residue = 'C', mod = 'C(carbamidomethyl)', formula = 'C5H8N2O2S1')
Amino_Acid(residue = 'Q', mod = 'Q', formula = 'C5H8N2O2')
Amino_Acid(residue = 'E', mod = 'E', formula = 'C5H7N1O3')
Amino_Acid(residue = 'G', mod = 'G', formula = 'C2H3N1O1')
Amino_Acid(residue = 'H', mod = 'H', formula = 'C6H7N3O1')
Amino_Acid(residue = 'I', mod = 'I', formula = 'C6H11N1O1')
Amino_Acid(residue = 'L', mod = 'L', formula = 'C6H11N1O1')
Amino_Acid(residue = 'K', mod = 'K', formula = 'C6H12N2O1')
Amino_Acid(residue = 'M', mod = 'M', formula = 'C5H9N1O1S1')
Amino_Acid(residue = 'M', mod = 'M(ox)', formula = 'C5H9N1O2S1')

Amino_Acid(residue = 'F', mod = 'F', formula = 'C9H9N1O1')
Amino_Acid(residue = 'P', mod = 'P', formula = 'C5H7N1O1')
Amino_Acid(residue = 'S', mod = 'S', formula = 'C3H5N1O2', losses = [LOSSES['OH2']])
Amino_Acid(residue = 'T', mod = 'T', formula = 'C4H7N1O2', losses = [LOSSES['OH2']])
Amino_Acid(residue = 'W', mod = 'W', formula = 'C11H10N2O1')
Amino_Acid(residue = 'Y', mod = 'Y', formula = 'C9H9N1O2')
Amino_Acid(residue = 'V', mod = 'V', formula = 'C5H9N1O1')
Amino_Acid(residue = 'N-term', mod = 'N-term', formula = '')
Amino_Acid(residue = 'C-term', mod = 'C-term', formula = 'OH2')

# Multiplex labeling reagents
Amino_Acid(residue = 'N-term', mod = 'TMT10plex', formula = 'C8[13C]4H20N[15N]O2')
Amino_Acid(residue = 'K', mod = 'K(TMT10plex)', formula = 'C14H33N3O3[13C]4[15N]')

# mod format for heavy-labeled peptide ordering
Amino_Acid(residue = 'K', mod = '[K_C13N15]', formula = '[13C]6H12[15N]2O1')
Amino_Acid(residue = 'R', mod = '[R_C13N15]', formula = '[13C]6H12[15N]4O1')
