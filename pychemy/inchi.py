import re
from decimal import Decimal
from inchi_converter import convert_inchi_to_formula
from adducts import positive_mode_adducts, negative_mode_adducts
import molmass


class Formula_Base(object):
  H_ADDUCT = 1.0073

  def __init__(self, formula=None):
    self.__formula = formula

  @property
  def formula(self):
    return molmass.Formula(self.__formula)

  def __str__(self):
    return self.formula

  @property
  def monoisotopic_mass(self):
    """
    >>> a = Formula_Base('C4H4Na2O4')
    >>> print a.monoisotopic_mass
    161.990497957
    """
    return self.formula.isotope.mass

  def positive_adduct(self, what):
    """
    >>> a = Formula_Base('C4H4Na2O4')
    >>> print a.positive_adduct('M+H')
    162.997797957
    """

    class StupidPythonClosure:
      mz = None
    o = StupidPythonClosure()

    def f(name, mul, off):
      if name == what:
        o.mz = self.monoisotopic_mass*mul+off
    positive_mode_adducts(f)
    return o.mz

  def negative_adduct(self, what):
    """
    >>> a = Formula_Base('C4H4Na2O4')
    >>> print a.negative_adduct('M-H')
    160.983197957
    """

    class StupidPythonClosure:
      mz = None
    o = StupidPythonClosure()

    def f(name, mul, off):
      if name == what:
        o.mz = self.monoisotopic_mass*mul+off
    negative_mode_adducts(f)
    return o.mz


class InChI(Formula_Base):
  """
  >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
  >>> print a.monoisotopic_mass
  161.990497957
  """

  def __init__(self, inchi, formula=None):
    if not inchi.startswith("InChI="):
      raise Exception("Invalid InChI %s, must start with InChI=" % (inchi,))
    self.__inchi = inchi
    self.__formula = formula
    inchi = inchi[6:]
    self.__layers = inchi.split('/')[1:]

  @property
  def inchi(self):
    """
    >>> a = InChI('InChI=1S/H2O/h1H2')
    >>> print a.inchi
    InChI=1S/H2O/h1H2
    """
    return self.__inchi

  @property
  def formula(self):
    """
    >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
    >>> print a.formula
    C4H4Na2O4
    >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2', 'H2O')
    >>> print a.formula
    H2O
    """
    if self.__formula is None:
      self.__formula = convert_inchi_to_formula(self.__inchi)
    return molmass.Formula(self.__formula)

  @property
  def spectrum(self):
    """
    >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
    >>> s = a.spectrum
    >>> print s
    Relative mass    Fraction %      Intensity
    161.99050         94.816904     100.000000
    162.99385          4.102055       4.326291
    162.99471          0.144473       0.152370
    162.99677          0.043621       0.046005
    163.99474          0.779393       0.821997
    163.99721          0.066550       0.070188
    163.99807          0.006250       0.006592
    164.00013          0.001887       0.001990
    164.00099          0.000017       0.000018
    164.00305          0.000004       0.000004
    164.99810          0.033719       0.035562
    164.99896          0.000891       0.000939
    165.00076          0.000838       0.000884
    165.00348          0.000008       0.000008
    165.00435          0.000001       0.000001
    165.99899          0.002402       0.002534
    166.00145          0.000547       0.000577
    166.00438          0.000004       0.000004
    >>> print ['%.4f %.6f' % (v[0],v[1]) for v in s.values()]
    ['161.9905 0.948169', '162.9939 0.041021', '162.9947 0.001445', '162.9968 0.000436', '163.9947 0.007794', '163.9972 0.000666', '163.9981 0.000063', '164.0001 0.000019', '164.0010 0.000000', '164.0031 0.000000', '164.9981 0.000337', '164.9990 0.000009', '165.0008 0.000008', '165.0035 0.000000', '165.0043 0.000000', '165.0064 0.000000', '165.9990 0.000024', '166.0015 0.000005', '166.0044 0.000000', '166.0052 0.000000', '166.0070 0.000000', '167.0053 0.000000', '167.0077 0.000000']
    """
    return self.formula.high_resolution_spectrum(3, minfract=1e-6)

if __name__ == '__main__':
  import doctest
  doctest.testmod()

