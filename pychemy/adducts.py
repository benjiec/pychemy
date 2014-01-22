
# adducts specified as (name, mz_mul, mz_off) tuples. based on monoisotopic
# neutral mz.

Positive_Adducts = (
  ('M+3H',           0.33,  1.0073),
  ('M+3H+Na',        0.33,  8.3346),
  ('M+H+2Na',        0.33, 15.7662),
  ('M+3Na',          0.33, 22.9892),
  ('M+2H',           0.50,  1.0073),
  ('M+H+NH4',        0.50,  9.5206),
  ('M+H+Na',         0.50, 11.9982),
  ('M+H+K',          0.50, 19.9852),
  ('M+ACN+2H',       0.50, 21.5256),
  ('M+2Na',          0.50, 22.9892),
  ('M+2ACN+2H',      0.50, 42.0338),
  ('M+3ACN+2H',      0.50, 62.5471),
  ('M+H',            1.00,  1.0073),
  ('M+NH4',          1.00, 18.0338),
  ('M+Na',           1.00, 22.9892),
  ('M+CH3OH+H',      1.00, 33.0335),
  ('M+K',            1.00, 38.9632),
  ('M+ACN+H',        1.00, 42.0338),
  ('M+2Na-H',        1.00, 44.9712),
  ('M+IsoProp+H',    1.00, 61.0653),
  ('M+ACN+Na',       1.00, 64.0157),
  ('M+2K-H',         1.00, 76.9190),
  ('M+DMSO+H',       1.00, 79.0212),
  ('M+2ACN+H',       1.00, 83.0604),
  ('M+IsoProp+Na+H', 1.00, 84.0551),
  ('2M+H',           2.00,  1.0073),
  ('2M+NH4',         2.00, 18.0338),
  ('2M+Na',          2.00, 22.9892),
  ('2M+3H2O+2H',     2.00, 28.0231),
  ('2M+K',           2.00, 38.9632),
  ('2M+ACN+H',       2.00, 42.0338),
  ('2M+ACN+Na',      2.00, 64.0158),
  ('M+H-H2O',        1.00, -17.0038),
)

Negative_Adducts = (
  ('M-3H',     0.33,  -1.0073),
  ('M-2H',     0.50,  -1.0073),
  ('M-H',      1.00,  -1.0073),
  ('M-H2O-H',  1.00, -19.0184),
  ('M+Na-2H',  1.00,  20.9747),
  ('M+Cl',     1.00,  34.9694),
  ('M+K-2H',   1.00,  36.9486),
  ('M+FA-H',   1.00,  44.9982),
  ('M+Hac-H',  1.00,  59.0139),
  ('M+Br',     1.00,  78.9189),
  ('M+TFA-H',  1.00, 112.9856),
  ('2M-H',     2.00,  -1.0073),
  ('2M+FA-H',  2.00,  44.9982),
  ('2M+Hac-H', 2.00,  59.0139),
  ('3M-H',     3.00,  -1.0073),
)

def _iterate_adducts(adduct_list, func):
  for adduct in adduct_list:
    func(adduct[0], adduct[1], adduct[2])

def positive_mode_adducts(func):
  """
  Iterates through positive mode adducts. Calls func for each adduct with
  adduct name, mass multiplier, and mass offset as arguments.

  >>> def a(name, mul, off):
  ...   if name == 'M+H': print 'x*%.2f+(%.4f)' % (mul, off)
  >>> positive_mode_adducts(a)
  x*1.00+(1.0073)
  """
  _iterate_adducts(Positive_Adducts, func)

def negative_mode_adducts(func):
  """
  Iterates through negative mode adducts. Calls func for each adduct with
  adduct name, mass multiplier, and mass offset as arguments.

  >>> def a(name, mul, off):
  ...   if name == 'M-H': print 'x*%.2f+(%.4f)' % (mul, off)
  >>> negative_mode_adducts(a)
  x*1.00+(-1.0073)
  """
  _iterate_adducts(Negative_Adducts, func)


if __name__ == '__main__':
  import doctest
  doctest.testmod()

