"""
Converts InChI to formula. Depends on/uses openbabel.
"""

import subprocess
import tempfile
import os


def convert_inchi_to_formula(inchi_string):
  output = None
  with tempfile.NamedTemporaryFile(mode='w+b', delete=False) as inf, \
       tempfile.NamedTemporaryFile(mode='w+b', delete=False) as outf:

    # get outfile name and make input file
    outf.close()
    inf.write(inchi_string)
    inf.close()
    inf.name

    # use openbabel to convert
    subprocess.check_output(['babel', '-iinchi', inf.name, '-oreport', outf.name], stderr=subprocess.STDOUT)
    f = open(outf.name, 'r')
    output = f.read().split('\n')
    f.close()

    # cleanup
    os.unlink(inf.name)
    os.unlink(outf.name)

  formula = None
  if output is not None:
    k = 'FORMULA: '
    f = [l.strip() for l in output if l.startswith(k)]
    if len(f) > 0:
      formula = f[0][len(k):]

  return formula 
