import unittest, tempfile, os
from pychemy.peptides import *
from pychemy.amino_acids import AMINO_ACIDS

class Peptides_Testing(unittest.TestCase):   


###############################

  def test_get_next_aa_with_valid_single_aa(self):
    for key in AMINO_ACIDS:
      if key.mod != 'N-term' and key.mod != 'C-term':
        self.assertEqual(get_next_aa(key.mod), AMINO_ACIDS[key.mod])

  def test_get_next_aa_with_valid_aa_sequence(self):
    for key in AMINO_ACIDS:
      if key.mod != 'N-term' and key.mod != 'C-term':
        print type(key.mod), key.mod
        self.assertEqual(get_next_aa(key.mod + 'A'), AMINO_ACIDS[key.mod])

  def test_get_next_aa_with_single_invalid_aa(self):
    self.assertRaises(Exception, get_next_aa, 'B')

  def test_get_next_aa_with_invalid_aa(self):
    self.assertRaises(Exception, get_next_aa, 'BA')

###############################
  
  def test_residues_from_sequence_with_single_aa(self):
    for key in AMINO_ACIDS:
      if key.mod != 'N-term' and key.mod != 'C-term':
        self.assertEqual(residues_from_sequence(key.mod), [key.mod])

  def test_residues_from_sequence_with_valid_sequence(self):
    for key in AMINO_ACIDS:
      if key.mod != 'N-term' and key.mod != 'C-term':
        self.assertEqual(residues_from_sequence(key.mod + 'A'), [key.mod, 'A'])

  def test_residues_from_sequence_with_complex_sequences(self):
    for key1 in AMINO_ACIDS:
      if key1.residue not in ('N-term', 'C-term'):
        for key2 in AMINO_ACIDS:
          if key2.residue not in ('N-term', 'C-term'):
            for key3 in AMINO_ACIDS: 
              if key3.residue not in ('N-term', 'C-term'):
                self.assertEqual(residues_from_sequence(key1.mod + key2.mod + key3.mod), [key1.mod, key2.mod, key3.mod])

  def test_residues_from_sequence_with_invalid_first_aa(self):
    self.assertRaises(Exception, residues_from_sequence, 'BA')

  def test_residues_from_sequence_with_invalid_internal_aa(self):
    self.assertRaises(Exception, residues_from_sequence, 'ABA')

  def test_residues_from_sequence_with_invalid_terminal_aa(self):
    self.assertRaises(Exception, residues_from_sequence, 'AB')
