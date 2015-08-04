import unittest, tempfile, os
from pychemy.peptides import *
from pychemy.amino_acids import AMINO_ACIDS

class Peptides_Testing(unittest.TestCase):   


###############################

  def test_get_next_AA_with_valid_single_AAs(self):
    self.assertEqual(get_next_AA('A'), AMINO_ACIDS['A'])
    self.assertEqual(get_next_AA('C'), AMINO_ACIDS['C'])
    self.assertEqual(get_next_AA('D'), AMINO_ACIDS['D'])
    self.assertEqual(get_next_AA('E'), AMINO_ACIDS['E'])
    self.assertEqual(get_next_AA('F'), AMINO_ACIDS['F'])
    self.assertEqual(get_next_AA('G'), AMINO_ACIDS['G'])
    self.assertEqual(get_next_AA('H'), AMINO_ACIDS['H'])
    self.assertEqual(get_next_AA('I'), AMINO_ACIDS['I'])
    self.assertEqual(get_next_AA('K'), AMINO_ACIDS['K'])
    self.assertEqual(get_next_AA('L'), AMINO_ACIDS['L'])
    self.assertEqual(get_next_AA('M'), AMINO_ACIDS['M'])
    self.assertEqual(get_next_AA('N'), AMINO_ACIDS['N'])
    self.assertEqual(get_next_AA('P'), AMINO_ACIDS['P'])
    self.assertEqual(get_next_AA('Q'), AMINO_ACIDS['Q'])
    self.assertEqual(get_next_AA('R'), AMINO_ACIDS['R'])
    self.assertEqual(get_next_AA('S'), AMINO_ACIDS['S'])
    self.assertEqual(get_next_AA('T'), AMINO_ACIDS['T'])
    self.assertEqual(get_next_AA('V'), AMINO_ACIDS['V'])
    self.assertEqual(get_next_AA('W'), AMINO_ACIDS['W'])
    self.assertEqual(get_next_AA('Y'), AMINO_ACIDS['Y'])

  def test_get_next_AA_with_valid_AAs(self):
    self.assertEqual(get_next_AA('AA'), AMINO_ACIDS['A'])
    self.assertEqual(get_next_AA('CA'), AMINO_ACIDS['C'])
    self.assertEqual(get_next_AA('DA'), AMINO_ACIDS['D'])
    self.assertEqual(get_next_AA('EA'), AMINO_ACIDS['E'])
    self.assertEqual(get_next_AA('FA'), AMINO_ACIDS['F'])
    self.assertEqual(get_next_AA('GA'), AMINO_ACIDS['G'])
    self.assertEqual(get_next_AA('HA'), AMINO_ACIDS['H'])
    self.assertEqual(get_next_AA('IA'), AMINO_ACIDS['I'])
    self.assertEqual(get_next_AA('KA'), AMINO_ACIDS['K'])
    self.assertEqual(get_next_AA('LA'), AMINO_ACIDS['L'])
    self.assertEqual(get_next_AA('MA'), AMINO_ACIDS['M'])
    self.assertEqual(get_next_AA('NA'), AMINO_ACIDS['N'])
    self.assertEqual(get_next_AA('PA'), AMINO_ACIDS['P'])
    self.assertEqual(get_next_AA('QA'), AMINO_ACIDS['Q'])
    self.assertEqual(get_next_AA('RA'), AMINO_ACIDS['R'])
    self.assertEqual(get_next_AA('SA'), AMINO_ACIDS['S'])
    self.assertEqual(get_next_AA('TA'), AMINO_ACIDS['T'])
    self.assertEqual(get_next_AA('VA'), AMINO_ACIDS['V'])
    self.assertEqual(get_next_AA('WA'), AMINO_ACIDS['W'])
    self.assertEqual(get_next_AA('YA'), AMINO_ACIDS['Y'])

  def test_get_next_AA_with_valid_single_mods(self):
    self.assertEqual(get_next_AA('C(carbamidomethyl)'), AMINO_ACIDS['C(carbamidomethyl)'])
    self.assertEqual(get_next_AA('M(ox)'), AMINO_ACIDS['M(ox)'])
    self.assertEqual(get_next_AA('[K_C13N15]'), AMINO_ACIDS['[K_C13N15]'])
    self.assertEqual(get_next_AA('[R_C13N15]'), AMINO_ACIDS['[R_C13N15]'])

  def test_get_next_AA_with_valid_mods(self):
    self.assertEqual(get_next_AA('C(carbamidomethyl)A'), AMINO_ACIDS['C(carbamidomethyl)'])
    self.assertEqual(get_next_AA('M(ox)A'), AMINO_ACIDS['M(ox)'])
    self.assertEqual(get_next_AA('[K_C13N15]A'), AMINO_ACIDS['[K_C13N15]'])
    self.assertEqual(get_next_AA('[R_C13N15]A'), AMINO_ACIDS['[R_C13N15]'])

  def test_get_next_AA_with_single_invalid_AA(self):
    self.assertRaises(Exception, get_next_AA, 'B')

  def test_get_next_AA_with_invalid_AA(self):
    self.assertRaises(Exception, get_next_AA, 'BA')

###############################
  
  def test_residues_from_sequence_with_single_AA(self):
    self.assertEqual(residues_from_sequence('A'), ['A'])
    self.assertEqual(residues_from_sequence('C'), ['C'])
    self.assertEqual(residues_from_sequence('D'), ['D'])
    self.assertEqual(residues_from_sequence('E'), ['E'])
    self.assertEqual(residues_from_sequence('F'), ['F'])
    self.assertEqual(residues_from_sequence('G'), ['G'])
    self.assertEqual(residues_from_sequence('H'), ['H'])
    self.assertEqual(residues_from_sequence('I'), ['I'])
    self.assertEqual(residues_from_sequence('K'), ['K'])
    self.assertEqual(residues_from_sequence('L'), ['L'])
    self.assertEqual(residues_from_sequence('M'), ['M'])
    self.assertEqual(residues_from_sequence('N'), ['N'])
    self.assertEqual(residues_from_sequence('P'), ['P'])
    self.assertEqual(residues_from_sequence('Q'), ['Q'])
    self.assertEqual(residues_from_sequence('R'), ['R'])
    self.assertEqual(residues_from_sequence('S'), ['S'])
    self.assertEqual(residues_from_sequence('T'), ['T'])
    self.assertEqual(residues_from_sequence('V'), ['V'])
    self.assertEqual(residues_from_sequence('W'), ['W'])
    self.assertEqual(residues_from_sequence('Y'), ['Y'])

  def test_residues_from_sequence_with_valid_sequence(self):
    self.assertEqual(residues_from_sequence('AA'), ['A','A'])
    self.assertEqual(residues_from_sequence('CA'), ['C','A'])
    self.assertEqual(residues_from_sequence('DA'), ['D','A'])
    self.assertEqual(residues_from_sequence('EA'), ['E','A'])
    self.assertEqual(residues_from_sequence('FA'), ['F','A'])
    self.assertEqual(residues_from_sequence('GA'), ['G','A'])
    self.assertEqual(residues_from_sequence('HA'), ['H','A'])
    self.assertEqual(residues_from_sequence('IA'), ['I','A'])
    self.assertEqual(residues_from_sequence('KA'), ['K','A'])
    self.assertEqual(residues_from_sequence('LA'), ['L','A'])
    self.assertEqual(residues_from_sequence('MA'), ['M','A'])
    self.assertEqual(residues_from_sequence('NA'), ['N','A'])
    self.assertEqual(residues_from_sequence('PA'), ['P','A'])
    self.assertEqual(residues_from_sequence('QA'), ['Q','A'])
    self.assertEqual(residues_from_sequence('RA'), ['R','A'])
    self.assertEqual(residues_from_sequence('SA'), ['S','A'])
    self.assertEqual(residues_from_sequence('TA'), ['T','A'])
    self.assertEqual(residues_from_sequence('VA'), ['V','A'])
    self.assertEqual(residues_from_sequence('WA'), ['W','A'])
    self.assertEqual(residues_from_sequence('YA'), ['Y','A'])

  def test_residues_from_sequence_with_invalid_first_AA(self):
    self.assertRaises(Exception, residues_from_sequence, 'BA')

  def test_residues_from_sequence_with_invalid_internal_AA(self):
    self.assertRaises(Exception, residues_from_sequence, 'ABA')

  def test_residues_from_sequence_with_invalid_terminal_AA(self):
    self.assertRaises(Exception, residues_from_sequence, 'AB')
