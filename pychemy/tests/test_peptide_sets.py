import unittest, tempfile, os
from pychemy.peptide_sets import *
import random

import numpy as np

AA_dict = {1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'L', 10:'K',
           11:'M', 12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y'}

def gen_protein(num_AA = 20):
  out = ''
  for i in range(num_AA - 1):
    out += AA_dict[round(0.5 + 20*random.random())]
  out += 'R'
  return out
   



class Peptide_Set_Testing(unittest.TestCase):   


###############################

  def test_check_AA_with_valid_sequence(self):
    self.assertEqual(check_AA('ACDEFGHIKLMNPQRSTVWY'), True)

  def test_check_AA_with_invalid_sequence(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AAB', check_AA, 'AAB')

###############################

  def test_check_proteins_with_valid_sequence(self):
    self.assertEqual(check_proteins(['ACDEFGHIKLMNPQRSTVWY','AA']), True)

  def test_check_proteins_with_invalid_sequence(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AAB', check_proteins, ['AA', 'AAB'])

###############################
    
  def test_get_peptides_with_valid_exaples(self):
    self.assertSetEqual(set(get_peptides('AAR')), set(['AAR']))
    self.assertSetEqual(set(get_peptides('AARAAK')), set(['AAR', 'AAK']))
    self.assertSetEqual(set(get_peptides('AARAAKAA')), set(['AAR', 'AAK']))

  def test_get_peptides_with_invalid_examples(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AAB', get_peptides, 'AAB')
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AARAAB', get_peptides, 'AARAAB')  

###############################

  def test_all_peptides_with_single_valid_protein(self):
    self.assertSetEqual(set(all_peptides(['AAR'])), set(['AAR']))
    self.assertSetEqual(set(all_peptides(['AARAAK'])), set(['AAR', 'AAK']))

  def test_all_peptides_with_multiple_proteins_valid(self):
    self.assertSetEqual(set(all_peptides(['AAR', 'AAK'])), set(['AAR', 'AAK']))
    self.assertSetEqual(set(all_peptides(['AARACR', 'AAKACK'])), set(['AAR', 'ACR', 'AAK', 'ACK']))
    self.assertSetEqual(set(all_peptides(['AARACRAA', 'AAKACK'])), set(['AAR', 'ACR', 'AAK', 'ACK']))

  def test_all_peptides_with_multiple_proteins_invalid(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AARAAB', all_peptides, ['AARAAB', 'AAKACK'])
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: AARAAB', all_peptides, ['AARAAB', 'AAKBCK'])

###############################

  def test_distinct_peptides_with_non_duplicated_peptides(self): 
    self.assertSetEqual(set(distinct_peptides(['AARACR', 'AAKACK'])), set(['AAR', 'ACR', 'AAK', 'ACK']))

  def test_distinct_peptides_with_duplicated_peptides(self): 
    self.assertSetEqual(set(distinct_peptides(['AARACR', 'AAKAAR'])), set(['AAR', 'ACR', 'AAK']))

###############################

  def test_unique_identifiers_with_single_protein_valid(self):
    temp = unique_identifiers(['AAR'])
    self.assertIn(['AAR', 'AAR'], temp)

    temp = unique_identifiers(['AARAAK'])
    self.assertIn(['AARAAK', 'AAR'], temp)
    self.assertIn(['AARAAK', 'AAK'], temp)

    temp = unique_identifiers(['AARAAKAA'])
    self.assertIn(['AARAAKAA', 'AAR'], temp)
    self.assertIn(['AARAAKAA', 'AAK'], temp)

  def test_unique_identifiers_with_single_protein_invalid(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABR', unique_identifiers, ['ABR'])
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABRAAK', unique_identifiers, ['ABRAAK'])
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABRAAKAA', unique_identifiers, ['ABRAAKAA'])
  
  def test_unique_identifiers_with_multiple_proteins_valid(self):
    temp = unique_identifiers(['AAR', 'AAK'])
    self.assertIn(['AAR', 'AAR'], temp)
    self.assertIn(['AAK', 'AAK'], temp)

    temp = unique_identifiers(['AARAAK', 'ACRACK'])
    self.assertIn(['AARAAK', 'AAR'], temp)
    self.assertIn(['AARAAK', 'AAK'], temp)
    self.assertIn(['ACRACK', 'ACR'], temp)
    self.assertIn(['ACRACK', 'ACK'], temp)

    temp = unique_identifiers(['AARAAKAA', 'ACRACKAA'])
    self.assertIn(['AARAAKAA', 'AAR'], temp)
    self.assertIn(['AARAAKAA', 'AAK'], temp)
    self.assertIn(['ACRACKAA', 'ACR'], temp)
    self.assertIn(['ACRACKAA', 'ACK'], temp)

  def test_unique_identifiers_with_multiple_proteins_invalid(self):
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABR', unique_identifiers, ['ABR', 'ACRACK'])
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABRAAK', unique_identifiers, ['ABRAAK', 'ACRACK'])
    self.assertRaisesRegexp(Exception, 'Invalid amino acid in seq: ABRAAKAA', unique_identifiers, ['ABRAAKAA', 'ACRACK'])

###############################

  def test_peptides_per_protein_with_single_protein(self):
    temp = peptides_per_protein(unique_identifiers(['AAR']))
    self.assertIn(['AAR', ['AAR']], temp)

    temp = peptides_per_protein(unique_identifiers(['AARAAK']))
    self.assertEqual('AARAAK', temp[0][0])
    self.assertEqual(len(temp[0][1]), 2)
    self.assertIn('AAR', temp[0][1])
    self.assertIn('AAK', temp[0][1])
    
    temp = peptides_per_protein(unique_identifiers(['AARAAKAA']))
    self.assertEqual('AARAAKAA', temp[0][0])
    self.assertEqual(len(temp[0][1]), 2)
    self.assertIn('AAR', temp[0][1])
    self.assertIn('AAK', temp[0][1])
    
  def test_peptides_per_protein_with_multiple_proteins_not_overlapping(self):
    temp = peptides_per_protein(unique_identifiers(['AAR', 'AAK']))
    self.assertEquals(len(temp), 2)
    self.assertIn(['AAR', ['AAR']], temp)
    self.assertIn(['AAK', ['AAK']], temp)

    temp = peptides_per_protein(unique_identifiers(['AARAAK', 'ACRACK']))
    self.assertEqual(len(temp), 2)
    for item in temp:
      self.assertEqual(len(item[1]), 2)
      for pep in item[1]:
        self.assertTrue(re.search(pep, item[0]))
    
  def test_peptides_per_protein_with_multiple_proteins_overlapping_properly_includes(self):
    temp = peptides_per_protein(unique_identifiers(['AARAAKCCR', 'ACRACKCCR']))
    self.assertEqual(len(temp), 2)
    for item in temp:
      self.assertEqual(len(item[1]), 2)
      for pep in item[1]:
        self.assertTrue(re.search(pep, item[0]))

  def test_peptides_per_protein_with_multiple_proteins_overlapping_properly_excludes(self):
    temp = peptides_per_protein(unique_identifiers(['AARAAKCCR', 'ACRACKCCR']))
    self.assertEqual(len(temp), 2)
    for item in temp:
      self.assertFalse('CCR' in item[1])

###############################
  
  def test_min_ui_count_with_one_protein(self):
    self.assertEqual(min_ui_count(['AAR']), 1)
    self.assertEqual(min_ui_count(['AARAAK']), 2)
    self.assertEqual(min_ui_count(['AARAAKAA']), 2)

  def test_min_ui_count_with_multiple_proteins_non_overlapping(self):
    self.assertEqual(min_ui_count(['AAR', 'AAK']), 1)
    self.assertEqual(min_ui_count(['AARAAK', 'ACR']), 1)
    self.assertEqual(min_ui_count(['AARAAK', 'ACRACK']), 2)
    
  def test_min_ui_count_with_multiple_proteins_overlapping(self):
    self.assertEqual(min_ui_count(['AARAAK', 'ACRAAK']), 1)
    self.assertEqual(min_ui_count(['AARAAK', 'ACRACKAAK']), 1)
    self.assertEqual(min_ui_count(['AARAAK', 'AARAAKACR']), 0)
    
###############################
  
  def test_balanced_sets_with_single_protein(self):
    self.assertEqual(balanced_sets(['AAR']), [['AAR']])

  def run_balanced_sets(self, num_iter = 1, prot_length = 100, num_prot = 100, max_set_size = 10, unique = 5):
    for i in range(num_iter):   
      test_proteins = [gen_protein(prot_length) for j in range(num_prot)]
      prot_count = 0
      for b in balanced_sets(test_proteins, max_set_size = max_set_size, unique = unique):
        self.assertTrue(len(b) <= max_set_size)
        for p in peptides_per_protein(unique_identifiers(b)):
          self.assertTrue(len(p[1]) >= unique or len(b) == 1)
          prot_count += 1
      self.assertEqual(prot_count, num_prot)
   
  def test_balanced_sets(self):
    self.run_balanced_sets()
    self.run_balanced_sets(max_set_size = 25)
    self.run_balanced_sets(unique = 2)

###############################

  def test_SVM_props_with_empty_sequence(self):
    f = peptide_flyability()
    self.assertTrue(np.array_equal(f.SVM_props(''), np.array([])))

  def test_SVM_props_with_nonempty_sequence(self):
    f = peptide_flyability()
    correct = np.array(
    [ 1.30000000e+01,   1.42764901e+03,   7.00000000e+00,   6.00000000e+00,
      2.00000000e+00,   4.00000000e+00,   0.00000000e+00,   4.00000000e+00,
      5.46153846e-02,   5.00000000e-01,  -6.38461538e-01,  -5.24615385e-01,
      8.89230769e+00,   1.60892308e+01,   1.54207692e+01,   1.00000000e+00,
      0.00000000e+00,   1.00000000e+00,   3.00000000e+00,   0.00000000e+00,
      0.00000000e+00,   0.00000000e+00,   1.00000000e+00,   0.00000000e+00,
      1.00000000e+00,   1.00000000e+00,   0.00000000e+00,   3.00000000e+00,
      0.00000000e+00,   0.00000000e+00,   1.00000000e+00,   1.00000000e+00,
      0.00000000e+00,   0.00000000e+00,   0.00000000e+00])
    test = f.SVM_props('SAMPLEPEPTIDE')

    for idx, item in enumerate(correct):
      self.assertAlmostEqual(test[idx], item, places=3)

###############################

  def test_SVM_score_without_complete_feature_vector(self):
    f = peptide_flyability()
    self.assertIsNone(f.SVM_score([]))
    self.assertIsNone(f.SVM_score([1]))
    self.assertIsNone(f.SVM_score([1,2,3]))

  def test_SVM_score_with_complete_feature_vector(self):
    f = peptide_flyability()
    test = f.SVM_score(np.divide(f.SVM_props('SAMPLEPEPTIDE') - f.STEPP_mean, f.STEPP_std))
    self.assertAlmostEqual(test, 1.6285437521)

###############################

  def test_ionization_prob(self):
    f = peptide_flyability()
    self.assertAlmostEqual(f.ionization_prob(0.209218277881), 0.798247804562)
    self.assertAlmostEqual(f.ionization_prob(1.6285437521), 1.00000003356)

###############################

  def test_calc_ionization_probs_with_empty_list(self):
    f = peptide_flyability()
    self.assertEqual([], f.ionization_probs([]))

  def test_calc_ionization_probs_with_non_empty_list(self):
    f = peptide_flyability()
    test = f.ionization_probs(['MASQASEK', 'DISLVQTPHK', 'VEVNEK'])

    self.assertEqual(len(test), 3)
    self.assertEqual(test[0]['seq'], 'MASQASEK')
    self.assertAlmostEqual(test[0]['prob'], 0.798247804562)
    self.assertAlmostEqual(test[0]['svm_score'], 0.209218277881)

    self.assertEqual(test[1]['seq'], 'DISLVQTPHK')
    self.assertAlmostEqual(test[1]['prob'], 0.985826563481)
    self.assertAlmostEqual(test[1]['svm_score'], 0.698591590494)

    self.assertEqual(test[2]['seq'], 'VEVNEK')
    self.assertAlmostEqual(test[2]['prob'],0.874711212861)
    self.assertAlmostEqual(test[2]['svm_score'],0.310199536208)

###############################


if __name__ == '__main__':
    unittest.main()
