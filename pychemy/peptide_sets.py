# Tools for grouping proteins into non-conflicting sets for
# multiplexed screening of homologues

import re

from pychemy.peptides import mass_from_sequence
import numpy as np
import csv
import os

AA = re.compile('[ACDEFGHILMNPQSTVWY]+[RK]')
AA_check = re.compile('[ACDEFGHILKMNPQRSTVWY]')


def check_AA(test_string):
  """Check that all entries in sequence represent valid amino acids"""
  if len(AA_check.findall(test_string)) == len(test_string):
    return True
  else:
    raise Exception("Invalid amino acid in seq: " + test_string)


def check_proteins(proteins):
  """
  Check that all proteins in list contain only valid amino acids
  input:
       proteins: list of protein sequences as strings ['protein_seq', ...]
  """
  for p in proteins:
    check_AA(p)
  return True


def get_peptides(protein_seq):
  """
  Creates list of peptides from protein sequence
  input:
         protein_seq: string containing protein sequence of amino acids
  output:
         list of strings containing peptide sequences from protein sequence
         ['peptide_seq', ...]
  """
  check_AA(protein_seq)
  return AA.findall(protein_seq)


def all_peptides(proteins):
  """
  Creates list of all peptides from protein set
  inputs:
          proteins: list of protein sequences as strings ['protein_seq', ...]
  output:
          list of strings containing peptide sequences from all protein
          sequences ['peptide_seq', ...]
  """
  return [p for peptides in [get_peptides(ps)
          for ps in proteins] for p in peptides]


def distinct_peptides(proteins):
  """
  Creates list of all distinct peptides from protein set
  inputs:
          proteins: list of protein sequences as strings ['protein_seq', ...]
  output:
          list of strings containing distinct peptide sequences from all
          protein sequences ['peptide_seq', ...]
  """
  return list(set(all_peptides(proteins)))


def unique_identifiers(proteins):
  """
  Creates list of unique identifiers for each protein from set
  inputs:
          proteins: list of protein sequences as strings ['protein_seq', ...]
  output:
          list of lists containing protein sequence and a unique identifying
          peptide sequence [['protein_seq', 'peptide_seq'], ...]
  """
  all_pep = [[ps[0], pep] for ps in enumerate(proteins)
             for pep in get_peptides(ps[1])]
  dp = distinct_peptides(proteins)
  ui = []
  for idx, seq in enumerate(dp):
    idx = [item[0] for item in all_pep if item[1] == seq]
    if len(idx) == 1:
      ui.append([proteins[idx[0]], seq])
  return ui


def peptides_per_protein(ui):
  """
  Creates list of unique identifier peptides grouped by protein
  input:
          ui: a set of unique identifiers such as that produced by
          unique_identifiers([proteins])
          [['protein_seq', 'peptide_seq'],.....]
  output:
          list of unique identifier peptides grouped by their protein
          [['protein_seq', ['peptide_seq', 'peptide_seq', ...]],.....]
  """
  proteins = list(set([item[0] for item in ui]))
  return [[p, [item[1] for item in ui if item[0] == p]] for p in proteins]


def min_ui_count(proteins):
  """
  Counts the minimum number of unique identifier peptides across all proteins
  in a set
  input:
         proteins: list of protein sequences as strings ['protein_seq', ...]
  output:
         minimum number of unique identifier peptides across all proteins in
         a set
  """
  temp = []
  for p in peptides_per_protein(unique_identifiers(proteins)):
    temp.append(len(p[1]))
  if len(proteins) > len(temp):
    return 0
  else:
    return min(temp)


def balanced_sets(proteins, max_set_size=10, unique=1):
  """
  Creates balanced sets of proteins based on unique peptides subject to
  contraints
  inputs:
          proteins:     list of protein sequences as strings
                        ['protein_seq', ...]
          max_set_size: maximum number of proteins permitted in a set
          unique:       minimum number of unique peptides required per protein
  output:
          sets of protein sequences that meet the set size and uniqueness
          requirements [['protein_seq', ...], ...]
  """
  check_proteins(proteins)
  sets = []
  for p in proteins:
    placed = False
    for idx, s in enumerate(sets):
      if len(s) < max_set_size and not placed:
        if min_ui_count(s + [p]) >= unique:
          placed = True
          sets[idx] += [p]
    if not placed:
     sets.append([p])
  return sets


##############################################################################
# Calculate peptide ionizability based on the method used by STEPP from PNNL #
##############################################################################

class peptide_flyability():
  def __init__(self):
    self.NUM_FEATURES = 35

    self.AA_ORDER = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(dir, "resources/STEPP/STEPP_NormFactor_mean.txt")) as tsv:
      self.STEPP_mean = np.array([float(line[0])
                                 for line in csv.reader(tsv, delimiter="\t")])

    with open(os.path.join(dir, "resources/STEPP/STEPP_NormFactor_std.txt")) as tsv:
      self.STEPP_std = np.array([float(line[0])
                                for line in csv.reader(tsv, delimiter="\t")])

    with open(os.path.join(dir, "resources/STEPP/STEPP_Weights.txt")) as tsv:
      self.STEPP_weights = np.array([float(line[0])
                                    for line in csv.reader(tsv,
                                                           delimiter="\t")])

    self.non_polar_hydrophobic = set(['A', 'F', 'G', 'I', 'L',
                                      'M', 'P', 'V', 'W', 'Y'])
    self.polar_hydrophillic = set(['C', 'D', 'E', 'H', 'K',
                                   'N', 'Q', 'R', 'S', 'T'])
    self.uncharged_polar_hydrophillic = set(['C', 'N', 'Q', 'S', 'T'])
    self.charged_polar_hydrophillic = set(['D', 'E', 'H', 'K', 'R'])
    self.postive_polar_hydrophillic = set(['R', 'H', 'K'])
    self.negative_polar_hydrophillic = set(['D', 'E'])
    self.eisenberg_hydrophobicity = {'A': 0.620, 'C': 0.290, 'D': -0.900,
                                     'E': -0.740, 'F': 1.190, 'G': 0.480,
                                     'H': -0.400, 'I': 1.380, 'K': -1.500,
                                     'L': 1.060, 'M': 0.640, 'N': -0.780,
                                     'P': 0.120, 'Q': -0.850, 'R': -2.530,
                                     'S': -0.180, 'T': -0.050, 'V': 1.080,
                                     'W': 0.810, 'Y': 0.260}
    self.hopp_woods_hydrophobicity = {'A': -0.500, 'C': -1.000, 'D': 3.000,
                                      'E':  3.000, 'F': -2.500, 'G': 0.000,
                                      'H': -0.500, 'I': -1.800, 'K': 3.000,
                                      'L': -1.800, 'M': -1.300, 'N':  0.200,
                                      'P': 0.000, 'Q': 0.200, 'R': 3.000,
                                      'S': 0.300, 'T': -0.400, 'V': -1.500,
                                      'W': -3.400, 'Y': -2.300}
    self.kyte_doolittle_hydrophobicity = {'A': 1.800, 'C': 2.500,
                                          'D': -3.500, 'E': -3.500,
                                          'F': 2.800, 'G': -0.400,
                                          'H': -3.200, 'I': 4.500,
                                          'K': -3.900, 'L': 3.800,
                                          'M': 1.900, 'N': -3.500,
                                          'P': -1.600, 'Q': -3.500,
                                          'R': -4.500, 'S': -0.800,
                                          'T': -0.700, 'V': 4.200,
                                          'W': -0.900, 'Y': -1.300}
    self.roseman_hydropathicity = {'A': 0.390, 'C': 0.250, 'D': -3.810,
                                   'E': -2.910, 'F': 2.270, 'G': 0.000,
                                   'H': -0.640, 'I': 1.820, 'K': -2.770,
                                   'L': 1.820, 'M': 0.960, 'N': -1.910,
                                   'P': 0.990, 'Q': -1.300, 'R': -3.950,
                                   'S': -1.240, 'T': -1.000, 'V': 1.300,
                                   'W': 2.130, 'Y': 1.470}
    self.grantham_polarity = {'A': 8.100, 'C': 5.500, 'D': 13.000,
                              'E': 12.300, 'F': 5.200, 'G': 9.000,
                              'H': 10.400, 'I':  5.200, 'K': 11.300,
                              'L': 4.900, 'M': 5.700, 'N': 11.600,
                              'P': 8.000, 'Q': 10.500, 'R': 10.500,
                              'S':  9.200, 'T': 8.600, 'V': 5.900,
                              'W': 5.400, 'Y': 6.200}
    self.zimmerman_polarity = {'A': 0.000, 'C': 1.480, 'D': 49.700,
                               'E': 49.900, 'F': 0.350, 'G': 0.000,
                               'H': 51.600, 'I':  0.130, 'K': 49.500,
                               'L': 0.130, 'M': 1.430, 'N': 3.380,
                               'P': 1.580, 'Q': 3.530, 'R': 52.000,
                               'S': 1.670, 'T': 1.660, 'V': 0.130,
                               'W': 2.100, 'Y': 1.610}
    self.zimmerman_bulkiness = {'A': 11.500, 'C': 13.460, 'D': 11.680,
                                'E': 13.570, 'F': 19.800, 'G': 3.400,
                                'H': 13.690, 'I': 21.400, 'K': 15.710,
                                'L': 21.400, 'M': 16.250, 'N': 12.820,
                                'P': 17.430, 'Q': 14.450, 'R': 14.280,
                                'V': 21.570, 'S': 9.470, 'T': 15.770,
                                'W': 21.670, 'Y': 18.030}

    self.zeta_pos = -0.3821  # K_pos
    self.sigma_pos = 0.3831  # scale_pos
    self.mu_pos = 0.0739     # location_pos

    self.zeta_neg = -0.1945  # K_neg
    self.sigma_neg = 0.3860  # scale_neg
    self.mu_neg = -0.4283    # location_neg

  def SVM_props(self, seq=''):
    """
    Creates feature vector containing properties used in SVM model
    inputs:
            seq:   string containing amino acid sequence
    output:
            numpy.array feature vector containing 35 properties of amino acid
            sequence
    FEATURE VECTOR
    1   Length
    2   Molecular weight
    3   Number of non-polar hydrophobic residues
    4   Number of polar hydrophilic residues
    5   Number of uncharged polar hydrophilic residues
    6   Number of charged polar hydrophilic residues
    7   Number of positively charged polar hydrophilic residues
    8   Number of negatively charged polar hydrophilic residues
    9   Hydrophobicity-Eisenberg scale (Eisenberg et al., 1984)
    10  Hydrophilicity-Hopp-Woods scale (Hopp and Woods, 1981)
    11  Hydrophobicity-Kyte-Doolittle (Kyte and Doolittle, 1982)
    12  Hydropathicity-Roseman scale (Roseman, 1988)
    13  Polarity-Grantham scale (Grantham, 1974)
    14  Polarity-Zimmerman scale (Zimmerman et al., 1968)
    15  Bulkiness (Zimmerman et al., 1968)
    16-35   Amino acid singlet counts in order: ACDEFGHIKLMNPQRSTVWY
    """
    if seq:
      props = [float(len(seq)),
               mass_from_sequence(seq),
               sum([1.0 for aa in seq if aa in self.non_polar_hydrophobic]),
               sum([1.0 for aa in seq if aa in self.polar_hydrophillic]),
               sum([1.0 for aa in seq
                    if aa in self.uncharged_polar_hydrophillic]),
               sum([1.0 for aa in seq
                    if aa in self.charged_polar_hydrophillic]),
               sum([1.0 for aa in seq
                    if aa in self.postive_polar_hydrophillic]),
               sum([1.0 for aa in seq
                    if aa in self.negative_polar_hydrophillic]),
               float(sum(map(lambda x: self.eisenberg_hydrophobicity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.hopp_woods_hydrophobicity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.kyte_doolittle_hydrophobicity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.roseman_hydropathicity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.grantham_polarity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.zimmerman_polarity[x],
                             seq)))/float(len(seq)),
               float(sum(map(lambda x: self.zimmerman_bulkiness[x],
                             seq)))/float(len(seq))]
      props += [sum([1 for aa in seq if aa == AA]) for AA in self.AA_ORDER]
      return np.array(props)
    else:
      return np.array([])

  def SVM_score(self, fv=np.array([])):
    """
    Calculates SVM score of a feature vector
    inputs:
            fv = numpy.array feature vector containing 35 properties of amino
            acid sequence
    output:
            SVM score
    """
    if len(fv) == self.NUM_FEATURES:
      dir = os.path.dirname(os.path.abspath(__file__))
      with open(os.path.join(dir, "resources/STEPP/STEPP_SupportVectors.txt")) as tsv:
        sv_idx = 0
        svm_score = 0
        for line in csv.reader(tsv, delimiter="\t"):
          sv = np.array([float(i) for i in line])
          kxy = np.dot(fv, sv)
          kxx = np.dot(fv, fv)
          kyy = np.dot(sv, sv)

          kxy /= np.sqrt(kxx * kyy)
          kxy += 10
          kxy = kxy ** 2

          kxy *= self.STEPP_weights[sv_idx]

          svm_score += kxy
          sv_idx += 1

        return svm_score
    else:
      return None

  def ionization_prob(self, svm_score):
    """
    Calculates ionization probability based on SVM score
    inputs:
            svm_score = SVM score output of calc_svm_score()
    output:
            ionization probability
    """

    def gevcdf(x, mu, sigma, zeta):
      z = (x-mu)/sigma
      t = np.real((1 + z * zeta + 0j) ** (-1/zeta))
      return np.exp(-t)

    def positive_probability(x):
      val = gevcdf(x, self.mu_pos, self.sigma_pos, self.zeta_pos)
      if np.isnan(val) or np.isinf(val):
        return 1
      else:
        return val

    def negative_probability(x):
      val = gevcdf(x, self.mu_neg, self.sigma_neg, self.zeta_neg)
      if np.isnan(val) or np.isinf(val):
        return 1
      else:
        return 1 - val

    pos_prob = positive_probability(svm_score)
    neg_prob = negative_probability(svm_score)
    divisor = pos_prob + neg_prob
    prob = pos_prob / divisor
    return prob if prob > 0 else 0

  def ionization_probs(self, peptides=[]):
    """
    Calculates ionization probability and SVM score for all peptides in a set
    inputs:
            peptides = array of peptide sequence strings
            ['SAMPLE', 'SAMPLE',...]
    output:
            array of dictionaries
            [{'seq': string, 'prob': float, 'svm_score': float}, ...]
    """
    # Normalize vector based on baseline mean and stdev
    norm_props = [np.divide(self.SVM_props(p) - self.STEPP_mean,
                            self.STEPP_std)
                  for p in peptides]
    out = []

    for idx, norm_prop in enumerate(norm_props):
      svm_score = self.SVM_score(norm_prop)
      if svm_score:
        out += [{'seq': peptides[idx],
                 'prob': self.ionization_prob(svm_score),
                 'svm_score': svm_score}]
      else:
        out += [{'seq': peptides[idx],
                 'prob': 0,
                 'svm_score': None}]
    return out
