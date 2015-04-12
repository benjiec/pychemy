import unittest, tempfile, os
from pychemy.chem_structure import chem_structure as cs
from pychemy.chem_structure import chem_graph as cg
from pychemy.elements import ELEMENTS


inchi1 = 'InChI=1/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H,(H4,6,7,8,9,10,11)/f/h8,10H,6H2'

# Aspirine
inchi2 = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'

class chem_structure_Testing(unittest.TestCase):


  ###############################
  # chem_structure

  def test_chem_structure_without_input(self):
    CS = cs()
    self.assertEqual(CS.inchi, None)

  def test_chem_structure_with_inchi_input(self):
    CS = cs(inchi = inchi1)
    self.assertEqual(CS.inchi, inchi1)
    self.assertEqual(CS.mol.NumAtoms(), 16)
    self.assertEqual(CS.mol.NumBonds(), 17)

  ###############################
  # chem_graph

  def test_graph_from_OBMol(self):
    CS = cs(inchi = inchi1)
    G = cg.graph_from_OBMol(CS.mol)

    self.assertEqual(G.number_of_nodes(), 16)

    self.assertEqual(G.node[1]['atom'], ELEMENTS['C'])
    self.assertEqual(G.node[2]['atom'], ELEMENTS['C'])
    self.assertEqual(G.node[3]['atom'], ELEMENTS['C'])
    self.assertEqual(G.node[4]['atom'], ELEMENTS['C'])
    self.assertEqual(G.node[5]['atom'], ELEMENTS['C'])
    self.assertEqual(G.node[6]['atom'], ELEMENTS['N'])
    self.assertEqual(G.node[7]['atom'], ELEMENTS['N'])
    self.assertEqual(G.node[8]['atom'], ELEMENTS['N'])
    self.assertEqual(G.node[9]['atom'], ELEMENTS['N'])
    self.assertEqual(G.node[10]['atom'], ELEMENTS['N'])
    self.assertEqual(G.node[11]['atom'], ELEMENTS['O'])
    self.assertEqual(G.node[12]['atom'], ELEMENTS['H'])
    self.assertEqual(G.node[13]['atom'], ELEMENTS['H'])
    self.assertEqual(G.node[14]['atom'], ELEMENTS['H'])
    self.assertEqual(G.node[15]['atom'], ELEMENTS['H'])
    self.assertEqual(G.node[16]['atom'], ELEMENTS['H'])

    for n in range(1,17):
      self.assertTrue('X' in G.node[n])
      self.assertTrue('Y' in G.node[n])

    for e in G.edges():
      self.assertTrue('order' in G[e[0]][e[1]])

    self.assertEqual(G.number_of_edges(), 17)
    edge_set = [(1, 12), (12,  1),
                (1,  7),  (7,  1),
                (1,  8),  (8,  1),
                (2,  7),  (7,  2),
                (3,  8),  (8,  3),
                (2,  3),  (3,  2),
                (2,  4),  (4,  2),
                (4, 11), (11,  4),
                (4, 10), (10,  4),
               (10, 16), (16, 10),
                (5, 10), (10,  5),
                (5,  6),  (6,  5),
                (5,  9),  (9,  5),
                (3,  9),  (9,  3),
                (8, 15), (15,  8),
                (6, 13), (13,  6),
                (6, 14), (14,  6)]

    for e in G.edges():
      self.assertIn(e,edge_set)
      edge_set.remove(e)
      edge_set.remove((e[1],e[0]))

    self.assertEqual(edge_set, [])

  def test_chem_graph_with_empty_input(self):
    CG = cg()
    self.assertEqual(CG.G.number_of_nodes(), 0)

  def test_chem_graph_from_OBMol(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    self.assertEqual(CG.G.number_of_nodes(), 16)
    self.assertEqual(CG.G.number_of_edges(), 17)

  def test_chem_formula(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    self.assertEqual(CG.chem_formula(),CS.chem_formula())

  def test_mass(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    self.assertTrue(abs(CG.mass() - CS.mass()) / CG.mass() < 0.001)

  def test_gen_frag_with_no_additional_steps(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    self.assertEqual(CG.gen_frag(0), [CG])

  def test_gen_frag_has_no_empty_molecules(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    for level in xrange(0, 4):
      frag = CG.gen_frag(level)
      nonempty_frag = filter(lambda f: len(f.chem_formula()) > 0, frag)
      self.assertEqual(frag, nonempty_frag)


  def test_gen_frag_with_additional_steps_uniqueness(self):
    CS = cs(inchi = inchi1)
    CG = cg(CS.mol)
    for level in range(1, 4):
      frag = CG.gen_frag(level)
      self.assertEqual(len(frag), len(set(frag)))
