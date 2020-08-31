import unittest

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from CFML_api import SpaceGroup

import numpy
import unittest

class TestSymmetryGroups(unittest.TestCase):
    '''
    Unittest for the geometry-related functions
    '''

    def test_print_description(self):
        self.group = SpaceGroup(1)
        self.group.print_description()

    def test_getters(self):
        self.group = SpaceGroup(1)
        numpy.testing.assert_almost_equal(self.group.lattice_translation, numpy.array([[0.], [0.], [0.]]))
        self.assertFalse(self.group.is_hexa)
        self.assertEqual(self.group.space_group_number, 1)
        self.assertEqual(self.group.space_group_symbol, "P 1")
        self.assertEqual(self.group.hall_symbol, "P 1")
        self.assertEqual(self.group.generalized_hall_symbol, "")
        #self.group.crystal_system
        #self.group.laue_class
        #self.group.point_group
        #self.group.extra_information
        #self.group.space_group_setting
        #self.group.space_group_lattice_type
        #self.group.space_group_lattice_type_symbol
        #self.group.number_of_lattice_points
        #self.group.bravais
        #self.group.centre
        #self.group.centric
        #self.group.inversion_centre
        #self.group.number_of_symops
        #self.group.multiplicity
        #self.group.operators_minimum_number
        #self.group.symmetry_operators
        #self.group.symmetry_operators_as_text
        #self.group.wyckoff_info
        #self.group.asymetric_unit
        #self.group.symmetry_operators[0].rotation_matrix
        #self.group.symmetry_operators[0].translation_matrix
        #self.group.wyckoff_info.num_orbit
        #self.group.wyckoff_info.orbits[0].multiplicity
        #self.group.wyckoff_info.orbits[0].site
        #self.group.wyckoff_info.orbits[0].orbit_number
        #self.group.wyckoff_info.orbits[0].str_orig
        #self.group.wyckoff_info.orbits[0].str_orbit[0]     
            
def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TestSymmetryGroups))
    return s

if __name__ == '__main__':
    unittest.main(verbosity=2)
            
        
