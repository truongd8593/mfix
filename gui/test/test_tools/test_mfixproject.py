import sys
import os
import unittest

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(myPath, os.pardir, os.pardir))

from gui.tools.mfixproject import Keyword, Equation, Project

class TestEquation(unittest.TestCase):
    """ Unit tests for equation """

    def test_mutiply(self):
        """ test multiplication """
        eq = Equation('2*5')

        self.assertEqual(float(eq), 10.0)

class TestKeyword(unittest.TestCase):
    """ Unit tests for keyword """

    def test_attr(self):
        """ test attributes """
        kw = Keyword('key', 'val', comment='comment', dtype=float, args=[1])

        self.assertEqual(kw.key, 'key')

    def test_int(self):
        """ test integer keyword """
        kw = Keyword('key', 10)

        assert int(kw) == 10

    def test_float(self):
        """ test float keyword """
        kw = Keyword('key', 1.123)

        self.assertEqual(float(kw), 1.123)

    def test_cmp_eq(self):
        """ test equal operator """
        kw = Keyword('key', -100.1)

        self.assertEqual(kw, -100.1)

    def test_cmp_lt(self):
        """ test less than operator """
        kw = Keyword('key', 100.2)

        self.assertLess(kw, 1000)

    def test_cmp_gt(self):
        """ test greater than operator """
        kw = Keyword('key', 100.3)

        self.assertGreater(kw, -100)

    def test_if(self):
        """ test if """
        kw = Keyword('key', True)

        if kw:
            assert True

    def test_equation(self):
        """ test equation detection """
        eq = Equation('2*5')
        kw = Keyword('key', eq)

        self.assertEqual(10.0, float(eq))
        self.assertEqual(10.0, float(kw))


class TestProject(unittest.TestCase):
    """ test the project """
    def test_parsemfixdat_keyvalue(self):
        """ parse a string """

        testString = """
                     key = 'value'
                     """

        project = Project(testString)

        self.assertEqual('value',
            project._keyword_dict['key'].value)

    def test_parsemfixdat_booleans(self):
        """ parse booleans """

        testString = """
                     truea = .True.
                     trueb = .T.
                     truec = .true.
                     trued = .t.
                     falsea = .False.
                     falseb = .F.
                     falsec = .false.
                     falsed = .f.
                     """

        project = Project(testString)

        expected = {'truea':  True,
                    'trueb':  True,
                    'truec':  True,
                    'trued':  True,
                    'falsea': False,
                    'falseb': False,
                    'falsec': False,
                    'falsed': False,
                    }

        self.assertDictEqual(expected, project._keyword_dict)

    def test_parsemfixdat_floats(self):
        """ parse floats """

        testString = """
                     keya = 3.0
                     keyb = 3D0
                     keyc = 3.
                     keyd = 3.E5
                     """

        project = Project(testString)

        expected = {'keya': 3.0,
                    'keyb': 3.0,
                    'keyc': 3.0,
                    'keyd': 3.0E5,
                    }

        for key in expected.keys():
            self.assertEqual(expected[key],
                             float(project._keyword_dict[key]))

    def test_parsemfixdat_exponentialnotation(self):
        """ parse exponentials """

        testString = """
                     keya = 3.0E-10
                     keyb = 1.1d10
                     """

        project = Project(testString)

        expected = {'keya': 3.0E-10,
                    'keyb': 1.1E10,
                    }

        self.assertDictEqual(expected, project._keyword_dict)

    #TODO: Fix
    @unittest.skip("Fix")
    def test_parsemfixdat_equations(self):
        """ parse equations """

        testString = """
                     ncy = 5 2*10 @(2*3)
                     """

        project = Project(testString)

        self.assertEqual('@(2*3)', str(project.variablegrid[4]['ncy']))

    def test_parsemfixdat_comments(self):
        """ ignore comments """

        testString = """
                     # comment before
                     keya = 'value' ! in line comment
                     keyb = 'value'# in line comment
                     ! comment after
                     """

        project = Project(testString)

        expected = {'keya': 'value',
                    'keyb': 'value'
                    }
        self.assertDictEqual(expected, project._keyword_dict)

    def test_parsemfixdat_commentedKeywords(self):
        """ ignore commented keywords """
        testString = """
                     ! key = 3.0
                     # key = 'test'
                     """

        project = Project(testString)

        self.assertDictEqual({}, project._keyword_dict)

    #TODO: Fix
    @unittest.skip("Fix")
    def test_parsemfixdat_expandshorthand(self):
        """ expand shorthand """
        testString = """
                     ic_ep_g = .4      1.0
                     leq_sweep = 3*'ISIS' 'JSJS'
                     ncx = 4*10
                     ncy = 5 2*10 @(2*3)
                     """

        newProject = Project(testString)

        # Check IC
        assert newProject.ics[1]['ic_ep_g'] == 0.4
        assert newProject.ics[2]['ic_ep_g'] == 1.0

        # Check linear EQ
        for i in range(1, 4):
            assert newProject.linearEq[i]['leq_sweep'] == 'ISIS'
        assert newProject.linearEq[4]['leq_sweep'] == 'JSJS'

        # Check Grid
        for i in range(1, 5):
            assert newProject.variablegrid[i]['ncx'] == 10

        assert newProject.variablegrid[1]['ncy'] == 5
        assert newProject.variablegrid[2]['ncy'] == 10
        assert newProject.variablegrid[3]['ncy'] == 10
        assert str(newProject.variablegrid[4]['ncy']) == '@(2*3)'

    def test_parsemfixdat_argswithmultiplevalues(self):
        """ parse args with multiple values"""
        testString = """
                     ic_ep_g(2) = .4      1.0
                     ic_rop_s(2,1) = 10 30
                     """

        newProject = Project(testString)

        # Check IC
        assert newProject.ics[2]['ic_ep_g'] == 0.4
        assert newProject.ics[3]['ic_ep_g'] == 1.0
        assert newProject.ics[2].solids[1]['ic_rop_s'] == 10
        assert newProject.ics[3].solids[1]['ic_rop_s'] == 30

    def test_parsemfixdat_ic(self):
        """ parse ics """
        testString = """
                     ic_x_w(1) = 2.0
                     ic_x_g(1,2) = 0.3
                     ic_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        assert newProject.ics[1].gasSpecies[2]['ic_x_g'] == 0.3
        assert newProject.ics[1]['ic_x_w'] == 2.0
        assert newProject.ics[1].solids[1]['ic_t_s'] == 300.0

    def test_parsemfixdat_bc(self):
        """ parse bcs """
        testString = """
                     bc_x_w(1) = 2.0
                     bc_x_g(1,2) = 0.3
                     bc_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        assert newProject.bcs[1].gasSpecies[2]['bc_x_g'] == 0.3
        assert newProject.bcs[1]['bc_x_w'] == 2.0
        assert newProject.bcs[1].solids[1]['bc_t_s'] == 300.0

    def test_parsemfixdat_ps(self):
        """ parse point sources """
        testString = """
                     ps_x_w(1) = 2.0
                     ps_x_g(1,2) = 0.3
                     ps_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        assert newProject.pss[1].gasSpecies[2]['ps_x_g'] == 0.3
        assert newProject.pss[1]['ps_x_w'] == 2.0
        assert newProject.pss[1].solids[1]['ps_t_s'] == 300.0
