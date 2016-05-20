import sys
import os
import unittest
import math

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(myPath, os.pardir, os.pardir))

from project import Keyword, Equation, Project

class TestEquation(unittest.TestCase):

    def test_add(self):
        eq = Equation('2+5')
        self.assertEqual(float(eq), 7.0)

    def test_add_eq(self):
        eq = Equation('2+5')
        self.assertEqual(float(eq)+10, 17.0)

    def test_sub(self):
        eq = Equation('2-5')
        self.assertEqual(float(eq), -3.0)

    def test_sub_eq(self):
        eq = Equation('2-5')
        self.assertEqual(float(eq) - 10, -13.0)

    def test_mul(self):
        eq = Equation('2*5')
        self.assertEqual(float(eq), 10.0)

    def test_mul_eq(self):
        eq = Equation('2*5')
        self.assertEqual(float(eq) * 2, 20.0)

    def test_div(self):
        eq = Equation('10/2')
        self.assertEqual(float(eq), 5.0)

    def test_div_eq(self):
        eq = Equation('10/2')
        self.assertAlmostEqual(float(eq) / 2.0, 5.0 / 2.0)

    def test_pow(self):
        eq = Equation('10**2')
        self.assertEqual(float(eq), 100)

    def test_variable(self):
        eq = Equation('10*pi / 15')
        self.assertAlmostEqual(float(eq), 10 * math.pi/15)


class TestKeyword(unittest.TestCase):
    def test_attr(self):
        kw = Keyword('key', 'val', comment='comment', dtype=float, args=[1])
        self.assertEqual(kw.key, 'key')

    def test_int(self):
        kw = Keyword('key', 10)
        assert int(kw) == 10

    def test_float(self):
        kw = Keyword('key', 1.123)
        self.assertEqual(float(kw), 1.123)

    def test_cmp_eq(self):
        kw = Keyword('key', -100.1)
        self.assertEqual(kw, -100.1)

    def test_cmp_lt(self):
        kw = Keyword('key', 100.2)
        self.assertLess(kw, 1000)

    def test_cmp_gt(self):
        kw = Keyword('key', 100.3)
        self.assertGreater(kw, -100)

    def test_if(self):
        kw = Keyword('key', True)
        if kw:
            assert True

    def test_eq(self):
        eq = Equation('2*5')
        kw = Keyword('key', eq)

        self.assertEqual(10.0, float(eq))
        self.assertEqual(10.0, float(kw))

    def test_update_value(self):
        kw = Keyword('key', True)

        for dtype, value in ((bool, True), (str, 'test'), (float, 10.0),
                             (int, 1), (float, None)):
            kw.dtype = dtype
            kw.updateValue(value)

            self.assertEqual(kw.value, value)

    def test_update_eq(self):
        kw = Keyword('key', Equation('2*4'))

        self.assertIsInstance(kw.value, Equation)

        kw.updateValue('10*2')
        self.assertIsInstance(kw.value, Equation)
        self.assertEqual(20.0, float(kw))

    def test_update_eq(self):
        kw = Keyword('key', 2.4)

        self.assertEqual(float, kw.dtype)

        kw.updateValue('10*20')

        self.assertIsInstance(kw.value, Equation)
        self.assertEqual(10.0*20, float(kw))

    def test_lower(self):
        kw = Keyword('key', 'CamelCase')
        self.assertEqual('camelcase', kw.lower())


class TestProject(unittest.TestCase):
    """ test the project """

    def setUp(self):
        self.project = Project()

    def test_parseKeywordLine_str(self):
        results = list(self.project.parseKeywordLine(
            "key = 'value'"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual('value', result[2])

    def test_parseKeywordLine_int(self):
        results = list(self.project.parseKeywordLine(
            "key = 1"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(1, result[2])

    def test_parseKeywordLine_float(self):
        results = list(self.project.parseKeywordLine(
            "key = 10.0"
            ))

        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(10.0, result[2])

    def test_parseKeywordLine_float_d(self):
        results = list(self.project.parseKeywordLine(
            "key = 10.0 10d 10.d 10D 10.D 10.0D" ))

        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10.0, result[2])

    def test_parseKeywordLine_float_e(self):
        results = list(self.project.parseKeywordLine(
            "key = 10.0 10e 10.e 10E 10.E 10.0E" ))

        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10.0, result[2])

    def test_parseKeywordLineShorthandEq(self):
        results = list(self.project.parseKeywordLine(
            "key = @(2*5) 2*10 10 @(4+6) 10" ))
        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10.0, float(result[2]))

    def test_parseKeywordLine_exp(self):
        results = list(self.project.parseKeywordLine(
            "key = 1.0e10 1.d10 1.0e+10 1.d+10"))

        self.assertEqual(len(results), 4)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(1.0E10, result[2])

    def test_parseKeywordLine_exp_neg(self):
        results = list(self.project.parseKeywordLine(
            "key = 1.0e-10 1.d-10"))

        self.assertEqual(len(results), 2)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(1.0E-10, result[2])

    def test_parseKeywordLine_bool(self):
        results = list(self.project.parseKeywordLine(
            "key = .T."
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(True, result[2])

    def test_parseKeywordLine_arg(self):
        results = list(self.project.parseKeywordLine(
            "key(3) = .T."
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([3], result[1])
        self.assertEqual(True, result[2])

    def test_parseKeywordLine_mul_arg(self):
        results = list(self.project.parseKeywordLine(
            "key(3,2) = .T."
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([3, 2], result[1])
        self.assertEqual(True, result[2])

    def test_parseKeywordLine_eq_multi(self):
        results = list(self.project.parseKeywordLine(
            "key = @(2*3)"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(2 * 3.0, float(result[2]))


    def test_parseKeywordLine_eq_mul_wspace(self):
        results = list(self.project.parseKeywordLine(
            "key = @( 2*3)"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(2 * 3.0, float(result[2]))

    def test_parseKeywordLine_eq_divide(self):
        results = list(self.project.parseKeywordLine(
            "key = @( 2/ 3)"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(2 / 3.0, float(result[2]))

    def test_parseKeywordLine_eq_pi(self):
        results = list(self.project.parseKeywordLine(
            "key = @( 2*pi)"
            ))
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual('key', result[0])
        self.assertEqual([], result[1])
        self.assertEqual(2 * math.pi, float(result[2]))

    def test_parseKeywordLine_shorthand_str(self):
        results = list(self.project.parseKeywordLine(
            "key = 4*'ISIS'"
            ))
        self.assertEqual(len(results), 4)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual('ISIS', result[2])

    def test_parseKeywordLine_shorthand_float(self):
        results = list(self.project.parseKeywordLine(
            "key = 4*6.7"
            ))
        self.assertEqual(len(results), 4)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(6.7, result[2])

    def test_parseKeywordLine_shorthand_float_other_val(self):
        results = list(self.project.parseKeywordLine(
            "key = 4*6.7 6.7 6.7"
            ))
        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(6.7, result[2])

    def test_parseKeywordLine_shorthand_eq(self):
        results = list(self.project.parseKeywordLine(
            "key = 4*6.7 @(1*6.7)"
            ))

        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(6.7, float(result[2]))

    def test_parseKeywordLine_shorthand_eq_2(self):
        results = list(self.project.parseKeywordLine(
            "key =  @(2*5) 10 2*10 @(2*5) @(2*5)"
            ))
        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10, float(result[2]))

    def test_parseKeywordLine_shorthand_eq_3(self):
        results = list(self.project.parseKeywordLine(
            "key =  @(5+ 5) 10 2*10 @(5+5) @ (5 + 5) "
            ))
        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10, float(result[2]))


    def test_parseKeywordLine_shorthand_bad(self):
        self.skipTest("what does 10*10*10 expand to?")
        results = list(self.project.parseKeywordLine(
            "key =  10*10*10"
            ))
        self.assertEqual(len(results), 6)
        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(10, float(result[2]))

    def test_parseKeywordLine_shorthand_eq_exp(self):
        results = list(self.project.parseKeywordLine(
            "key = 4*6.7 @(1*6.7)"
            ))

        for i, result in enumerate(results):
            self.assertEqual('key', result[0])
            self.assertEqual([i + 1], result[1])
            self.assertEqual(6.7, float(result[2]))

    def test_parseKeywordLine_twokeys(self):
        results = list(self.project.parseKeywordLine(
            "key0 = 10 key1 = 10"
            ))

        for i, result in enumerate(results):
            self.assertEqual('key' + str(i), result[0])
            self.assertEqual([], result[1])
            self.assertEqual(10, result[2])

    def test_parsemfixdat_keyvalue(self):
        testString = """
                     key = 'value'
                     """

        project = Project(testString)

        self.assertEqual('value',
            project._keyword_dict['key'].value)

    def test_parsemfixdat_booleans(self):
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

        expected = {'truea': True,
                    'trueb': True,
                    'truec': True,
                    'trued': True,
                    'falsea': False,
                    'falseb': False,
                    'falsec': False,
                    'falsed': False,
                    }

        self.assertDictEqual(expected, project._keyword_dict)

    def test_parsemfixdat_floats(self):
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
        testString = """
                     keya = 3.0E-10
                     keyb = 1.1d10
                     """

        project = Project(testString)

        expected = {'keya': 3.0E-10,
                    'keyb': 1.1E10,
                    }

        self.assertDictEqual(expected, project._keyword_dict)

    def test_parsemfixdat_equations(self):
        testString = """
                     ncy = 5 2*10 @(2*3)
                     """

        project = Project(testString)

        self.assertEqual('@(2*3)', str(project.variablegrid[4]['ncy']))

    def test_parsemfixdat_comments(self):
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
        testString = """
                     ! key = 3.0
                     # key = 'test'
                     """

        project = Project(testString)

        self.assertDictEqual({}, project._keyword_dict)

    def test_parsemfixdat_expandshorthand(self):
        testString = """
                     ic_ep_g = .4      1.0
                     leq_sweep = 3*'ISIS' 'JSJS'
                     ncx = 4*10
                     ncy = 5 2*10 @(2*3)
                     """

        newProject = Project(testString)

        # Check IC
        self.assertEqual(0.4, newProject.ics[1]['ic_ep_g'])
        self.assertEqual(1.0, newProject.ics[2]['ic_ep_g'])

        # Check linear EQ
        for i in range(1, 4):
            self.assertEqual('ISIS', newProject.linearEq[i]['leq_sweep'])
        self.assertEqual('JSJS', newProject.linearEq[4]['leq_sweep'])

        # Check Grid
        for i in range(1, 5):
            self.assertEqual(10, newProject.variablegrid[i]['ncx'])

        self.assertEqual(5, newProject.variablegrid[1]['ncy'])
        self.assertEqual(10, newProject.variablegrid[2]['ncy'])
        self.assertEqual(10, newProject.variablegrid[3]['ncy'])
        self.assertEqual('@(2*3)', str(newProject.variablegrid[4]['ncy']))

    def test_parsemfixdat_argswithmultiplevalues(self):
        testString = """
                     ic_ep_g(2) = .4      1.0
                     ic_rop_s(2,1) = 10 30
                     """

        newProject = Project(testString)

        # Check IC
        self.assertEqual(0.4, newProject.ics[2]['ic_ep_g'])
        self.assertEqual(1.0, newProject.ics[3]['ic_ep_g'])
        self.assertEqual(10, newProject.ics[2].solids[1]['ic_rop_s'])
        self.assertEqual(30, newProject.ics[3].solids[1]['ic_rop_s'])

    def test_parsemfixdat_ic(self):
        testString = """
                     ic_x_w(1) = 2.0
                     ic_x_g(1,2) = 0.3
                     ic_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        self.assertEqual(0.3, newProject.ics[1].gasSpecies[2]['ic_x_g'])
        self.assertEqual(2.0, newProject.ics[1]['ic_x_w'])
        self.assertEqual(300.0, newProject.ics[1].solids[1]['ic_t_s'])

    def test_parsemfixdat_bc(self):
        testString = """
                     bc_x_w(1) = 2.0
                     bc_x_g(1,2) = 0.3
                     bc_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        self.assertEqual(0.3, newProject.bcs[1].gasSpecies[2]['bc_x_g'])
        self.assertEqual(2.0, newProject.bcs[1]['bc_x_w'])
        self.assertEqual(300.0, newProject.bcs[1].solids[1]['bc_t_s'])

    def test_parsemfixdat_ps(self):
        testString = """
                     ps_x_w(1) = 2.0
                     ps_x_g(1,2) = 0.3
                     ps_t_s(1,1) = 300.0
                     """

        newProject = Project(testString)

        self.assertEqual(0.3, newProject.pss[1].gasSpecies[2]['ps_x_g'])
        self.assertEqual(2.0, newProject.pss[1]['ps_x_w'])
        self.assertEqual(300.0, newProject.pss[1].solids[1]['ps_t_s'])
