import unittest
from goph420_lab02 import rootFunctions as rt


class TestNewtonRaphson(unittest.TestCase):

    def test_Newton_Raphson(self):
        x0 = 10
        fx = rt.function_Test()
        dfdx = rt.function_dfdx_Test()
        result = rt.root_newton_raphson(x0, fx, dfdx)
        self.assertAlmostEqual(result, -1, 5)

        x0 = 10
        fx = rt.function_Test_2()
        dfdx = rt.function_dfdx_Test_2()
        result = rt.root_newton_raphson(x0, fx, dfdx)
        self.assertAlmostEqual(result, 1, 5)

        x0 = 150
        result = rt.root_newton_raphson(x0, fx, dfdx)
        self.assertAlmostEqual(result, 100, 5)

    def test_Newton_Raphson_raises(self):

        x0 = 'A'
        fx = rt.function_Test()
        dfdx = rt.function_dfdx_Test()
        self.assertRaises(ValueError, rt.root_newton_raphson, x0, fx, dfdx)

        x0 = 10
        fx = 1
        dfdx = rt.function_dfdx_Test()
        self.assertRaises(TypeError, rt.root_newton_raphson, x0, fx, dfdx)

        x0 = 10
        fx = rt.function_Test()
        dfdx = 1
        self.assertRaises(TypeError, rt.root_newton_raphson, x0, fx, dfdx)


if __name__ == '__main__':
    unittest.main()
