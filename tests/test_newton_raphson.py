import unittest
from goph420_lab02 import integration as Itg


class TestNewtonRaphson(unittest.TestCase):

    def test_Newton_Raphson(self):
        x0 = 10
        fx = Itg.function_Test()
        dfdx = Itg.function_dfdx_Test()
        result = Itg.root_newton_raphson(x0, fx, dfdx)
        self.assertAlmostEqual(result, -1, 5)

    def test_Newton_Raphson_raises(self):

        x0 = 'A'
        fx = Itg.function_Test()
        dfdx = Itg.function_dfdx_Test()
        self.assertRaises(ValueError, Itg.root_newton_raphson, x0, fx, dfdx)

        x0 = 10
        fx = 1
        dfdx = Itg.function_dfdx_Test()
        self.assertRaises(TypeError, Itg.root_newton_raphson, x0, fx, dfdx)

        x0 = 10
        fx = Itg.function_Test()
        dfdx = 1
        self.assertRaises(TypeError, Itg.root_newton_raphson, x0, fx, dfdx)


if __name__ == '__main__':
    unittest.main()
