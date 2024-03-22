import unittest
from goph420_lab02 import rootFunctions as rgr


class TestSecantModified(unittest.TestCase):

    def test_Secant_Modified(self):
        x0 = 5
        dx = 0.1
        fx = rgr.function_Test()
        result = rgr.root_secant_modified(x0, dx, fx)
        self.assertAlmostEqual(result, -1, 5)

    def test_Secant_Modified_raises(self):

        x0 = 'A'
        dx = 1
        fx = rgr.function_Test()
        self.assertRaises(ValueError, rgr.root_secant_modified, x0, dx, fx)

        x0 = 10
        dx = 'A'
        fx = rgr.function_Test()
        self.assertRaises(ValueError, rgr.root_secant_modified, x0, dx, fx)

        x0 = 10
        dx = 1
        fx = 1
        self.assertRaises(TypeError, rgr.root_secant_modified, x0, dx, fx)


if __name__ == '__main__':
    unittest.main()
