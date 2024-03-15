import unittest
from goph420_lab02 import integration as Itg


class TestSecantModified(unittest.TestCase):

    def test_Secant_Modified(self):
        x0 = 5
        dx = 0.1
        fx = Itg.function_Test()
        result = Itg.root_secant_modified(x0, dx, fx)
        self.assertAlmostEqual(result, -1, 5)

    def test_Secant_Modified_raises(self):

        x0 = 'A'
        dx = 1
        fx = Itg.function_Test()
        self.assertRaises(ValueError, Itg.root_secant_modified, x0, dx, fx)

        x0 = 10
        dx = 'A'
        fx = Itg.function_Test()
        self.assertRaises(ValueError, Itg.root_secant_modified, x0, dx, fx)

        x0 = 10
        dx = 1
        fx = 1
        self.assertRaises(TypeError, Itg.root_secant_modified, x0, dx, fx)


if __name__ == '__main__':
    unittest.main()
