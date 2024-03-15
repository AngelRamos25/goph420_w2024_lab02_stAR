import unittest
from goph420_lab02 import integration as Itg


class TestGaussLegendre(unittest.TestCase):

    def test_Newton_Raphson(self):
        x0 = 1
        fx = Itg.functions_Test()

        result = Itg.root_newton_raphson(x0, fx, dfdx,)
        self.assertEqual(result, 60)

    def test_raises_gauss(self):
        # Test raises:
        f = Itg.GaussDist()

        # Raise for len(lims) != 2
        lims = [0]
        npts = 3
        self.assertRaises(ValueError, Itg.integrate_gauss, f, lims, npts)


if __name__ == '__main__':
    unittest.main()
