import unittest
import mongomock

from abcd import ABCD


class Mongo(unittest.TestCase):

    @classmethod
    @mongomock.patch(servers=(('localhost', 27017),))
    def setUpClass(cls):
        url = 'mongodb://localhost'
        abcd = ABCD.from_url(url)
        abcd.print_info()

        cls.abcd = abcd

    @classmethod
    def tearDownClass(cls):
        cls.abcd.destroy()

    def test_thing(self):
        print(self.abcd.info())

    def test_push(self):
        from io import StringIO
        from ase.io import read

        xyz = StringIO("""2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000 
            Si       0.00000000       0.00000000       0.00000000 
            """)

        atoms = read(xyz, format='xyz')
        atoms.set_cell([1, 1, 1])

        self.abcd.destroy()
        self.abcd.push(atoms)
        new = list(self.abcd.get_atoms())[0]

        assert atoms == new
        self.abcd.destroy()


if __name__ == '__main__':
    unittest.main(verbosity=1, exit=False)
