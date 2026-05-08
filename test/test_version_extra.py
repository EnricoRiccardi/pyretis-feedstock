import unittest
from importlib import reload
import pyretis.version


class TestVersionExtra(unittest.TestCase):
    def test_version_not_release(self):
        # Save original value
        original_release = pyretis.version.RELEASE
        original_version = pyretis.version.VERSION

        try:
            # Mock RELEASE to False
            pyretis.version.RELEASE = False
            # Reload module to trigger the 'if not RELEASE' block
            reload(pyretis.version)
            # Check if VERSION is now GIT_VERSION
            self.assertEqual(pyretis.version.VERSION,
                             pyretis.version.GIT_VERSION)
        finally:
            # Restore original value and reload again
            pyretis.version.RELEASE = original_release
            reload(pyretis.version)
            self.assertEqual(pyretis.version.VERSION, original_version)


if __name__ == '__main__':
    unittest.main()
