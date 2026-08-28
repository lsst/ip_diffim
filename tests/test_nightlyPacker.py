import unittest

import lsst.utils.tests

try:
    from lsst.daf.butler import DimensionUniverse, DataCoordinate
    HAVE_DAF_BUTLER = True
except ImportError:
    HAVE_DAF_BUTLER = False

from lsst.ip.diffim.nightly import (
    SkyMapNightDimensionPacker,
    SkyMapNightDimensionPackerConfig,
    SkyMapNightIdGeneratorConfig,
)


@unittest.skipUnless(HAVE_DAF_BUTLER, "daf_butler not setup")
class SkyMapNightDimensionPackerTestCase(lsst.utils.tests.TestCase):
    """Tests for the SkyMapNightDimensionPacker.

    Only cover behaviour that's distinct from SkyMapDimensionPacker:
    - extra dimensions day_obs and instrument
    - night axis (n_nights / day_zero)
    - updated dimension validation.
    """

    def setUp(self):
        self.universe = DimensionUniverse()

        # Fixed dataId: {skymap, instrument}, with a skymap record so we can
        # derive n_tracts and n_patches (just like the skymap tests).
        self.fixed = DataCoordinate.from_full_values(
            self.universe.conform(["skymap", "instrument"]),
            values=("unimportant", "TestCam"),
        ).expanded(
            records={
                "skymap": self.universe["skymap"].RecordClass(
                    name="unimportant",
                    tract_max=5,
                    patch_nx_max=3,
                    patch_ny_max=3,
                )
            }
        )

    # ------------------------------------------------------------------
    # Round-trip tests that exercise day_obs + instrument
    # ------------------------------------------------------------------

    def test_with_day_obs_no_filter(self):
        """Pack/unpack with {tract,patch,day_obs,instrument} dimensions."""
        dimensions = self.universe.conform(
            ["skymap", "tract", "patch", "day_obs", "instrument"]
        )
        dataId = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=2,
            patch=6,
            day_obs=20200102,
            universe=self.universe,
        )

        packer = SkyMapNightDimensionPacker(
            self.fixed,
            dimensions=dimensions,
            n_nights=10,
        )
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)

        unpacked = packer.unpack(packedId)
        self.assertEqual(unpacked, dataId)

    def test_with_day_obs_and_filter(self):
        """Pack/unpack with {tract,patch,band,day_obs,instrument} dimensions."""
        dimensions = self.universe.conform(
            ["skymap", "tract", "patch", "band", "day_obs", "instrument"]
        )
        dataId = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=2,
            patch=6,
            band="g",
            day_obs=20200103,
            universe=self.universe,
        )

        packer = SkyMapNightDimensionPacker(
            self.fixed,
            dimensions=dimensions,
            n_nights=10,
        )
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)

        unpacked = packer.unpack(packedId)
        self.assertEqual(unpacked, dataId)

    # ------------------------------------------------------------------
    # Night axis (n_nights / day_zero) tests
    # ------------------------------------------------------------------

    def test_night_range_respected(self):
        """day_obs must be within [day_zero, day_zero + n_nights)."""
        dimensions = self.universe.conform(
            ["skymap", "tract", "patch", "band", "day_obs", "instrument"]
        )
        packer = SkyMapNightDimensionPacker(
            self.fixed,
            dimensions=dimensions,
            n_nights=2,  # only 2020-01-01 and 2020-01-02 are allowed
        )

        # These should be fine.
        dataId1 = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=0,
            patch=0,
            band="g",
            day_obs=20200101,
            universe=self.universe,
        )
        dataId2 = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=0,
            patch=0,
            band="g",
            day_obs=20200102,
            universe=self.universe,
        )

        packed1 = packer.pack(dataId1)
        packed2 = packer.pack(dataId2)
        self.assertNotEqual(packed1, packed2)

        # This is one day beyond the allowed range; should fail.
        dataId3 = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=0,
            patch=0,
            band="g",
            day_obs=20200103,
            universe=self.universe,
        )
        with self.assertRaises(ValueError):
            packer.pack(dataId3)

    # ------------------------------------------------------------------
    # from_config with nights
    # ------------------------------------------------------------------

    def test_from_config_with_nights(self):
        """from_config should honor n_nights and still round-trip day_obs."""
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=2,
            patch=6,
            band="g",
            day_obs=20200105,
            universe=self.universe,
        )
        config = SkyMapNightDimensionPackerConfig()
        config.n_tracts = 5
        config.n_patches = 9
        config.n_bands = 3
        # Use a small, explicit band mapping so this isn't relying on defaults.
        config.bands = {"r": 0, "g": 1, "i": 2}
        config.n_nights = 7

        packer = SkyMapNightDimensionPacker.from_config(data_id, config=config)
        packed_id = packer.pack(data_id)
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packed_id), data_id)

    def test_from_config_no_bands_with_nights(self):
        """from_config should also work in a bandless configuration."""
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=1,
            patch=4,
            day_obs=20200110,
            universe=self.universe,
        )
        config = SkyMapNightDimensionPackerConfig()
        config.n_tracts = 5
        config.n_patches = 9
        # No bands configured; this should create a dimensions set
        # without "band" and still pack/unpack correctly.
        config.n_bands = 0
        config.n_nights = 100

        packer = SkyMapNightDimensionPacker.from_config(data_id, config=config)
        packed_id = packer.pack(data_id)
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packed_id), data_id)


@unittest.skipUnless(HAVE_DAF_BUTLER, "daf_butler not setup")
class SkyMapNightIdGeneratorConfigTestCase(lsst.utils.tests.TestCase):
    """Sanity checks that the nightly IdGenerator config is wired correctly."""

    def setUp(self):
        self.universe = DimensionUniverse()
        self.fixed = DataCoordinate.from_full_values(
            self.universe.conform(["skymap", "instrument"]),
            values=("unimportant", "TestCam"),
        ).expanded(
            records={
                "skymap": self.universe["skymap"].RecordClass(
                    name="unimportant",
                    tract_max=5,
                    patch_nx_max=3,
                    patch_ny_max=3,
                )
            }
        )

    def test_packer_field_type_and_config(self):
        cfg = SkyMapNightIdGeneratorConfig()

        # Field is present
        self.assertTrue(hasattr(cfg, "packer"))

        # Nightly-specific option is exposed on the packer config
        self.assertTrue(hasattr(cfg.packer, "n_nights"))


    def test_make_dimension_packer_uses_nightly_packer(self):
        """_make_dimension_packer should return a SkyMapNightDimensionPacker."""
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            instrument=self.fixed["instrument"],
            tract=0,
            patch=0,
            day_obs=20200101,
            universe=self.universe,
        )
        cfg = SkyMapNightIdGeneratorConfig()
        cfg.packer.n_tracts = 5
        cfg.packer.n_patches = 9
        cfg.packer.n_nights = 100

        packer = cfg._make_dimension_packer(data_id)
        self.assertIsInstance(packer, SkyMapNightDimensionPacker)

        # Quick pack/unpack sanity check.
        packed_id = packer.pack(
            DataCoordinate.standardize(
                skymap=self.fixed["skymap"],
                instrument=self.fixed["instrument"],
                tract=1,
                patch=2,
                day_obs=20200105,
                universe=self.universe,
            )
        )
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
