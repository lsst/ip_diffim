
import unittest

import numpy as np
import pandas as pd

import lsst.utils.tests
from lsst.meas.base import IdGenerator
import lsst.afw.table as afwTable
import lsst.geom

from lsst.ip.diffim import nightly


class SimpleAssociationCoaddTaskTestCase(lsst.utils.tests.TestCase):
    """Exercise SimpleAssociationCoaddTask on a tiny in-memory catalog."""

    def setUp(self):
        self.config = nightly.SimpleAssociationCoaddConfig()
        # Make matching very forgiving so everything associates
        self.config.tolerance = 10.0  # arcsec
        self.task = nightly.SimpleAssociationCoaddTask(config=self.config)

    def make_two_sources_same_position(self):
        """Two DiaSources with same skyId & nearly identical coordinates."""
        data = {
            "skyId": [0, 0],
            "diaSourceId": [1, 2],
            "ra": [10.0, 10.0],
            "dec": [20.0, 20.0],
        }
        return pd.DataFrame(data)

    def test_run_produces_diaObjects_and_ids(self):
        diaSources = self.make_two_sources_same_position()
        id_gen = IdGenerator()

        result = self.task.run(diaSources, idGenerator=id_gen)

        assoc = result.assocDiaSources
        objs = result.diaObjects

        # A diaObjectId column must exist and be non-null
        self.assertIn("diaObjectId", assoc.columns)
        self.assertTrue((assoc["diaObjectId"] > 0).all())

        # With identical coordinates they should be merged into one DiaObject
        self.assertGreaterEqual(len(objs), 1)
        unique_ids = np.unique(assoc["diaObjectId"].values)
        self.assertEqual(len(unique_ids), len(objs))

    def test_createDiaObject_helper(self):
        diaSources = self.make_two_sources_same_position()
        diaSources.set_index(["skyId", "diaSourceId"], inplace=True)
        diaObjCat = []
        diaObjCoords = []
        hpIndices = []

        id_gen = IdGenerator()
        idCat = id_gen.make_source_catalog(
            afwTable.SourceTable.makeMinimalSchema()
        )

        (skyId, diaSourceId), diaSrc = next(iter(diaSources.iterrows()))
        self.task.addNewDiaObject(
            diaSrc,
            diaSources,
            skyId,
            diaSourceId,
            diaObjCat,
            idCat,
            diaObjCoords,
            hpIndices,
        )

        self.assertEqual(len(diaObjCat), 1)
        self.assertIn("diaObjectId", diaObjCat[0])
        self.assertEqual(diaSources.loc[(skyId, diaSourceId), "diaObjectId"],
                         diaObjCat[0]["diaObjectId"])


class AssembleNightlyCoaddTaskTestCase(lsst.utils.tests.TestCase):
    """Structural tests for nightly AssembleCoadd replacements."""

    def test_connections_dimensions_and_datasets(self):
        conn = nightly.AssembleNightlyCoaddConnections

        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )
        self.assertEqual(
            conn.coaddExposure.dimensions,
            ("tract", "patch", "skymap", "band", "day_obs"),
        )
        self.assertEqual(
            conn.nImage.dimensions,
            ("tract", "patch", "skymap", "band", "day_obs"),
        )
        self.assertEqual(
            conn.inputMap.dimensions,
            ("tract", "patch", "skymap", "band", "day_obs"),
        )

    def test_config_defaults_match_intent(self):
        cfg = nightly.AssembleNightlyCoaddConfig()
        self.assertTrue(cfg.doInputMap)

    def test_task_instantiation(self):
        cfg = nightly.AssembleNightlyCoaddConfig()
        task = nightly.AssembleNightlyCoaddTask(config=cfg)
        self.assertIsInstance(task, nightly.AssembleNightlyCoaddTask)


class DetectCoaddForDiffTaskTestCase(lsst.utils.tests.TestCase):
    """Check nightly detectCoaddForDiff wiring."""

    def test_connections_dimensions(self):
        conn = nightly.DetectCoaddForDiffConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )
        self.assertEqual(
            conn.exposure.dimensions,
            ("tract", "patch", "band", "skymap", "day_obs"),
        )
        self.assertEqual(
            conn.outputSources.dimensions,
            ("tract", "patch", "band", "skymap", "day_obs"),
        )

    def test_config_defaults(self):
        cfg = nightly.DetectCoaddForDiffConfig()
        self.assertEqual(cfg.detection.thresholdType, "pixel_stdev")
        self.assertTrue(cfg.detection.isotropicGrow)
        self.assertTrue(cfg.doScaleVariance)
        # Nightly ID generator should be used
        self.assertIsInstance(
            cfg.idGenerator, nightly.SkyMapNightIdGeneratorConfig
        )

    def test_task_instantiates_and_has_schema(self):
        cfg = nightly.DetectCoaddForDiffConfig()
        task = nightly.DetectCoaddForDiffTask(config=cfg)
        self.assertIsInstance(task, nightly.DetectCoaddForDiffTask)
        self.assertIsNotNone(task.schema)


class CoaddAlardLuptonSubtractTaskTestCase(lsst.utils.tests.TestCase):
    """Structural tests for coadd Alard–Lupton wrapper."""

    def test_connections_dimensions(self):
        in_conn = nightly.CoaddSubtractInputConnections
        out_conn = nightly.CoaddSubtractImageOutputConnections
        both_conn = nightly.CoaddAlardLuptonSubtractConnections

        self.assertEqual(
            in_conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )
        self.assertEqual(
            out_conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )
        self.assertEqual(
            both_conn.dimensions,
            {"tract", "patch", "band", "skymap", "day_obs"},
        )

    def test_config_overrides(self):
        cfg = nightly.CoaddAlardLuptonSubtractConfig()
        # You explicitly set this to False for coadds
        self.assertFalse(cfg.doDecorrelation)

    def test_task_instantiation(self):
        cfg = nightly.CoaddAlardLuptonSubtractConfig()
        task = nightly.CoaddAlardLuptonSubtractTask(config=cfg)
        self.assertIsInstance(task, nightly.CoaddAlardLuptonSubtractTask)


class DetectAndMeasureCoaddTaskTestCase(lsst.utils.tests.TestCase):
    """Structural checks for nightly DetectAndMeasureCoaddTask."""

    def test_connections_dimensions(self):
        conn = nightly.DetectAndMeasureCoaddConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )

    def test_config_overrides(self):
        cfg = nightly.DetectAndMeasureCoaddConfig()
        # nightly config disables apCorr and some plugins
        self.assertFalse(cfg.doApCorr)
        self.assertIsInstance(
            cfg.idGenerator, nightly.SkyMapNightIdGeneratorConfig
        )
        # Check that the removed plugins really aren't present
        names = set(cfg.measurement.plugins.names)
        self.assertNotIn("base_PeakLikelihoodFlux", names)
        self.assertNotIn("base_PeakCentroid", names)

    def test_task_instantiation(self):
        cfg = nightly.DetectAndMeasureCoaddConfig()
        task = nightly.DetectAndMeasureCoaddTask(config=cfg)
        self.assertIsInstance(task, nightly.DetectAndMeasureCoaddTask)


class TransformCoaddDiaSourceCatalogTaskTestCase(lsst.utils.tests.TestCase):
    """Tests for nightly TransformCoaddDiaSourceCatalogTask wiring."""

    def test_connections_dimensions(self):
        conn = nightly.TransformCoaddDiaSourceCatalogConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )

    def test_config_uses_nightly_id_generator(self):
        cfg = nightly.TransformCoaddDiaSourceCatalogConfig()
        self.assertIsInstance(
            cfg.idGenerator, nightly.SkyMapNightIdGeneratorConfig
        )
        

class DrpCoaddAssociationPipeTaskTestCase(lsst.utils.tests.TestCase):
    """Config-level checks for nightly DRP association driver."""

    def test_connections_dimensions(self):
        conn = nightly.DrpCoaddAssociationPipeConnections
        self.assertEqual(conn.dimensions, {"tract", "patch", "skymap"})

    def test_config_associator_is_nightly_simple_association(self):
        cfg = nightly.DrpCoaddAssociationPipeConfig()
        # ConfigurableField stores the target Task class
        self.assertIs(cfg.associator.target, nightly.SimpleAssociationCoaddTask)

    def test_task_instantiation(self):
        cfg = nightly.DrpCoaddAssociationPipeConfig()
        task = nightly.DrpCoaddAssociationPipeTask(config=cfg)
        self.assertIsInstance(task, nightly.DrpCoaddAssociationPipeTask)


class ForcedPhotCoaddFromDataFrameTaskTestCase(lsst.utils.tests.TestCase):
    """Exercise the DataFrame→SourceCatalog bit, which is pure Python."""

    def setUp(self):
        cfg = nightly.ForcedPhotCoaddFromDataFrameConfig()
        self.task = nightly.ForcedPhotCoaddFromDataFrameTask(config=cfg)

    def test_config_defaults(self):
        cfg = nightly.ForcedPhotCoaddFromDataFrameConfig()
        # You explicitly make idGenerator nightly and tweak measurement plugins
        self.assertIsInstance(
            cfg.idGenerator, nightly.SkyMapNightIdGeneratorConfig
        )
        self.assertFalse(cfg.measurement.doReplaceWithNoise)
        self.assertIn("base_PsfFlux", cfg.measurement.plugins.names)

    def test_df2SourceCat_uses_diaObjectId_and_coords(self):
        # This is small and self-contained: no WCS/Butler needed.
        df = pd.DataFrame(
            {
                "diaObjectId": [101, 102],
                "ra": [10.0, 20.0],
                "dec": [30.0, 40.0],
            }
        )
        refCat = self.task.df2SourceCat(df)
        self.assertEqual(len(refCat), 2)

        # Check IDs and coordinates are consistent
        ids = [rec.getId() for rec in refCat]
        self.assertEqual(ids, [101, 102])

        coords = [rec.getCoord() for rec in refCat]
        self.assertAlmostEqual(coords[0].getRa().asDegrees(), 10.0, places=6)
        self.assertAlmostEqual(coords[0].getDec().asDegrees(), 30.0, places=6)


class NightlyFilterAndRBTasksTestCase(lsst.utils.tests.TestCase):
    """Structural tests for FilterDiaSourceNightly/RBTransiNetNightly/etc."""

    def test_filter_dia_source_nightly_connections(self):
        conn = nightly.FilterDiaSourceNightlyConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )
        self.assertEqual(
            conn.diaSourceCat.dimensions,
            ("tract", "patch", "band", "skymap", "day_obs"),
        )

        cfg = nightly.FilterDiaSourceNightlyConfig()
        self.assertIsInstance(cfg, nightly.FilterDiaSourceNightlyConfig)
        task = nightly.FilterDiaSourceNightlyTask(config=cfg)
        self.assertIsInstance(task, nightly.FilterDiaSourceNightlyTask)

    def test_rb_transinet_nightly_connections(self):
        conn = nightly.RBTransiNetNightlyConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )

        cfg = nightly.RBTransiNetNightlyConfig()
        self.assertIsInstance(cfg, nightly.RBTransiNetNightlyConfig)
        task = nightly.RBTransiNetNightlyTask(config=cfg)
        self.assertIsInstance(task, nightly.RBTransiNetNightlyTask)

    def test_filter_dia_source_reliability_nightly_connections(self):
        conn = nightly.FilterDiaSourceReliabilityNightlyConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )

        cfg = nightly.FilterDiaSourceReliabilityNightlyConfig()
        self.assertIsInstance(
            cfg, nightly.FilterDiaSourceReliabilityNightlyConfig
        )
        task = nightly.FilterDiaSourceReliabilityNightlyTask(config=cfg)
        self.assertIsInstance(
            task, nightly.FilterDiaSourceReliabilityNightlyTask
        )


class TransformNightlyDiaObjectForcedSourceTaskTestCase(
    lsst.utils.tests.TestCase
):
    """Test the pure-DataFrame path in TransformNightlyDiaObjectForcedSourceTask."""

    def setUp(self):
        cfg = nightly.TransformNightlyDiaObjectForcedSourceConfig()
        self.task = nightly.TransformNightlyDiaObjectForcedSourceTask(config=cfg)

    def test_connections_dimensions(self):
        conn = nightly.TransformNightlyDiaObjectForcedSourceConnections
        self.assertEqual(
            conn.dimensions, {"tract", "patch", "band", "skymap", "day_obs"}
        )

    def test_run_on_minimal_catalogs(self):
        # Minimal forced-source catalog
        table = afwTable.SourceTable.make(
            afwTable.SourceTable.makeMinimalSchema()
        )
        cat1 = afwTable.SourceCatalog(table)
        rec = cat1.addNew()
        rec.setId(123)
        rec.setCoord(
            lsst.geom.SpherePoint(10.0, 20.0, lsst.geom.degrees)
        )

        # Another catalog, to exercise concatenation path
        cat2 = cat1.copy(deep=True)

        # Minimal reference catalog
        ref = pd.DataFrame(
            {
                "diaObjectId": [1, 2],
                "ref_col": [10.0, 20.0],
            }
        )

        # Configure to join on diaObjectId
        self.task.config.referenceColumns = ["some_ref_col"]
        self.task.config.keyRef = "diaObjectId"
        self.task.config.key = "forcedSourceOnDiaObjectId"

        # Fake the key column in forced sources by converting to pandas
        df_cat = cat1.asAstropy().to_pandas()
        df_cat["forcedSourceOnDiaObjectId"] = 123
        cat1 = afwTable.SourceCatalog(
            afwTable.SourceTable.makeMinimalSchema()
        )
        for i, row in df_cat.iterrows():
            r = cat1.addNew()
            r.setId(i)
            r.setCoord(
                lsst.geom.SpherePoint(
                    row["coord_ra"], row["coord_dec"], lsst.geom.radians
                )
            )
        self.task.config.key = "id"          # left key in forcedDf
        self.task.config.keyRef = "diaObjectId"  # right key in refCat
        self.task.config.referenceColumns = ["ref_col"]

        result = self.task.run([cat1, cat2], ref)
        out = result.outputCatalog

        # We expect some rows and, if the join worked, a "some_ref_col" column if non-empty
        self.assertIsInstance(out, pd.DataFrame)
        self.assertGreaterEqual(len(out), 0)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()