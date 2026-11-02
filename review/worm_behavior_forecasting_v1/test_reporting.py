import unittest
import numpy as np
from summarise import contrast_stats
class ReportingTests(unittest.TestCase):
 def test_fixed_fraction_reduction_exact(self):
  s=contrast_stats([1,2,3],[.5,1,1.5]);self.assertAlmostEqual(s['relative_MSE_reduction_percent'],50);self.assertAlmostEqual(s['CI95_low_percent'],50);self.assertAlmostEqual(s['CI95_high_percent'],50)
 def test_animal_cluster_not_rows(self):
  s=contrast_stats([1,1,2],[.9,.8,2.2],['date1','date1','date2']);self.assertEqual(s['n_animals'],3);self.assertEqual(s['n_bootstrap_clusters'],2)
 def test_negative_results_retained(self):
  s=contrast_stats([1,2,3],[2,4,6]);self.assertAlmostEqual(s['relative_MSE_reduction_percent'],-100);self.assertEqual(s['animals_improved'],0);self.assertLess(s['CI95_high_percent'],0)
 def test_degenerate_losses_rejected(self):
  self.assertRaises(ValueError,contrast_stats,[0,0],[1,1])
if __name__=='__main__':unittest.main()
