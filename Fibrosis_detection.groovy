/**
 * Fibrosis quantification — Masson's trichrome sections
 * QuPath version: 0.5.1
 *
 * Workflow:
 *   1. Detect tissue and exclude background
 *   2. Exclude epicardial, endocardial, and large connective tissue regions
 *      via automatic 700 µm inward erosion + low-resolution fibrosis classifier
 *   3. *** MANUAL REVIEW — inspect and adjust excluded regions as needed ***
 *   4. Apply fibrosis pixel classifier to remaining myocardial tissue
 *
 * Run steps 1–2 first. After manual review (step 3), uncomment and run step 4.
 *
 * Required pixel classifiers (must be present in the QuPath project):
 *   "Tissue identification"          — tissue vs. background (14.7 µm/pixel,
 *                                      channel average threshold 225, sigma 5,
 *                                      Gaussian prefilter)
 *   "LowRes_Fibrosis_290824"         — large connective tissue regions (3.6 µm/pixel,
 *                                      RGB, Gaussian features)
 *   "Cardiomyocyte_fibrosis_210824"  — fibrosis vs. cardiomyocytes (0.91 µm/pixel,
 *                                      RGB, Gaussian + weighted deviation)
 *
 * Output measurements added to tissue annotations:
 *   Cardiomyocyte_fibrosis_210824: [class] area µm²
 *   Fibrosis % = Fibrosis area / (Fibrosis area + Cardiomyocyte area)
 *
 * See supplemental methods for classifier training details.
 */

import qupath.lib.gui.commands.Commands
import qupath.lib.roi.RoiTools.CombineOp

// ── Step 1: Tissue detection ──────────────────────────────────────────────────
// Minimum tissue fragment: 1 × 10⁶ µm²; minimum hole to fill: 200,000 µm²

createAnnotationsFromPixelClassifier("Tissue identification", 1000000.0, 200000.0, "SPLIT", "DELETE_EXISTING")
runPlugin('qupath.lib.plugins.objects.FillAnnotationHolesPlugin', '{}')

// ── Step 2: Exclude border connective tissue ──────────────────────────────────
// Erodes tissue annotation inward by 700 µm to isolate the border region.
// A low-resolution fibrosis classifier detects large connective tissue objects
// (> 4,000 µm²) within that border zone. These are then subtracted from the
// tissue annotation to exclude epicardial and endocardial fibrosis.

// Erode inward 700 µm; keep only the border ring; rename temporarily
selectObjectsByClassification("Tissue")
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin',
    '{"radiusMicrons":-700.0,"lineCap":"ROUND","removeInterior":true,"constrainToParent":false}')
def borderClass = getPathClass("Whole_heart")
getSelectedObjects().forEach { it.setPathClass(borderClass) }

// Detect large connective tissue objects in the remaining (inner) tissue area
selectObjectsByClassification("Tissue")
createAnnotationsFromPixelClassifier("LowRes_Fibrosis_290824", 4000.0, 0.0)

// Remove the border ring annotation; restore its class to "Tissue"
selectObjectsByClassification("Tissue")
clearSelectedObjects(true)
clearSelectedObjects()
selectObjectsByClassification("Whole_heart")
getSelectedObjects().forEach { it.setPathClass(getPathClass("Tissue")) }

// Subtract detected connective tissue regions from the tissue annotation
def fibAnnotation  = getAnnotationObjects().find { it.getPathClass() == getPathClass("Fibrosis") }
def tissAnnotation = getAnnotationObjects().find { it.getPathClass() == getPathClass("Tissue") }
getCurrentHierarchy().getSelectionModel().setSelectedObject(fibAnnotation, true)
getCurrentHierarchy().getSelectionModel().setSelectedObject(tissAnnotation, true)
Commands.combineSelectedAnnotations(getCurrentImageData(), CombineOp.SUBTRACT)
fireHierarchyUpdate()

print """
=======================================================================
 MANUAL REVIEW REQUIRED (Step 3)
 Inspect the excluded regions and adjust annotations where needed.
 When satisfied, uncomment Step 4 below and re-run the script.
=======================================================================
"""

// ── Step 4: Fibrosis quantification ──────────────────────────────────────────
// Uncomment after completing manual review.

// selectObjectsByClassification("Tissue")
// addPixelClassifierMeasurements("Cardiomyocyte_fibrosis_210824", "Cardiomyocyte_fibrosis_210824")