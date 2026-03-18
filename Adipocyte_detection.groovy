/**
 * Adipose tissue quantification and adipocyte size analysis
 * QuPath version: 0.5.x
 *
 * Requires the following pixel classifiers to be present in the QuPath project:
 *   - Tissue_detection_131124
 *   - Adipocyte_131124
 *   - Adipocyte_area_determination_131124
 *   - Adipocyte_membrane_131124
 *
 * Parameters:
 *   createAnnotationsFromPixelClassifier: min object size 10,000 µm², min hole size 3,000 µm²
 *   createDetectionsFromPixelClassifier:  min object size 70 µm²,     min hole size 30 µm²
 *
 * See supplemental methods (Section 1) for classifier training details.
 */

selectAnnotations();
addPixelClassifierMeasurements("Tissue_detection_131124", "Tissue_detection_131124")
addPixelClassifierMeasurements("Adipocyte_131124", "Adipocyte_131124")
createAnnotationsFromPixelClassifier("Adipocyte_area_determination_131124", 10000.0, 3000.0, "SELECT_NEW")
createDetectionsFromPixelClassifier("Adipocyte_membrane_131124", 70.0, 30.0, "SPLIT", "DELETE_EXISTING")
selectDetections();
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER")
