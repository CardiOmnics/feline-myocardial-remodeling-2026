/**
 * Nucleus detection for myocardial nuclear morphometry
 * QuPath version: 0.5.1
 *
 * Adapted from the CellposeBuilder template script.
 * @original-author Olivier Burri, BIOP, EPFL
 * @source https://github.com/BIOP/qupath-extension-cellpose
 * @license Apache-2.0
 *
 * Detects nuclei in HE-stained myocardial sections using a custom-trained
 * Cellpose model. The hematoxylin channel is isolated via color deconvolution
 * before detection. Detected nuclei are returned as annotations with shape
 * and intensity measurements.
 *
 * Requirements:
 *   QuPath Cellpose extension — https://github.com/BIOP/qupath-extension-cellpose
 *   Custom nucleus model     — available in this repository under models/
 *
 * Usage:
 *   1. Open an HE image in QuPath with color deconvolution stains defined
 *   2. Select the annotation(s) to run detection within
 *   3. Update pathModel below to the local path of the downloaded model
 *   4. Run the script
 *
 * Output:
 *   Annotations with shape measurements (area, circularity, etc.)
 *   and hematoxylin channel intensity measurements per nucleus.
 */

import qupath.ext.biop.cellpose.Cellpose2D

// ── Configuration ─────────────────────────────────────────────────────────────
// Update this path to where you saved the model file from the repository
def pathModel = 'models/cellpose_nucleus_model'

// ── Detection ─────────────────────────────────────────────────────────────────
def stains    = getCurrentImageData().getColorDeconvolutionStains()
def imageData = getCurrentImageData()

def cellpose = Cellpose2D.builder(pathModel)
        .pixelSize(0.2)                                                    // Detection resolution (µm/pixel)
        .preprocess(
            ImageOps.Channels.deconvolve(stains),                          // Color deconvolution
            ImageOps.Channels.extract(0)                                   // Extract hematoxylin channel
        )
        .measureShape()                                                    // Shape measurements per nucleus
        .measureIntensity()                                                // Intensity measurements per nucleus
        .createAnnotations()                                               // Return nuclei as annotations
        .build()

def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("Cellpose", "Please select a parent object!")
    return
}

cellpose.detectObjects(imageData, pathObjects)

println 'Nucleus detection done'