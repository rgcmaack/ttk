### USAGE
### viscousFingering.py <mode> <#threads>
### modes == msc / msswall / mssfine / mssfast / msssplit

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

datasetname = 'viscousFingering'
filepath = '/data/viscousFingering.vti'
imgSavePath = 'img/'

computeMode = 'Walls'
if sys.argv[1] == 'msswall':
    computeMode = 'Walls'
elif sys.argv[1] == 'mssfine':
    computeMode = 'Separate Basins fine'
elif sys.argv[1] == 'mssfast':
    computeMode = 'Separate Basins fast'
elif sys.argv[1] == 'mssslpit':
    computeMode = 'Experiment'

# setup layouts and views used in the visualization
view = CreateView('RenderView')
view.ViewSize = [1706, 1177]
view.AxesGrid = 'GridAxes3DActor'
view.CenterOfRotation = [0.0, 0.0, 5.0]
view.StereoType = 'Crystal Eyes'
view.CameraPosition = [15.634257889126577, -14.647614283003584, -8.437428609796228]
view.CameraFocalPoint = [-4.932682195184598, 4.350956154788873, 9.88388518478105]
view.CameraViewUp = [-0.4347021842541925, 0.33540424692072573, -0.8357858590286003]
view.CameraFocalDisk = 1.0
view.CameraParallelScale = 8.660254037844387
view.EnableRayTracing = 1
view.BackEnd = 'OSPRay raycaster'
view.Shadows = 1
view.AmbientSamples = 16
view.LightScale = 0.8
view.OSPRayMaterialLibrary = GetMaterialLibrary()

layout = CreateLayout(name='Layout')
layout.AssignView(0, view)
layout.SetSize(1920, 1080)

SetActiveView(view)

# ---- setup the data processing pipelines

# create a new 'Image Reader'
data = XMLImageDataReader(registrationName='viscousFingering.vti', FileName=[filepath])
data.CellArrayStatus = ['vtkGhostType']
data.PointArrayStatus = ['concentration', 'velocity', 'vtkValidPointMask', 'vtkGhostType']
data.TimeArray = 'None'

# create a new 'TTK TopologicalSimplificationByPersistence'
data = TTKTopologicalSimplificationByPersistence(registrationName='TTKTopologicalSimplificationByPersistence1', Input=data)
data.InputArray = ['POINTS', 'concentration']
data.PersistenceThreshold = 100.0

if sys.argv[1] == 'msc':
    # create a new 'TTK MorseSmaleComplex'
    msc = TTKMorseSmaleComplex(registrationName='msc', Input=data)
    msc.ScalarField = ['POINTS', 'concentration']
    msc.OffsetField = ['POINTS', 'concentration']
    msc.CriticalPoints = 0
    msc.Ascending1Separatrices = 0
    msc.Descending1Separatrices = 0
    msc.SaddleConnectors = 0
    msc.Ascending2Separatrices = 1
    msc.Descending2Separatrices = 1

    if len(sys.argv) > 2:
        msc.UseAllCores = 0
        msc.ThreadNumber = int(sys.argv[2])

    clip1 = Clip(registrationName='Clip1', Input=OutputPort(msc,2))
    clip1.ClipType = 'Box'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['CELLS', 'NumberOfCriticalPointsOnBoundary']
    clip1.Value = 25.5

    # init the 'Box' selected for 'ClipType'
    clip1.ClipType.Position = [-5.0, -5.0, 0.0]
    clip1.ClipType.Length = [10.0, 10.0, 10.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [0.0, 0.0, 5.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'view'
    # ----------------------------------------------------------------

    # show data from tTKMorseSmaleComplex1
    tTKMorseSmaleComplex1Display = Show(msc, view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1Display.Representation = 'Surface'
    tTKMorseSmaleComplex1Display.ColorArrayName = [None, '']
    tTKMorseSmaleComplex1Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1Display.OSPRayScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1Display.ScaleFactor = 1.0000000953674317
    tTKMorseSmaleComplex1Display.SelectScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1Display.GlyphTableIndexArray = 'CellDimension'
    tTKMorseSmaleComplex1Display.GaussianRadius = 0.05000000476837158
    tTKMorseSmaleComplex1Display.SetScaleArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.OpacityArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # show data from tTKMorseSmaleComplex1_2
    tTKMorseSmaleComplex1_2Display = Show(OutputPort(msc, 1), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1_2Display.Representation = 'Surface'
    tTKMorseSmaleComplex1_2Display.ColorArrayName = [None, '']
    tTKMorseSmaleComplex1_2Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1_2Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1_2Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1_2Display.OSPRayScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1_2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1_2Display.ScaleFactor = 1.0000000953674317
    tTKMorseSmaleComplex1_2Display.SelectScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1_2Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1_2Display.GlyphTableIndexArray = 'CellDimension'
    tTKMorseSmaleComplex1_2Display.GaussianRadius = 0.05000000476837158
    tTKMorseSmaleComplex1_2Display.SetScaleArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1_2Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.OpacityArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1_2Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1_2Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1_2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1_2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # show data from tTKMorseSmaleComplex1_3
    tTKMorseSmaleComplex1_3Display = Show(OutputPort(msc, 3), view, 'UniformGridRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1_3Display.Representation = 'Outline'
    tTKMorseSmaleComplex1_3Display.ColorArrayName = [None, '']
    tTKMorseSmaleComplex1_3Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1_3Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1_3Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1_3Display.OSPRayScaleArray = 'AscendingManifold'
    tTKMorseSmaleComplex1_3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1_3Display.ScaleFactor = 0.9999999999960001
    tTKMorseSmaleComplex1_3Display.SelectScaleArray = 'AscendingManifold'
    tTKMorseSmaleComplex1_3Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1_3Display.GlyphTableIndexArray = 'AscendingManifold'
    tTKMorseSmaleComplex1_3Display.GaussianRadius = 0.0499999999998
    tTKMorseSmaleComplex1_3Display.SetScaleArray = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleComplex1_3Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.OpacityArray = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleComplex1_3Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1_3Display.PolarAxes = 'PolarAxesRepresentation'
    tTKMorseSmaleComplex1_3Display.ScalarOpacityUnitDistance = 0.13638195335133454
    tTKMorseSmaleComplex1_3Display.OpacityArrayName = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleComplex1_3Display.SliceFunction = 'Plane'
    tTKMorseSmaleComplex1_3Display.Slice = 63

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1_3Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1_3Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    tTKMorseSmaleComplex1_3Display.SliceFunction.Origin = [-2.000000165480742e-11, -2.000000165480742e-11, 4.99999999998]

    # show data from clip1
    clip1Display = Show(clip1, view, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'SeparatrixId'
    separatrixIdLUT = GetColorTransferFunction('SeparatrixId')
    separatrixIdLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 9368.0, 0.865003, 0.865003, 0.865003, 18736.0, 0.705882, 0.0156863, 0.14902]
    separatrixIdLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'SeparatrixId'
    separatrixIdPWF = GetOpacityTransferFunction('SeparatrixId')
    separatrixIdPWF.Points = [0.0, 0.0, 0.5, 0.0, 18736.0, 1.0, 0.5, 0.0]
    separatrixIdPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = ['CELLS', 'SeparatrixId']
    clip1Display.LookupTable = separatrixIdLUT
    clip1Display.SelectTCoordArray = 'None'
    clip1Display.SelectNormalArray = 'None'
    clip1Display.SelectTangentArray = 'None'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.SelectScaleArray = 'None'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'None'
    clip1Display.GaussianRadius = 0.05
    clip1Display.SetScaleArray = [None, '']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = [None, '']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityFunction = separatrixIdPWF
    clip1Display.ScalarOpacityUnitDistance = 0.06005740053042452
    clip1Display.OpacityArrayName = ['CELLS', 'NumberOfCriticalPointsOnBoundary']

else:
    # create a new 'TTK MorseSmaleSegmentationPL'
    msspl = TTKMorseSmaleSegmentationPL(registrationName='msspl', Input=data)
    msspl.InputArray = ['POINTS', 'concentration']
    msspl.a2SeparaticiesMode = computeMode
    if len(sys.argv) > 2:
        msspl.UseAllCores = 0
        msspl.ThreadNumber = int(sys.argv[2])

    # show data from msspl_2
    tTKMorseSmaleSegmentationPL1_1Display = Show(OutputPort(msspl, 1), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1_1Display.Representation = 'Surface'
    tTKMorseSmaleSegmentationPL1_1Display.ColorArrayName = [None, '']
    tTKMorseSmaleSegmentationPL1_1Display.SelectTCoordArray = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.SelectNormalArray = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.SelectTangentArray = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_1Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.ScaleFactor = -2.0000000000000002e+298
    tTKMorseSmaleSegmentationPL1_1Display.SelectScaleArray = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.GlyphType = 'Arrow'
    tTKMorseSmaleSegmentationPL1_1Display.GlyphTableIndexArray = 'None'
    tTKMorseSmaleSegmentationPL1_1Display.GaussianRadius = -1e+297
    tTKMorseSmaleSegmentationPL1_1Display.SetScaleArray = [None, '']
    tTKMorseSmaleSegmentationPL1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_1Display.OpacityArray = [None, '']
    tTKMorseSmaleSegmentationPL1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_1Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # show data from tTKMorseSmaleSegmentationPL1_2
    tTKMorseSmaleSegmentationPL1_2Display = Show(OutputPort(msspl, 2), view, 'GeometryRepresentation')
    
    # get color transfer function/color map for 'MSCIds'
    mSCIdsLUT = GetColorTransferFunction('MSCIds')
    mSCIdsLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 13.155999999999999, 0.0, 0.0, 0.360784313725, 26.22, 0.0, 1.0, 1.0, 39.467999999999996, 0.0, 0.501960784314, 0.0, 52.532, 1.0, 1.0, 0.0, 65.688, 1.0, 0.380392156863, 0.0, 78.844, 0.419607843137, 0.0, 0.0, 92.0, 0.878431372549, 0.301960784314, 0.301960784314]
    mSCIdsLUT.ColorSpace = 'RGB'
    mSCIdsLUT.ScalarRangeInitialized = 1.0
    
    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1_2Display.Representation = 'Surface'
    tTKMorseSmaleSegmentationPL1_2Display.ColorArrayName = ['CELLS', 'MSCIds']
    tTKMorseSmaleSegmentationPL1_2Display.LookupTable = mSCIdsLUT
    tTKMorseSmaleSegmentationPL1_2Display.SelectTCoordArray = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.SelectNormalArray = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.SelectTangentArray = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_2Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.SelectScaleArray = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.GlyphType = 'Arrow'
    tTKMorseSmaleSegmentationPL1_2Display.GlyphTableIndexArray = 'None'
    tTKMorseSmaleSegmentationPL1_2Display.GaussianRadius = 0.05
    tTKMorseSmaleSegmentationPL1_2Display.SetScaleArray = [None, '']
    tTKMorseSmaleSegmentationPL1_2Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_2Display.OpacityArray = [None, '']
    tTKMorseSmaleSegmentationPL1_2Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_2Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_2Display.PolarAxes = 'PolarAxesRepresentation'
    
    # show data from tTKMorseSmaleSegmentationPL1_3
    tTKMorseSmaleSegmentationPL1_3Display = Show(OutputPort(msspl, 3), view, 'UniformGridRepresentation')
    
    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1_3Display.Representation = 'Outline'
    tTKMorseSmaleSegmentationPL1_3Display.ColorArrayName = [None, '']
    tTKMorseSmaleSegmentationPL1_3Display.SelectTCoordArray = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.SelectNormalArray = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.SelectTangentArray = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.OSPRayScaleArray = 'AscendingManifold'
    tTKMorseSmaleSegmentationPL1_3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.ScaleFactor = 0.9999999999960001
    tTKMorseSmaleSegmentationPL1_3Display.SelectScaleArray = 'AscendingManifold'
    tTKMorseSmaleSegmentationPL1_3Display.GlyphType = 'Arrow'
    tTKMorseSmaleSegmentationPL1_3Display.GlyphTableIndexArray = 'AscendingManifold'
    tTKMorseSmaleSegmentationPL1_3Display.GaussianRadius = 0.0499999999998
    tTKMorseSmaleSegmentationPL1_3Display.SetScaleArray = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleSegmentationPL1_3Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.OpacityArray = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleSegmentationPL1_3Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_3Display.PolarAxes = 'PolarAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_3Display.ScalarOpacityUnitDistance = 0.13638195335133454
    tTKMorseSmaleSegmentationPL1_3Display.OpacityArrayName = ['POINTS', 'AscendingManifold']
    tTKMorseSmaleSegmentationPL1_3Display.SliceFunction = 'Plane'
    tTKMorseSmaleSegmentationPL1_3Display.Slice = 63
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleSegmentationPL1_3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 21887.0, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleSegmentationPL1_3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 21887.0, 1.0, 0.5, 0.0]
    
    # init the 'Plane' selected for 'SliceFunction'
    tTKMorseSmaleSegmentationPL1_3Display.SliceFunction.Origin = [-2.000000165480742e-11, -2.000000165480742e-11, 4.99999999998]
    
    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------
    
    # get opacity transfer function/opacity map for 'MSCIds'
    mSCIdsPWF = GetOpacityTransferFunction('MSCIds')
    mSCIdsPWF.Points = [0.0, 0.0, 0.5, 0.0, 92.0, 1.0, 0.5, 0.0]
    mSCIdsPWF.ScalarRangeInitialized = 1

if __name__ == '__main__':
    ResetCamera(view=view)
    SaveScreenshot(imgSavePath + datasetname + '_' + sys.argv[1] +'.png', layout, SaveAllViews=1, ImageResolution=[1920,1080])
