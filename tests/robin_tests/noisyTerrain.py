### USAGE
### noisyTerrain.py <mode> <#threads>
### mode = msc / mss

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

datasetname = 'noisyTerrain'
filepath = '/data/noisyTerrain.vtu'
imgSavePath = 'img/'

# setup layouts and views used in the visualization
view = CreateView('RenderView')
view.ViewSize = [1920, 1080]
view.AxesGrid = 'GridAxes3DActor'
view.OrientationAxesVisibility = 0
view.CenterOfRotation = [0.0, 0.0, -0.7406161576509476]
view.StereoType = 'Crystal Eyes'
view.CameraPosition = [0.0, 0.0, 9.75045632631921]
view.CameraFocalPoint = [0.0, 0.0, -17.83267207840769]
view.CameraFocalDisk = 1.0
view.CameraParallelScale = 7.139038954651633
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
data = XMLUnstructuredGridReader(registrationName=datasetname, FileName=[filepath])
data.PointArrayStatus = ['TextureCoordinates', 'RandomPointScalars', 'Sine', 'DistanceField', 'Blend']
data.TimeArray = 'None'

# create a new 'TTK TopologicalSimplificationByPersistence'
data = TTKTopologicalSimplificationByPersistence(registrationName='simpdata', Input=data)
data.InputArray = ['POINTS', 'Blend']
data.PersistenceThreshold = 0.55

if sys.argv[1] == 'msc':
    # create a new 'TTK MorseSmaleComplex'
    msc = TTKMorseSmaleComplex(registrationName='msc', Input=data)
    msc.ScalarField = ['POINTS', 'Blend']
    msc.OffsetField = ['POINTS', 'Blend']
    if len(sys.argv) > 2:
        msc.UseAllCores = 0
        msc.ThreadNumber = int(sys.argv[2])
   
    tTKMorseSmaleComplex1Display = Show(msc, view, 'GeometryRepresentation')

    # get color transfer function/color map for 'CellDimension'
    cellDimensionLUT = GetColorTransferFunction('CellDimension')
    cellDimensionLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.0, 0.865003, 0.865003, 0.865003, 2.0, 0.705882, 0.0156863, 0.14902]
    cellDimensionLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1Display.Representation = 'Surface'
    tTKMorseSmaleComplex1Display.ColorArrayName = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1Display.LookupTable = cellDimensionLUT
    tTKMorseSmaleComplex1Display.PointSize = 7.0
    tTKMorseSmaleComplex1Display.RenderPointsAsSpheres = 1
    tTKMorseSmaleComplex1Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1Display.OSPRayScaleArray = 'Blend'
    tTKMorseSmaleComplex1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1Display.SelectScaleArray = 'Blend'
    tTKMorseSmaleComplex1Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1Display.GlyphTableIndexArray = 'Blend'
    tTKMorseSmaleComplex1Display.GaussianRadius = 0.05
    tTKMorseSmaleComplex1Display.SetScaleArray = ['POINTS', 'Blend']
    tTKMorseSmaleComplex1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.OpacityArray = ['POINTS', 'Blend']
    tTKMorseSmaleComplex1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1Display.ScaleTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1Display.OpacityTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # show data from tTKMorseSmaleComplex1_1
    tTKMorseSmaleComplex1_1Display = Show(OutputPort(msc, 1), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1_1Display.Representation = 'Surface'
    tTKMorseSmaleComplex1_1Display.ColorArrayName = [None, '']
    tTKMorseSmaleComplex1_1Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1_1Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1_1Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1_1Display.OSPRayScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_1Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1_1Display.SelectScaleArray = 'CellDimension'
    tTKMorseSmaleComplex1_1Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1_1Display.GlyphTableIndexArray = 'CellDimension'
    tTKMorseSmaleComplex1_1Display.GaussianRadius = 0.05
    tTKMorseSmaleComplex1_1Display.SetScaleArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_1Display.OpacityArray = ['POINTS', 'CellDimension']
    tTKMorseSmaleComplex1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_1Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1_1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1_1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1_1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

    # show data from tTKMorseSmaleComplex1_2
    tTKMorseSmaleComplex1_2Display = Show(OutputPort(msc, 2), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1_2Display.Representation = 'Surface'
    tTKMorseSmaleComplex1_2Display.ColorArrayName = [None, '']
    tTKMorseSmaleComplex1_2Display.SelectTCoordArray = 'None'
    tTKMorseSmaleComplex1_2Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1_2Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1_2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1_2Display.ScaleFactor = -0.2
    tTKMorseSmaleComplex1_2Display.SelectScaleArray = 'None'
    tTKMorseSmaleComplex1_2Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1_2Display.GlyphTableIndexArray = 'None'
    tTKMorseSmaleComplex1_2Display.GaussianRadius = -0.01
    tTKMorseSmaleComplex1_2Display.SetScaleArray = [None, '']
    tTKMorseSmaleComplex1_2Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.OpacityArray = [None, '']
    tTKMorseSmaleComplex1_2Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_2Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1_2Display.PolarAxes = 'PolarAxesRepresentation'

    # show data from tTKMorseSmaleSegmentationPL1_3
    tTKMorseSmaleComplex1_3Display = Show(OutputPort(msc, 3), view, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'MorseSmaleManifold'
    morseSmaleManifoldLUT = GetColorTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 5.433999999999999, 0.0, 0.0, 0.360784313725, 10.829999999999998, 0.0, 1.0, 1.0, 16.302, 0.0, 0.501960784314, 0.0, 21.697999999999997, 1.0, 1.0, 0.0, 27.131999999999998, 1.0, 0.380392156863, 0.0, 32.566, 0.419607843137, 0.0, 0.0, 38.0, 0.878431372549, 0.301960784314, 0.301960784314]
    morseSmaleManifoldLUT.ColorSpace = 'RGB'
    morseSmaleManifoldLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'MorseSmaleManifold'
    morseSmaleManifoldPWF = GetOpacityTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldPWF.Points = [0.0, 0.0, 0.5, 0.0, 38.0, 1.0, 0.5, 0.0]
    morseSmaleManifoldPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    tTKMorseSmaleComplex1_3Display.Representation = 'Surface'
    tTKMorseSmaleComplex1_3Display.ColorArrayName = ['POINTS', 'MorseSmaleManifold']
    tTKMorseSmaleComplex1_3Display.LookupTable = morseSmaleManifoldLUT
    tTKMorseSmaleComplex1_3Display.SelectTCoordArray = 'TextureCoordinates'
    tTKMorseSmaleComplex1_3Display.SelectNormalArray = 'None'
    tTKMorseSmaleComplex1_3Display.SelectTangentArray = 'None'
    tTKMorseSmaleComplex1_3Display.OSPRayScaleArray = 'Blend'
    tTKMorseSmaleComplex1_3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleComplex1_3Display.ScaleFactor = 0.1
    tTKMorseSmaleComplex1_3Display.SelectScaleArray = 'Blend'
    tTKMorseSmaleComplex1_3Display.GlyphType = 'Arrow'
    tTKMorseSmaleComplex1_3Display.GlyphTableIndexArray = 'Blend'
    tTKMorseSmaleComplex1_3Display.GaussianRadius = 0.005
    tTKMorseSmaleComplex1_3Display.SetScaleArray = ['POINTS', 'Blend']
    tTKMorseSmaleComplex1_3Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.OpacityArray = ['POINTS', 'Blend']
    tTKMorseSmaleComplex1_3Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleComplex1_3Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleComplex1_3Display.PolarAxes = 'PolarAxesRepresentation'
    tTKMorseSmaleComplex1_3Display.ScalarOpacityFunction = morseSmaleManifoldPWF
    tTKMorseSmaleComplex1_3Display.ScalarOpacityUnitDistance = 0.02599714950829805
    tTKMorseSmaleComplex1_3Display.OpacityArrayName = ['POINTS', 'Blend']

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleComplex1_3Display.ScaleTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleComplex1_3Display.OpacityTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get opacity transfer function/opacity map for 'CellDimension'
    cellDimensionPWF = GetOpacityTransferFunction('CellDimension')
    cellDimensionPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
    cellDimensionPWF.ScalarRangeInitialized = 1

    # get opacity transfer function/opacity map for 'MorseSmaleManifold'
    morseSmaleManifoldPWF = GetOpacityTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldPWF.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]
    morseSmaleManifoldPWF.ScalarRangeInitialized = 1
else:
    # create a new 'TTK MorseSmaleSegmentationPL'
    msspl = TTKMorseSmaleSegmentationPL(registrationName='msspl', Input=data)
    msspl.InputArray = ['POINTS', 'Blend']
    msspl.ComputeSaddles = 1
    if len(sys.argv) > 2:
        msspl.UseAllCores = 0
        msspl.ThreadNumber = int(sys.argv[2])

    # create a new 'TTK MorseSmaleSegmentationPL'
    tTKMorseSmaleSegmentationPL1Display = Show(msspl, view, 'GeometryRepresentation')

    # get color transfer function/color map for 'CriticalityIndex'
    criticalityIndexLUT = GetColorTransferFunction('CriticalityIndex')
    criticalityIndexLUT.InterpretValuesAsCategories = 1
    criticalityIndexLUT.AnnotationsInitialized = 1
    criticalityIndexLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.5, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902]
    criticalityIndexLUT.ScalarRangeInitialized = 1.0
    criticalityIndexLUT.Annotations = ['0', '', '1', '\x01', '3', '\x03']
    criticalityIndexLUT.IndexedColors = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
    criticalityIndexLUT.IndexedOpacities = [1.0, 1.0, 1.0]

    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1Display.Representation = 'Surface'
    tTKMorseSmaleSegmentationPL1Display.ColorArrayName = ['POINTS', 'Criticality Index']
    tTKMorseSmaleSegmentationPL1Display.LookupTable = criticalityIndexLUT
    tTKMorseSmaleSegmentationPL1Display.PointSize = 7.0
    tTKMorseSmaleSegmentationPL1Display.RenderPointsAsSpheres = 1
    tTKMorseSmaleSegmentationPL1Display.SelectTCoordArray = 'None'
    tTKMorseSmaleSegmentationPL1Display.SelectNormalArray = 'None'
    tTKMorseSmaleSegmentationPL1Display.SelectTangentArray = 'None'
    tTKMorseSmaleSegmentationPL1Display.OSPRayScaleArray = 'Blend_Order'
    tTKMorseSmaleSegmentationPL1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleSegmentationPL1Display.SelectScaleArray = 'Blend_Order'
    tTKMorseSmaleSegmentationPL1Display.GlyphType = 'Arrow'
    tTKMorseSmaleSegmentationPL1Display.GlyphTableIndexArray = 'Blend_Order'
    tTKMorseSmaleSegmentationPL1Display.GaussianRadius = 0.05
    tTKMorseSmaleSegmentationPL1Display.SetScaleArray = ['POINTS', 'Blend_Order']
    tTKMorseSmaleSegmentationPL1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1Display.OpacityArray = ['POINTS', 'Blend_Order']
    tTKMorseSmaleSegmentationPL1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleSegmentationPL1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleSegmentationPL1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9006000.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleSegmentationPL1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9006000.0, 1.0, 0.5, 0.0]

    # show data from tTKMorseSmaleSegmentationPL1_1
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

    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1_2Display.Representation = 'Surface'
    tTKMorseSmaleSegmentationPL1_2Display.ColorArrayName = [None, '']
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
    tTKMorseSmaleSegmentationPL1_3Display = Show(OutputPort(msspl, 3), view, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'MorseSmaleManifold'
    morseSmaleManifoldLUT = GetColorTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 5.433999999999999, 0.0, 0.0, 0.360784313725, 10.829999999999998, 0.0, 1.0, 1.0, 16.302, 0.0, 0.501960784314, 0.0, 21.697999999999997, 1.0, 1.0, 0.0, 27.131999999999998, 1.0, 0.380392156863, 0.0, 32.566, 0.419607843137, 0.0, 0.0, 38.0, 0.878431372549, 0.301960784314, 0.301960784314]
    morseSmaleManifoldLUT.ColorSpace = 'RGB'
    morseSmaleManifoldLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'MorseSmaleManifold'
    morseSmaleManifoldPWF = GetOpacityTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldPWF.Points = [0.0, 0.0, 0.5, 0.0, 38.0, 1.0, 0.5, 0.0]
    morseSmaleManifoldPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    tTKMorseSmaleSegmentationPL1_3Display.Representation = 'Surface'
    tTKMorseSmaleSegmentationPL1_3Display.ColorArrayName = ['POINTS', 'MorseSmaleManifold']
    tTKMorseSmaleSegmentationPL1_3Display.LookupTable = morseSmaleManifoldLUT
    tTKMorseSmaleSegmentationPL1_3Display.SelectTCoordArray = 'TextureCoordinates'
    tTKMorseSmaleSegmentationPL1_3Display.SelectNormalArray = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.SelectTangentArray = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.OSPRayScaleArray = 'Blend'
    tTKMorseSmaleSegmentationPL1_3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.SelectOrientationVectors = 'None'
    tTKMorseSmaleSegmentationPL1_3Display.ScaleFactor = 0.1
    tTKMorseSmaleSegmentationPL1_3Display.SelectScaleArray = 'Blend'
    tTKMorseSmaleSegmentationPL1_3Display.GlyphType = 'Arrow'
    tTKMorseSmaleSegmentationPL1_3Display.GlyphTableIndexArray = 'Blend'
    tTKMorseSmaleSegmentationPL1_3Display.GaussianRadius = 0.005
    tTKMorseSmaleSegmentationPL1_3Display.SetScaleArray = ['POINTS', 'Blend']
    tTKMorseSmaleSegmentationPL1_3Display.ScaleTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.OpacityArray = ['POINTS', 'Blend']
    tTKMorseSmaleSegmentationPL1_3Display.OpacityTransferFunction = 'PiecewiseFunction'
    tTKMorseSmaleSegmentationPL1_3Display.DataAxesGrid = 'GridAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_3Display.PolarAxes = 'PolarAxesRepresentation'
    tTKMorseSmaleSegmentationPL1_3Display.ScalarOpacityFunction = morseSmaleManifoldPWF
    tTKMorseSmaleSegmentationPL1_3Display.ScalarOpacityUnitDistance = 0.02599714950829805
    tTKMorseSmaleSegmentationPL1_3Display.OpacityArrayName = ['POINTS', 'Blend']

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tTKMorseSmaleSegmentationPL1_3Display.ScaleTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tTKMorseSmaleSegmentationPL1_3Display.OpacityTransferFunction.Points = [-34.46813385051024, 0.0, 0.5, 0.0, 4.843487179864361, 1.0, 0.5, 0.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for criticalityIndexLUT in view view
    criticalityIndexLUTColorBar = GetScalarBar(criticalityIndexLUT, view)
    criticalityIndexLUTColorBar.WindowLocation = 'Upper Right Corner'
    criticalityIndexLUTColorBar.Title = 'Criticality Index'
    criticalityIndexLUTColorBar.ComponentTitle = ''

    # get color legend/bar for morseSmaleManifoldLUT in view view
    morseSmaleManifoldLUTColorBar = GetScalarBar(morseSmaleManifoldLUT, view)
    morseSmaleManifoldLUTColorBar.Title = 'MorseSmaleManifold'
    morseSmaleManifoldLUTColorBar.ComponentTitle = ''

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get opacity transfer function/opacity map for 'CriticalityIndex'
    criticalityIndexPWF = GetOpacityTransferFunction('CriticalityIndex')
    criticalityIndexPWF.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]
    criticalityIndexPWF.ScalarRangeInitialized = 1

    # get opacity transfer function/opacity map for 'MorseSmaleManifold'
    morseSmaleManifoldPWF = GetOpacityTransferFunction('MorseSmaleManifold')
    morseSmaleManifoldPWF.Points = [0.0, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]
    morseSmaleManifoldPWF.ScalarRangeInitialized = 1

if __name__ == '__main__':
    ResetCamera(view=view)
    SaveScreenshot(imgSavePath + datasetname + '_' + sys.argv[1] + '.png', layout, SaveAllViews=1, ImageResolution=[1920,1080])