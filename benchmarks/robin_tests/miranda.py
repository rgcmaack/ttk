### USAGE
### miranda.py <mode> <extend> <#threads>
### modes == msc / msswall / mssfine / mssfast / msssplit

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

datasetname = 'miranda'
filepath = '/data/miranda_1024x1024x1024_float32.raw'

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
view.ViewSize = [1405, 821]
view.AxesGrid = 'GridAxes3DActor'
view.CenterOfRotation = [127.5, 127.5, 127.5]
view.StereoType = 'Crystal Eyes'
view.CameraPosition = [-122.06393563707445, 892.766610259801, 713.4196031447459]
view.CameraFocalPoint = [127.49999999999996, 127.50000000000001, 127.5]
view.CameraViewUp = [0.36088842401332133, 0.638293441419574, -0.6799566368923377]
view.CameraFocalDisk = 1.0
view.CameraParallelScale = 260.12445969769175
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
data = ImageReader(registrationName='miranda_1024x1024x1024_float32.raw', FileNames=[filepath])
data.DataScalarType = 'float'
data.DataByteOrder = 'LittleEndian'
data.DataExtent = [0, 1023, 0, 1023, 0, 1023]

extent = int(sys.argv[2])

# create a new 'Extract Subset'
data = ExtractSubset(registrationName='subdata', Input=data)
data.VOI = [0, extent, 0, extent, 0, extent]

# create a new 'TTK TopologicalSimplificationByPersistence'
data = TTKTopologicalSimplificationByPersistence(registrationName='simpdata', Input=data)
data.InputArray = ['POINTS', 'ImageFile']
data.PersistenceThreshold = 0.1

if sys.argv[1] == 'msc':
    # create a new 'TTK MorseSmaleComplex'
    msc = TTKMorseSmaleComplex(registrationName='msc', Input=data)
    msc.ScalarField = ['POINTS', 'ImageFile']
    msc.OffsetField = ['POINTS', 'ImageFile']
    msc.CriticalPoints = 0
    msc.Ascending1Separatrices = 0
    msc.Descending1Separatrices = 0
    msc.SaddleConnectors = 0
    msc.Ascending2Separatrices = 1
    msc.Descending2Separatrices = 1

    if len(sys.argv) > 3:
        msc.UseAllCores = 0
        msc.ThreadNumber = int(sys.argv[3])

    # show data from msspl_2
    msc_d = Show(OutputPort(msc, 2), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    msc_d.Representation = 'Surface'
    msc_d.ColorArrayName = [None, '']
    msc_d.SelectTCoordArray = 'None'
    msc_d.SelectNormalArray = 'None'
    msc_d.SelectTangentArray = 'None'
    msc_d.OSPRayScaleFunction = 'PiecewiseFunction'
    msc_d.SelectOrientationVectors = 'None'
    msc_d.ScaleFactor = 25.5
    msc_d.SelectScaleArray = 'None'
    msc_d.GlyphType = 'Arrow'
    msc_d.GlyphTableIndexArray = 'None'
    msc_d.GaussianRadius = 1.2750000000000001
    msc_d.SetScaleArray = [None, '']
    msc_d.ScaleTransferFunction = 'PiecewiseFunction'
    msc_d.OpacityArray = [None, '']
    msc_d.OpacityTransferFunction = 'PiecewiseFunction'
    msc_d.DataAxesGrid = 'GridAxesRepresentation'
    msc_d.PolarAxes = 'PolarAxesRepresentation'
else:
    # create a new 'TTK MorseSmaleSegmentationPL'
    msspl = TTKMorseSmaleSegmentationPL(registrationName='msspl', Input=data)
    msspl.InputArray = ['POINTS', 'ImageFile']
    msspl.a2SeparaticiesMode = computeMode
    if len(sys.argv) > 3:
        msspl.UseAllCores = 0
        msspl.ThreadNumber = int(sys.argv[3])

    # show data from msspl_2
    msspl_d = Show(OutputPort(msspl, 2), view, 'GeometryRepresentation')

    # trace defaults for the display properties.
    msspl_d.Representation = 'Surface'
    msspl_d.ColorArrayName = [None, '']
    msspl_d.SelectTCoordArray = 'None'
    msspl_d.SelectNormalArray = 'None'
    msspl_d.SelectTangentArray = 'None'
    msspl_d.OSPRayScaleFunction = 'PiecewiseFunction'
    msspl_d.SelectOrientationVectors = 'None'
    msspl_d.ScaleFactor = 25.5
    msspl_d.SelectScaleArray = 'None'
    msspl_d.GlyphType = 'Arrow'
    msspl_d.GlyphTableIndexArray = 'None'
    msspl_d.GaussianRadius = 1.2750000000000001
    msspl_d.SetScaleArray = [None, '']
    msspl_d.ScaleTransferFunction = 'PiecewiseFunction'
    msspl_d.OpacityArray = [None, '']
    msspl_d.OpacityTransferFunction = 'PiecewiseFunction'
    msspl_d.DataAxesGrid = 'GridAxesRepresentation'
    msspl_d.PolarAxes = 'PolarAxesRepresentation'

# ---- setup the visualization

# show data from subdata
data_d = Show(data, view, 'UniformGridRepresentation')

lut = GetColorTransferFunction('ImageFile')
lut.RGBPoints = [0.999999463558197, 0.0, 0.0, 0.423499, 1.1024206972986161, 0.0, 0.119346, 0.529237, 1.2048426966392696, 0.0, 0.238691, 0.634976, 1.3072639303796887, 0.0, 0.346852, 0.68788, 1.4096851641201078, 0.0, 0.45022, 0.718141, 1.512106397860527, 0.0, 0.553554, 0.664839, 1.6145283972011804, 0.0, 0.651082, 0.519303, 1.7169494778215528, 0.115841, 0.72479, 0.352857, 1.8193709412420422, 0.326771, 0.781195, 0.140187, 1.921792098422438, 0.522765, 0.798524, 0.0284624, 2.0242140977630916, 0.703162, 0.788685, 0.00885756, 2.1266353315035103, 0.845118, 0.751133, 0.0, 2.22905656524393, 0.955734, 0.690825, 0.0, 2.3314777989843485, 0.995402, 0.567916, 0.0618524, 2.4338997983250024, 0.987712, 0.403398, 0.164851, 2.5311999320983887, 0.980407, 0.247105, 0.262699]
lut.ColorSpace = 'Lab'
lut.ScalarRangeInitialized = 1.0

pfw = GetOpacityTransferFunction('ImageFile')
pfw.Points = [0.999999463558197, 0.0, 0.5, 0.0, 2.5311999320983887, 1.0, 0.5, 0.0]
pfw.ScalarRangeInitialized = 1

# trace defaults for the display properties.
data_d.Representation = 'Surface'
data_d.ColorArrayName = ['POINTS', 'ImageFile']
data_d.LookupTable = lut
data_d.Opacity = 0.1
data_d.SelectTCoordArray = 'None'
data_d.SelectNormalArray = 'None'
data_d.SelectTangentArray = 'None'
data_d.OSPRayScaleArray = 'ImageFile'
data_d.OSPRayScaleFunction = 'PiecewiseFunction'
data_d.SelectOrientationVectors = 'None'
data_d.ScaleFactor = 25.5
data_d.SelectScaleArray = 'ImageFile'
data_d.GlyphType = 'Arrow'
data_d.GlyphTableIndexArray = 'ImageFile'
data_d.GaussianRadius = 1.2750000000000001
data_d.SetScaleArray = ['POINTS', 'ImageFile']
data_d.ScaleTransferFunction = 'PiecewiseFunction'
data_d.OpacityArray = ['POINTS', 'ImageFile']
data_d.OpacityTransferFunction = 'PiecewiseFunction'
data_d.DataAxesGrid = 'GridAxesRepresentation'
data_d.PolarAxes = 'PolarAxesRepresentation'
data_d.ScalarOpacityUnitDistance = 1.7320508075688774
data_d.ScalarOpacityFunction = pfw
data_d.OpacityArrayName = ['POINTS', 'ImageFile']
data_d.IsosurfaceValues = [1.7655996978282928]
data_d.SliceFunction = 'Plane'
data_d.Slice = 127

data_d.ScaleTransferFunction.Points = [0.999999463558197, 0.0, 0.5, 0.0, 2.5311999320983887, 1.0, 0.5, 0.0]
data_d.OpacityTransferFunction.Points = [0.999999463558197, 0.0, 0.5, 0.0, 2.5311999320983887, 1.0, 0.5, 0.0]
data_d.SliceFunction.Origin = [127.5, 127.5, 127.5]

if __name__ == '__main__':
    ResetCamera(view=view)
    SaveScreenshot('/data/img/' + datasetname + '_' + sys.argv[1] + '_' + sys.argv[2] + '.png', layout, SaveAllViews=1, ImageResolution=[1920,1080])
