/// \ingroup vtk
/// \class ttkmorseSmaleComplexOrder
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::morseSmaleComplexOrder module.
///
/// This VTK filter uses the ttk::morseSmaleComplexOrder module to compute a
/// a MSC variant using the order field defined on the input vtkDataSet.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkPolyData)
/// \param Output1 Output 1-separatrices (vtkPolyData)
/// \param Output2 Output 2-separatrices (vtkPolyData)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// \sa ttk::morseSmaleComplexOrder
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkmorseSmaleComplexOrderModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

/* NAME
 *   ttkmorseSmaleComplexOrder
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <morseSmaleComplexOrder.h>

class vtkPolyData;

using namespace ttk;

class TTKMORSESMALECOMPLEXORDER_EXPORT ttkmorseSmaleComplexOrder
  : public ttkAlgorithm,
    protected ttk::morseSmaleComplexOrder
{
public:
  static ttkmorseSmaleComplexOrder *New();
  vtkTypeMacro(ttkmorseSmaleComplexOrder, ttkAlgorithm);

  vtkSetMacro(ComputeCriticalPoints, bool);
  vtkGetMacro(ComputeCriticalPoints, bool);

  vtkSetMacro(ComputeAscendingSeparatrices1, bool);
  vtkGetMacro(ComputeAscendingSeparatrices1, bool);

  vtkSetMacro(ComputeDescendingSeparatrices1, bool);
  vtkGetMacro(ComputeDescendingSeparatrices1, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkSetMacro(ComputeAscendingSeparatrices2, bool);
  vtkGetMacro(ComputeAscendingSeparatrices2, bool);

  vtkSetMacro(ComputeDescendingSeparatrices2, bool);
  vtkGetMacro(ComputeDescendingSeparatrices2, bool);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

protected:
  ttkmorseSmaleComplexOrder();
  ~ttkmorseSmaleComplexOrder() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputArray,
               vtkPolyData *const outputCriticalPoints,
               vtkPolyData *const outputSeparatrices1,
               vtkPolyData *const outputSeparatrices2,
               const triangulationType &triangulation);

private:
  bool ComputeCriticalPoints{true};
  bool ComputeAscendingSeparatrices1{true};
  bool ComputeDescendingSeparatrices1{true};
  bool ComputeSaddleConnectors{true};
  bool ComputeAscendingSeparatrices2{false};
  bool ComputeDescendingSeparatrices2{false};
  bool ComputeAscendingSegmentation{true};
  bool ComputeDescendingSegmentation{true};
  bool ComputeFinalSegmentation{true};

  // critical points
  std::vector<std::array<float, 3>> criticalPoints_points{};
  std::vector<char> criticalPoints_points_cellDimensions{};
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds{};
};