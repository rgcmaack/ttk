#/// \ingroup vtk
/// \class ttkMorseSmaleComplexQuasi
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MorseSmaleComplexQuasi module.
///
/// This VTK filter uses the ttk::MorseSmaleComplexQuasi module to compute a
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
/// \sa ttk::MorseSmaleComplexQuasi
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMorseSmaleComplexQuasiModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

/* NAME
 *   ttkMorseSmaleComplexQuasi
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <MorseSmaleComplexQuasi.h>

class vtkPolyData;

using namespace ttk;

class TTKMORSESMALECOMPLEXQUASI_EXPORT ttkMorseSmaleComplexQuasi
  : public ttkAlgorithm,
    protected ttk::MorseSmaleComplexQuasi
{
public:
  static ttkMorseSmaleComplexQuasi *New();
  vtkTypeMacro(ttkMorseSmaleComplexQuasi, ttkAlgorithm);

  ttkSetEnumMacro(SeparatricesManifold, SEPARATRICES_MANIFOLD);
  vtkGetEnumMacro(SeparatricesManifold, SEPARATRICES_MANIFOLD);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

  vtkSetMacro(ComputeSeparatrices, bool);
  vtkGetMacro(ComputeSeparatrices, bool);

  vtkSetMacro(Fast1Separatrices, bool);
  vtkGetMacro(Fast1Separatrices, bool);

  vtkSetMacro(Fast2Separatrices, bool);
  vtkGetMacro(Fast2Separatrices, bool);

protected:
  ttkMorseSmaleComplexQuasi();
  ~ttkMorseSmaleComplexQuasi() override;

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
  bool ComputeSeparatrices{true};
  bool ComputeAscendingSegmentation{true};
  bool ComputeDescendingSegmentation{true};
  bool ComputeFinalSegmentation{true};
  bool Fast1Separatrices{false};
  bool Fast2Separatrices{false};
  SEPARATRICES_MANIFOLD SeparatricesManifold
    {SEPARATRICES_MANIFOLD::MORSESMALE};

  // critical points
  std::vector<std::array<float, 3>> criticalPoints_points{};
  std::vector<char> criticalPoints_points_cellDimensions{};
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds{};

  // 1-separatrices data
  std::vector<float> separatrices1_points{};
  std::vector<char> separatrices1_points_smoothingMask{};
  std::vector<char> separatrices1_points_cellDimensions{};
  std::vector<ttk::SimplexId> separatrices1_points_cellIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_connectivity{};
  std::vector<ttk::SimplexId> separatrices1_cells_sourceIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_destinationIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds{};
  std::vector<char> separatrices1_cells_separatrixTypes{};
  std::vector<char> separatrices1_cells_isOnBoundary{};

  // 2-separatrices data
  std::vector<float> separatrices2_points{};
  std::vector<ttk::SimplexId> separatrices2_cells_connectivity{};
  std::vector<ttk::SimplexId> separatrices2_cells_mscIds{};
  std::vector<ttk::SimplexId> separatrices2_cells_caseTypes{};
};