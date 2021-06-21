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
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/morseSmaleComplexOrder/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
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

class TTKMORSESMALECOMPLEXORDER_EXPORT ttkmorseSmaleComplexOrder
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::morseSmaleComplexOrder // and we inherit from the base class
{
private:
  std::string OutputArrayName{"Descending Manifold"};

public:
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  static ttkmorseSmaleComplexOrder *New();
  vtkTypeMacro(ttkmorseSmaleComplexOrder, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkmorseSmaleComplexOrder();
  ~ttkmorseSmaleComplexOrder() override;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
