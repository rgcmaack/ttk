/// \ingroup base
/// \class ttk::morseSmaleComplexOrder
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.

#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <unordered_map>

namespace ttk {

  /**
   * The morseSmaleComplexOrder class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class morseSmaleComplexOrder : virtual public Debug {

  public:
    morseSmaleComplexOrder();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeDescendingManifold(dataType *outputData,
                        const dataType *inputData,
                        const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;
      int iterations = 1;
      SimplexId numMaxima = 0;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute MSC variant
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Descending Manifold",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_);

        /* compute the Descending Maifold iterating over each vertex, searching
         * the biggest neighbor and compressing its path to its maximum */
        SimplexId nVertices = triangulation->getNumberOfVertices();
        std::unordered_map<SimplexId, SimplexId> maxima;
        
        //std::vector<SimplexId> descendingManifold(nVertices, -1);

        bool globalChange = true; // paths compressed by any thread

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel num_threads(this->threadNumber_) \
                shared(globalChange, iterations, maxima, numMaxima)
        #endif
        {
          bool localChange = false; // paths compressed by the local thread

          #ifdef TTK_ENABLE_OPENMP
          #pragma omp for schedule(dynamic)
          #endif
          // find the biggest neighbor for each vertex
          for(SimplexId i = 0; i < nVertices; i++) {
            SimplexId numNeighbors =
              triangulation->getVertexNeighborNumber(i);
            SimplexId neighborId;
            outputData[i] = i;
            bool hasBiggerNeighbor = false;
            // check all neighbors
            for(SimplexId n = 0; n < numNeighbors; n++) {
              triangulation->getVertexNeighbor(i, n, neighborId);
              if(inputData[neighborId] > inputData[i]) {
                outputData[i] = neighborId;
                hasBiggerNeighbor = true;
              }
            } 
            if(!hasBiggerNeighbor) {
              #ifdef TTK_ENABLE_OPENMP
              #pragma omp critical
              #endif      
              maxima.insert( {i, numMaxima} );
              
              #ifdef TTK_ENABLE_OPENMP
              #pragma omp atomic update
              #endif     
              numMaxima += 1;
            }
          }

          #ifdef TTK_ENABLE_OPENMP
          #pragma omp barrier
          #endif

          #ifdef TTK_ENABLE_OPENMP
          #pragma omp single
          {
            this->printMsg("Computed Maxima",
              0.1, // progress
              localTimer.getElapsedTime(), this->threadNumber_);
          }
          #endif

          // compress paths until no changes occur
          while(globalChange) {
            #ifdef TTK_ENABLE_OPENMP
            #pragma omp barrier
            #endif

            #ifdef TTK_ENABLE_OPENMP
            #pragma omp single
            #endif
            {
              globalChange = false;
            }

            localChange = false;

            #ifdef TTK_ENABLE_OPENMP
            #pragma omp single
            {
              std::string msgdbg = "iteration: " + std::to_string(iterations);
              this->printMsg(msgdbg,
                0.12, // progress
                localTimer.getElapsedTime(), this->threadNumber_);
            }
            #endif

            #ifdef TTK_ENABLE_OPENMP
            #pragma omp for schedule(dynamic)
            #endif
            for(SimplexId i = 0; i < nVertices; i++) {
              // maxima or pointing to maxima
              if(i == (SimplexId)outputData[i]
              || outputData[i] == outputData[(SimplexId)outputData[i]]) {
                continue;
              }
              // path compression
              else {
                outputData[i] = outputData[(SimplexId)outputData[i]];
                localChange = true;
              }
            }
          
            #ifdef TTK_ENABLE_OPENMP
            #pragma omp atomic update
            #endif
              globalChange |= localChange;
            
            #ifdef TTK_ENABLE_OPENMP
            #pragma omp single
            #endif
            {
              iterations += 1;
            }
            
            #ifdef TTK_ENABLE_OPENMP
            #pragma omp barrier
            #endif
          }

          #ifdef TTK_ENABLE_OPENMP
          #pragma omp single
          {
            this->printMsg("Compressed Paths",
              0.1, // progress
              localTimer.getElapsedTime(), this->threadNumber_);
          }
          #endif

          #ifdef TTK_ENABLE_OPENMP
          #pragma omp for schedule(dynamic, 4)
          #endif
          for(SimplexId i = 0; i < nVertices; i++) {
            outputData[i] = maxima[outputData[i]];
          }
        }

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computed Descending Manifold",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator

        const std::string maxMsg = "#Maxima: " + std::to_string(numMaxima);
        this->printMsg(maxMsg, 1, globalTimer.getElapsedTime() );

        const std::string itMsg = "#Iterations: " + std::to_string(iterations);
        this->printMsg(itMsg, 1, globalTimer.getElapsedTime() );

        this->printMsg("Completed", 1, globalTimer.getElapsedTime() );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  }; // morseSmaleComplexOrder class

} // namespace ttk
