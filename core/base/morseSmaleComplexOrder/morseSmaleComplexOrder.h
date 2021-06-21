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

        // vertices that may still be compressed
        std::vector<SimplexId>* activeVertices
          = new std::vector<SimplexId>();
        SimplexId nActiveVertices;

#ifndef TTK_ENABLE_OPENMP
        // find maxima and intialize vector of not fully compressed vertices
        for(SimplexId i = 0; i < nVertices; i++) {
          SimplexId neighborId;
          SimplexId numNeighbors =
            triangulation->getVertexNeighborNumber(i);
          bool hasBiggerNeighbor = false;

          outputData[i] = i;

          // check all neighbors
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation->getVertexNeighbor(i, n, neighborId);
            if(inputData[neighborId] > inputData[(SimplexId)outputData[i]]) {
              outputData[i] = neighborId;
              hasBiggerNeighbor = true;
            }
          }
          if(hasBiggerNeighbor) {
            activeVertices->push_back(i);
          }
          else {
            maxima.insert( {i, numMaxima++} );
          }
        }

        nActiveVertices = activeVertices->size();
        std::vector<SimplexId>* newActiveVert;

        // compress paths until no changes occur
        while(nActiveVertices > 0) {
          newActiveVert = new std::vector<SimplexId>();

          std::string msgit = "Iteration: " + std::to_string(iterations) +
            "(" + std::to_string(nActiveVertices) + "/" +
            std::to_string(nVertices) + ")";
          this->printMsg(msgit, 1.0 - ((double)nActiveVertices / nVertices),
            localTimer.getElapsedTime(), this->threadNumber_);
            
          for(SimplexId i = 0; i < nActiveVertices; i++) {
            SimplexId v = activeVertices->at(i);
            dataType &vo = outputData[v];
            // compress path
            vo = outputData[(SimplexId)vo];

            if(vo != outputData[(SimplexId)vo]) {
              newActiveVert->push_back(v);
            }
          }

          iterations += 1;
          delete activeVertices;
          activeVertices = &*newActiveVert;

          nActiveVertices = activeVertices->size();
        }

        delete activeVertices;

        for(SimplexId i = 0; i < nVertices; i++) {
          outputData[i] = maxima[outputData[i]];
        }
#else
/*        bool globalChange = true; // paths compressed by any thread

        #pragma omp parallel num_threads(this->threadNumber_) \
                shared(globalChange, iterations, maxima, numMaxima, \
                activeVertices, nActiveVertices)
        {
          bool localChange = false; // paths compressed by the local thread
          std::vector<SimplexId> lActiveVertices;

          #pragma omp for schedule(dynamic)
          // find the biggest neighbor for each vertex
          for(SimplexId i = 0; i < nVertices; i++) {
            SimplexId neighborId;
            SimplexId numNeighbors =
              triangulation->getVertexNeighborNumber(i);
            bool hasBiggerNeighbor = false;

            outputData[i] = i;

            // check all neighbors
            for(SimplexId n = 0; n < numNeighbors; n++) {
              triangulation->getVertexNeighbor(i, n, neighborId);
              if(inputData[neighborId] > outputData[i]) {
                outputData[i] = neighborId;
                hasBiggerNeighbor = true;
              }
            }

            if(hasBiggerNeighbor) {
              if(outputData[i] != outputData[outputData[i]]) {
                lActiveVertices.push_back(i);
              }
            }
            else {
              #pragma omp critical
              maxima.insert( {i, numMaxima} );

              #pragma omp atomic update
              numMaxima += 1;
            }
          }

          #pragma omp critical
          activeVertices->insert(activeVertices->end(),
            lActiveVertices.begin(), lActiveVertices.end());

          lActiveVertices.clear();

          #pragma omp barrier

          #pragma omp single
          {
            this->printMsg("Computed Maxima",
              0.1, // progress
              localTimer.getElapsedTime(), this->threadNumber_);
          }

          // compress paths until no changes occur
          while(globalChange) {
            #pragma omp barrier

            #pragma omp single
            {
              globalChange = false;
              nActiveVertices = activeVertices->size();
            }

            localChange = false;

            #pragma omp single 
            {
              std::string msgdbg = "iteration: " + std::to_string(iterations);
              this->printMsg(msgdbg, 0.12,
                localTimer.getElapsedTime(), this->threadNumber_);
            }

            #pragma omp for schedule(dynamic)
            for(SimplexId i = 0; i < nActiveVertices; i++) {
              SimplexId v = activeVertices->at(i);


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
          
            #pragma omp atomic update
              globalChange |= localChange;
            
            #pragma omp single {
              iterations += 1;
            }
            
            #pragma omp barrier
          }

          #pragma omp single
          {
            this->printMsg("Compressed Paths",
              0.1, // progress
              localTimer.getElapsedTime(), this->threadNumber_);
          }

          #pragma omp for schedule(dynamic, 4)
          for(SimplexId i = 0; i < nVertices; i++) {
            outputData[i] = maxima[outputData[i]];
          }
        }*/
#endif

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
