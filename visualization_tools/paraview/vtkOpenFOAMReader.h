/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOpenFOAMReader - reads a dataset in OpenFOAM format
// .SECTION Description
// vtkOpenFOAMReader creates an multiblock dataset. It reads a controlDict
// file, mesh information, and time dependent data.  The controlDict file
// contains timestep information. The polyMesh folders contain mesh information
// The time folders contain transient data for the cells  Each folder can
// contain any number of data files.

// .SECTION Thanks
// Thanks to Terry Jordan of SAIC at the National Energy
// Technology Laboratory who developed this class.
// Please address all comments to Terry Jordan (terry.jordan@sa.netl.doe.gov)

#ifndef __vtkOpenFOAMReader_h
#define __vtkOpenFOAMReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkstd/vector"
#include "vtkstd/string"
#include "vtkstd/map"

typedef struct
{
  int faceIndex;
  bool neighborFace;
} face;

class vtkUnstructuredGrid;
class vtkPoints;
class vtkIntArray;
class vtkFloatArray;
class vtkDoubleArray;
class vtkDataArraySelection;

class VTK_IO_EXPORT vtkOpenFOAMReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  //METHODS FOR PARAVIEW
  static vtkOpenFOAMReader *New();
  vtkTypeRevisionMacro(vtkOpenFOAMReader, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  //PARAVIEW MACROS
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  vtkGetMacro(NumberOfTimeSteps, int);
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
  vtkGetVector2Macro(TimeStepRange, int);
  vtkSetVector2Macro(TimeStepRange, int);

  //CELL ARRAY METHODS
  int GetNumberOfCellArrays(void);
  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();


protected:
  vtkOpenFOAMReader();
  ~vtkOpenFOAMReader();
  int RequestData(vtkInformation *,
    vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation *,
    vtkInformationVector **, vtkInformationVector *);

private:
  char * FileName;
  int NumberOfTimeSteps;
  int TimeStep;
  int TimeStepRange[2];
  double * Steps;
  bool RequestInformationFlag;
  int StartFace;
  //BTX
  vtkstd::string Path;
  vtkstd::string PathPrefix;
  vtkstd::vector< vtkstd::string > TimeStepData;
  vtkDataArraySelection * CellDataArraySelection;
  vtkstd::vector< vtkstd::vector<int> > FacePoints;
  vtkstd::vector< vtkstd::vector<int> > FacesOwnerCell;
  vtkstd::vector< vtkstd::vector<int> > FacesNeighborCell;
  vtkstd::vector< vtkstd::vector<face> > FacesOfCell;
  vtkPoints * Points;
  vtkIdType NumCells;
  vtkIdType NumFaces;
  vtkIntArray * FaceOwner;
  //vtkIntArray * FaceNeighbor;
  vtkstd::vector< vtkstd::string > PolyMeshPointsDir;
  vtkstd::vector< vtkstd::string > PolyMeshFacesDir;
  vtkIdType NumPoints;
  vtkstd::vector< int > SizeOfBoundary;
  vtkstd::vector< vtkstd::string > BoundaryNames;
  vtkstd::vector< vtkstd::string > PointZoneNames;
  vtkstd::vector< vtkstd::string > FaceZoneNames;
  vtkstd::vector< vtkstd::string > CellZoneNames;
  int NumBlocks;
  void CombineOwnerNeigbor();
  vtkUnstructuredGrid * MakeInternalMesh();
  double ControlDictDataParser(vtkstd::string);
  void ReadControlDict ();
  void GetPoints (int);
  void ReadFacesFile (vtkstd::string);
  void ReadOwnerFile(vtkstd::string);
  void ReadNeighborFile(vtkstd::string);
  void PopulatePolyMeshDirArrays();
  vtkstd::string GetDataType(vtkstd::string, vtkstd::string);
  vtkDoubleArray * GetInternalVariableAtTimestep( vtkstd::string, int);
  vtkDoubleArray * GetBoundaryVariableAtTimestep(int, vtkstd::string, int,
                                                 vtkUnstructuredGrid *);
  vtkstd::vector< vtkstd::string > GatherBlocks(vtkstd::string, int);
  vtkUnstructuredGrid * GetBoundaryMesh(int, int);
  vtkUnstructuredGrid * GetPointZoneMesh(int, int);
  vtkUnstructuredGrid * GetFaceZoneMesh(int, int);
  vtkUnstructuredGrid * GetCellZoneMesh(int, int);
  void CreateDataSet(vtkMultiBlockDataSet *);
  //ETX
};

#endif
