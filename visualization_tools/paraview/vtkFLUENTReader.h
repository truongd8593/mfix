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
// .NAME vtkFLUENTReader - reads a dataset in Fluent file format
// .SECTION Description
// vtkFLUENTReader creates an unstructured grid dataset. It reads .cas and
// .dat files stored in FLUENT native format.
//
// .SECTION Thanks
// Thanks to Brian W. Dotson & Terry E. Jordan (Department of Energy, National
// Energy Technology Laboratory) & Douglas McCorkle (Iowa State University)
// who developed this class.
// Please address all comments to Brian Dotson (brian.dotson@netl.doe.gov) &
// Terry Jordan (terry.jordan@sa.netl.doe.gov)
//
// Please address all comments to Brian Dotson (Brian.Dotson@netl.doe.gov)

// .SECTION See Also
// vtkGAMBITReader

#ifndef __vtkFLUENTReader_h
#define __vtkFLUENTReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"

#include "vtkstd/map"
#include "vtkstd/vector"
#include "vtkstd/set"
#include <fstream>
#include <sstream>

class vtkDataArraySelection;
class vtkPoints;
class vtkTriangle;
class vtkTetra;
class vtkQuad;
class vtkHexahedron;
class vtkPyramid;
class vtkWedge;
class vtkConvexPointSet;

class VTK_IO_EXPORT vtkFLUENTReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkFLUENTReader *New();
  vtkTypeRevisionMacro(vtkFLUENTReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the file name of the Fluent case file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the total number of cells. The number of cells is only valid after a
  // successful read of the data file is performed.
  vtkGetMacro(NumberOfCells,int);
  //vtkGetMacro(NumberOfCellArrays, int);

  int GetNumberOfCellArrays(void);
  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();

protected:
  vtkFLUENTReader();
  ~vtkFLUENTReader();
  int RequestInformation(vtkInformation *, 
    vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, 
    vtkInformationVector *);
  vtkDataArraySelection* CellDataArraySelection;
  char * FileName;
  int NumberOfCells;
  int NumberOfCellArrays;
  int                    OpenCaseFile(const char *filename);
  int                    OpenDataFile(const char *filename);
  int                    GetCaseChunk ();
  void                   GetNumberOfCellZones();
  int                    GetCaseIndex();
  void                   LoadVariableNames();
  int                    GetDataIndex();
  int                    GetDataChunk();

  void                   ParseCaseFile();
  int                    GetDimension();
  void                   GetLittleEndianFlag();
  void                   GetNodesAscii();
  void                   GetNodesSinglePrecision();
  void                   GetNodesDoublePrecision();
  void                   GetCellsAscii();
  void                   GetCellsBinary();
  void                   GetFacesAscii();
  void                   GetFacesBinary();
  void                   GetPeriodicShadowFacesAscii();
  void                   GetPeriodicShadowFacesBinary();
  void                   GetCellTreeAscii();
  void                   GetCellTreeBinary();
  void                   GetFaceTreeAscii();
  void                   GetFaceTreeBinary();
  void                   GetInterfaceFaceParentsAscii();
  void                   GetInterfaceFaceParentsBinary();
  void                   GetNonconformalGridInterfaceFaceInformationAscii();
  void                   GetNonconformalGridInterfaceFaceInformationBinary();
  void                   CleanCells();
  void                   PopulateCellNodes();
  int                    GetCaseBufferInt(int ptr);
  float                  GetCaseBufferFloat(int ptr);
  double                 GetCaseBufferDouble(int ptr);
  void                   PopulateTriangleCell(int i);
  void                   PopulateTetraCell(int i);
  void                   PopulateQuadCell(int i);
  void                   PopulateHexahedronCell(int i);
  void                   PopulatePyramidCell(int i);
  void                   PopulateWedgeCell(int i);
  void                   PopulatePolyhedronCell(int i);
  void                   ParseDataFile();
  int                    GetDataBufferInt(int ptr);
  float                  GetDataBufferFloat(int ptr);
  double                 GetDataBufferDouble(int ptr);
  void                   GetData(int dataType);

  //
  //  Variables
  //
  //BTX
  ifstream FluentCaseFile;
  ifstream FluentDataFile;
  vtkstd::string CaseBuffer;
  vtkstd::string DataBuffer;

  vtkPoints           *Points;
  vtkTriangle         *Triangle;
  vtkTetra            *Tetra;
  vtkQuad             *Quad;
  vtkHexahedron       *Hexahedron;
  vtkPyramid          *Pyramid;
  vtkWedge            *Wedge;
  vtkConvexPointSet   *ConvexPointSet;

  struct Cell 
    {
    int type;
    int zone;
    vtkstd::vector<int> faces;
    int parent;
    int child;
    vtkstd::vector<int> nodes;
    };

  struct Face
    {
    int type;
    int zone;
    vtkstd::vector<int> nodes;
    int c0;
    int c1;
    int periodicShadow;
    int parent;
    int child;
    int interfaceFaceParent;
    int interfaceFaceChild;
    int ncgParent;
    int ncgChild;
    };

  struct ScalarDataChunk
    {
    int subsectionId;
    int zoneId;
    vtkstd::vector<double> scalarData;
    };

  struct VectorDataChunk
    {
    int subsectionId;
    int zoneId;
    vtkstd::vector<double> iComponentData;
    vtkstd::vector<double> jComponentData;
    vtkstd::vector<double> kComponentData;
    };


  vtkstd::vector< Cell > Cells;
  vtkstd::vector< Face > Faces;
  vtkstd::map< int, vtkstd::string > VariableNames;
  vtkstd::vector< int >  CellZones;
  vtkstd::vector< ScalarDataChunk > ScalarDataChunks;
  vtkstd::vector< VectorDataChunk > VectorDataChunks;

  vtkstd::vector< vtkstd::vector<int> > SubSectionZones;
  vtkstd::vector< int > SubSectionIds;
  vtkstd::vector< int > SubSectionSize;

  vtkstd::vector< vtkstd::string > ScalarVariableNames;
  vtkstd::vector< int > ScalarSubSectionIds;
  vtkstd::vector< vtkstd::string > VectorVariableNames;
  vtkstd::vector< int > VectorSubSectionIds;

  int LittleEndianFlag;
  int GridDimension;
  int DataPass;
  int NumberOfScalars;
  int NumberOfVectors;
  //ETX
};
#endif
