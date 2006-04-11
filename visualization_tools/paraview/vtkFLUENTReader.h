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
// Thanks to Brian W. Dotson (Department of Energy, National Energy 
//       Technology Laboratory) who developed this class.
//
// Please address all comments to Brian Dotson (Brian.Dotson@netl.doe.gov)

// .SECTION See Also
// vtkGAMBITReader

#ifndef __vtkFLUENTReader_h
#define __vtkFLUENTReader_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkIntArray;
class vtkFloatArray;
class vtkIdTypeArray;
class vtkDataArraySelection;
class vtkPoints;
class vtkDoubleArray;
class vtkTriangle;
class vtkQuad;
class vtkTetra;
class vtkPyramid;
class vtkWedge;
class vtkHexahedron;
class VTK_IO_EXPORT vtkFLUENTReader : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkFLUENTReader *New();
  vtkTypeRevisionMacro(vtkFLUENTReader,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the file name of the Fluent case file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the total number of cells. The number of cells is only valid after a
  // successful read of the data file is performed.
  vtkGetMacro(NumberOfCells,int);
  vtkGetMacro(NumberOfCellArrays, int);

  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name);

protected:
  vtkFLUENTReader();
  ~vtkFLUENTReader();
  int RequestInformation(vtkInformation *, 
    vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, 
    vtkInformationVector *);
  int BinaryFile;
  int NumberOfCells;
  int NumberOfCellComponents;
  int NumberOfCellFields;
  int NumberOfCellArrays;
  ifstream *FileStream;
  ifstream *DataFileStream;
  int ObjectsFlag;

  vtkDataArraySelection* CellDataArraySelection;
  //BTX
  struct DataInfo
    {
    long FileOffset; // offset in binary file
    int VectorLength;  // number of components in the node or cell variable
    float Minimum[3]; // pre-calculated data minima (max size 3 for vectors)
    float Maximum[3]; // pre-calculated data maxima (max size 3 for vectors)
    };
  //ETX
  DataInfo *CellDataInfo;

private:

  void ReadFile(vtkUnstructuredGrid *output);
  void CreateVTKObjects(void);
  void DeleteVTKObjects(void);
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();
  void GetCellDataRange(int cellComp, int index, float *min, float *max);
  int OpenCaseAndDataFiles(void);
  void ParseCaseFile(void);
  void MakeFaceTreeParentTable(void);
  void MakeCellTreeParentTable(void);
  void LoadFaceKidFlags(void);
  void LoadFaceParentFlags(void);
  void LoadInterfaceFaceChildFlags(void);
  void LoadNCGFaceChildFlags(void);
  void BuildCells(void);
  void LoadCellParentFlags(void);
  void LoadCellNumberOfFaces(void);
  void LoadCellFaces(void);
  void RemoveExtraFaces(void);
  void ParseDataFile(void);
  void InitializeVariableNames(void);
  int GetCaseIndex(int ix);
  int ExecuteCaseTask(int task, int file_index);
  int GetDataIndex(int ix);
  int ExecuteDataTask(int task, int file_index);
  int GetComment(int ix);
  int GetHeader(int ix);
  int GetGridDimension(int ix);
  int GetMachineConfiguration(int ix);
  int GetVariablesASCII(int ix);
  int GetCortexVariablesASCII(int ix);
  int GetDomainVariablesASCII(int ix);
  int GetCellsASCII(int ix);
  int GetFacesASCII(int ix);
  int GetNodesASCII(int ix);
  int GetFaceParentsASCII(int ix);
  int GetNCG1InformationASCII(int ix);
  int GetNCG2InformationASCII(int ix);
  int GetNodeFlagsASCII(int ix);
  int Command54(int ix);
  int GetZoneSectionsASCII(int ix);
  int GetPeriodicShadowFacesASCII(int ix);
  int GetGridSizeASCII(int ix);
  int GetPartitionASCII(int ix);
  int GetCellTreeASCII(int ix);
  int GetFaceTreeASCII(int ix);

  int GetVariablesSinglePrecision(int ix);
  int GetCortexVariablesSinglePrecision(int ix);
  int GetDomainVariablesSinglePrecision(int ix);
  int GetCellsSinglePrecision(int ix);
  int GetFacesSinglePrecision(int ix );
  int GetNodesSinglePrecision(int ix );
  int GetFaceParentsSinglePrecision(int ix );
  int GetNCG1InformationSinglePrecision(int ix );
  int GetNCG2InformationSinglePrecision(int ix );
  int GetNodeFlagsSinglePrecision(int ix );
  int GetZoneSectionsSinglePrecision(int ix );
  int GetPeriodicShadowFacesSinglePrecision(int ix );
  int GetGridSizeSinglePrecision(int ix );
  int GetPartitionSinglePrecision(int ix );
  int GetCellTreeSinglePrecision(int ix );
  int GetFaceTreeSinglePrecision(int ix );
  int GetVariablesDoublePrecision(int ix );
  int GetCortexVariablesDoublePrecision(int ix );
  int GetDomainVariablesDoublePrecision(int ix );
  int GetCellsDoublePrecision(int ix );
  int GetFacesDoublePrecision(int ix );
  int GetNodesDoublePrecision(int ix );
  int GetFaceParentsDoublePrecision(int ix );
  int GetNCG1InformationDoublePrecision(int ix );
  int GetNCG2InformationDoublePrecision(int ix );
  int GetNodeFlagsDoublePrecision(int ix );
  int GetZoneSectionsDoublePrecision(int ix );
  int GetPeriodicShadowFacesDoublePrecision(int ix );
  int GetGridSizeDoublePrecision(int ix );
  int GetPartitionDoublePrecision(int ix );
  int GetCellTreeDoublePrecision(int ix );
  int GetFaceTreeDoublePrecision(int ix );
  int GetDataComment(int ix );
  int GetDataHeader(int ix );
  int GetDataGridDimension(int ix );
  int GetDataMachineConfiguration(int ix );
  int GetDataGridSizeASCII(int ix );
  int GetDataVariablesASCII(int ix );
  int GetUnknownASCII313(int ix );
  int GoToNextRightParenData(int ix );
  int GoToNextLeftParenData(int ix );
  int GoToNextSectionASCIIData(int ix );
  int GoToNextSectionSinglePrecisionData(int ix, char buf[] );
  int GoToNextSectionDoublePrecisionData(int ix, char buf[] );
  int GetDataASCII(int ix );
  int GetDataSinglePrecision(int ix );
  int GetDataDoublePrecision(int ix );
  int GetUnknownSinglePrecision2301(int ix );
  int GetUnknownSinglePrecision2302(int ix );
  int GetUnknownSinglePrecision2313(int ix );
  int GetUnknownDoublePrecision3301(int ix );
  int GetUnknownDoublePrecision3302(int ix );
  int GetUnknownDoublePrecision3313(int ix );
  int GetUnknownASCII301(int ix );
  int GetUnknownASCII302(int ix );
  int GetUnknownASCII303(int ix );
  int GetDataUnknownASCII(int ix );
  void GetStringToNextRightParenData(int ix, char buf[]  );
  int IsCellZoneId(int zi );
  int IsNewVariable(int ssid );
  int GetVariableIndex(int ssid );
  int GetBinaryIntegerData(int ix );
  float GetBinaryFloatData(int ix );
  double GetBinaryDoubleData(int ix );
  int IsASCIICharacterHexDigit(int ix );
  int GoToNextASCIIHexDigit(int ix );
  int GoToNextRightParen(int ix );
  int GoToNextLeftParen(int ix );
  int GoToNextEOL(int ix );
  int GoToNextSectionASCII(int ix );
  void GetStringToNextRightParen(int ix, char buf[]  );
  void GetStringToNextRightParenOrEOL(int ix, char buf[]  );
  void GetMixedCellTypes(int ix, int fi, int li );
  int GoToNextSectionSinglePrecision(int ix, char buf[] );
  int GoToNextSectionDoublePrecision(int ix, char buf[] );
  int GetBinaryInteger(int ix );
  float GetBinaryFloat(int ix );
  double GetBinaryDouble(int ix );
  int GetAsciiInteger(int ix);
  int GoPastAsciiInteger(int ix);
  int GoToNextEOLData(int ix);
  void GetStringToNextRightParenOrEOLData(int ix, char buf[] );

  char *FileName;
  char *DataFileName;
  char *CaseFileBuffer;
  char *DataFileBuffer;
  int CaseFileBufferLength;
  int DataFileBufferLength;
  int GridDimension;
  int NumberOfNodes;
  int NumberOfFaces;
  int NumberOfFaceParents;
  int NumberOfPeriodicShadowFaces;
  int NumberOfCellZones;
  int NumberOfVariables;
  int LittleEndianFlag;
  int NumberOfFaceTrees;
  int NumberOfFaceTreeKids;
  int NumberOfFaceTreeParents;
  int LastFaceTreeParent;
  int NumberOfCellTrees;
  int NumberOfCellTreeKids;
  int NumberOfCellTreeParents;
  int LastCellTreeParent;
  int NumberOfNCGFaceHeaders;
  int NumberOfNCGFaces;
  int NumberOfNCGNodeHeaders;
  int NumberOfNCGNodes;
  int DataPass;
  int NumberOfFaceParentChildren;

  int *VectorLength;
  float *Minimum;
  float *Maximum;

  vtkPoints *Points;
  vtkIntArray *CellTypes;
  vtkIntArray *CellFaces;
  vtkIntArray *CellFacesClean;
  vtkIntArray *FaceTypes;
  vtkIntArray *FaceNodes;
  vtkIntArray *FaceCells;
  vtkIntArray *FaceParents;
  vtkIntArray *PeriodicShadowFaces;
  vtkIntArray *FaceTreesNumberOfKids;
  vtkIntArray *FaceTreesKids;
  vtkIntArray *FaceTreesKidsIndex;
  vtkIntArray *CellTreesNumberOfKids;
  vtkIntArray *CellTreesKids;
  vtkIntArray *CellTreesKidsIndex;
  vtkIntArray *FaceTreeParentFaceId0;
  vtkIntArray *FaceTreeParentFaceId1;
  vtkIntArray *FaceTreeParentZoneId;
  vtkIntArray *FaceTreeChildZoneId;
  vtkIntArray *FaceTreeParentTable;
  vtkIntArray *CellTreeParentCellId0;
  vtkIntArray *CellTreeParentCellId1;
  vtkIntArray *CellTreeParentZoneId;
  vtkIntArray *CellTreeChildZoneId;
  vtkIntArray *CellTreeParentTable;
  vtkIntArray *NCGFaceKidId;
  vtkIntArray *NCGFaceParentId;
  vtkIntArray *NCGFaceNumberOfFaces;
  vtkIntArray *NCGFaceChild;
  vtkIntArray *NCGFaceParent;
  vtkIntArray *NCGNodeZoneId;
  vtkIntArray *NCGNodeNumberOfNodesNCG;
  vtkDoubleArray *NCGNodes;
  vtkIntArray *NCGNodeIds;
  vtkIntArray *CellNumberOfFaces;
  vtkIntArray *FaceKidFlags;
  vtkIntArray *FaceParentFlags;
  vtkIntArray *CellIndex;
  vtkIntArray *InterfaceFaceChildFlags;
  vtkIntArray *FaceParentsChildren;
  vtkIntArray *NCGFaceChildFlags;
  vtkIntArray *CellParentFlags;
  vtkTriangle *ATriangle;
  vtkQuad *AQuad;
  vtkTetra *ATetra;
  vtkPyramid *APyramid;
  vtkWedge *AWedge;
  vtkHexahedron *AHexahedron;
  vtkIntArray *CellZones;
  vtkIntArray *VariableIds;
  vtkIntArray *VariableSizes;
  vtkDoubleArray **CellData;
  char *VariableNames[1500];
  vtkUnstructuredGrid *Mesh;

  vtkFLUENTReader(const vtkFLUENTReader&);  // Not implemented.
  void operator=(const vtkFLUENTReader&);  // Not implemented.
};
#endif
